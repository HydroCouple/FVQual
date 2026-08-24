// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   internalmode.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  InternalModeSolver — the 3-D hydrostatic baroclinic mode and the
 *         mode-splitting driver.
 *
 * This is the heart of the model. It advances the three-dimensional velocity,
 * the layer geometry, and the transported scalars, and it drives the external
 * mode as a subcycled inner loop.
 *
 * **Why splitting.** Surface gravity waves are fast and, for a reservoir
 * water-quality question, uninteresting: they set a time step of order a
 * second while the physics that matters — stratification, baroclinic
 * circulation, transport — evolves over minutes to seasons. Mode splitting
 * confines the fast wave to the cheap two-dimensional external solver and
 * lets the expensive three-dimensional state advance at the internal-wave and
 * advective limit, typically tens of seconds. Nothing here needs a linear
 * solve except the vertical column sweeps, which are tridiagonal.
 *
 * **The consistency problem splitting creates.** The two modes can disagree:
 * the depth-integrated internal transport need not equal what the external
 * mode actually moved, and if a scalar is transported on inconsistent fluxes,
 * a uniform tracer stops being uniform and mass stops being conserved. Two
 * mechanisms keep them honest:
 *
 * 1. *Mode adjustment.* After the internal momentum update, the depth-mean of
 *    the three-dimensional velocity is replaced by the external mode's
 *    depth-averaged velocity. The internal mode supplies the shear; the
 *    external mode owns the transport.
 * 2. *Flux adjustment on the time-averaged transport.* Layer fluxes are
 *    rescaled so that their sum across each face equals the volume the
 *    external mode actually moved through it, averaged over the subcycling
 *    interval.
 *
 * Together these make the vertical velocity diagnosed from continuity vanish
 * at the free surface to round-off — which is exactly the residual reported as
 * `surfaceFluxResidual`, and exactly what the verification suite checks.
 *
 * **Order of operations in a step**, each stage depending on the last:
 * density from the scalars, mixing coefficients from the resulting
 * stratification, baroclinic pressure from the density field, the external
 * mode, internal momentum, the new layer geometry, mode adjustment,
 * mass-consistent fluxes, vertical velocity, and finally scalar transport on
 * those fluxes.
 *
 * \threadsafety One solver per thread; step() mutates the mesh and state it
 * was initialized with.
 */

#ifndef FVQual_HYDRO_INTERNALMODE_H
#define FVQual_HYDRO_INTERNALMODE_H

#include "fvqual/fvqual_export.h"
#include "fvqual/hydro/barotropicsolver.h"
#include "fvqual/hydro/pressureclosure.h"
#include "fvqual/kinetics/reactionnetwork.h"
#include "fvqual/forcing/sources.h"
#include "fvqual/forcing/outletgroups.h"
#include "fvqual/forcing/structures.h"
#include "fvqual/forcing/transfers.h"
#include "fvqual/transport/euleriantransport.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/state/statestore.h"

#include <memory>
#include <span>
#include <string>
#include <vector>

namespace FVQual::Hydro {

/*!
 * \brief Physical and numerical settings of the internal mode. SI units.
 */
struct FVQual_EXPORT InternalModeOptions
{
    double gravity          = 9.80665; //!< Gravitational acceleration (m/s²).
    double referenceDensity = 1000.0;  //!< Boussinesq reference density (kg/m³).

    /*!
     * \brief Coriolis parameter (1/s). Zero for reservoirs small enough to
     *        ignore rotation; set to 2 Ω sin(latitude) for large lakes.
     */
    double coriolis = 0.0;

    // ── Vertical mixing (Pacanowski-Philander form) ───────────────────────

    /*!
     * \brief Neutral (unstratified) eddy viscosity (m²/s).
     */
    double neutralViscosity = 1.0e-2;

    /*!
     * \brief Background viscosity and diffusivity retained even in strong
     *        stratification (m²/s). The diffusivity floor is deliberately
     *        small: an over-generous floor erodes the thermocline, which is
     *        the classic way a reservoir model loses its stratification.
     */
    double backgroundViscosity   = 1.0e-6;
    double backgroundDiffusivity = 1.0e-7;

    /*!
     * \brief Mixing applied where the column is statically unstable (m²/s) —
     *        convective adjustment.
     */
    double convectiveMixing = 1.0e-1;

    // ── Horizontal mixing ─────────────────────────────────────────────────

    double horizontalViscosity   = 0.0; //!< Constant eddy viscosity (m²/s).
    double horizontalDiffusivity = 0.0; //!< Constant eddy diffusivity (m²/s).

    /*!
     * \brief Advective discretization for transported scalars.
     *
     * MUSCL with flux correction by default. First-order upwind remains
     * selectable because it is the honest baseline the accuracy tests measure
     * against — and because a bounded, dependable scheme is worth having when
     * something else is under suspicion.
     */
    Transport::AdvectionScheme advectionScheme = Transport::AdvectionScheme::MusclFct;

    // ── Forcing ───────────────────────────────────────────────────────────

    double windStressX = 0.0; //!< Surface wind stress, x (Pa).
    double windStressY = 0.0; //!< Surface wind stress, y (Pa).

    /*!
     * \brief Net surface heat flux into the water (W/m²), positive warming.
     *
     * A prescribed net flux, not a bulk meteorological formulation: the
     * latent/sensible/radiative decomposition and its met-data plumbing
     * arrive with the heat-transport work in a later phase.
     */
    double surfaceHeatFlux = 0.0;

    /*!
     * \brief Specific heat capacity of water (J/(kg·K)).
     */
    double heatCapacity = 4186.0;

    double manningN = 0.03; //!< Bed roughness for the bottom stress.

    // ── Time stepping ─────────────────────────────────────────────────────

    /*!
     * \brief Courant number for the internal step.
     */
    double cflNumber = 0.5;

    double minTimeStep = 1.0e-3; //!< Floor on the internal step (s).
    double maxTimeStep = 300.0;  //!< Ceiling on the internal step (s).

    /*!
     * \brief Salinity species id used by the equation of state; empty to
     *        treat the water as fresh.
     */
    std::string salinitySpeciesId;
};

/*!
 * \brief Per-step diagnostics of the internal mode.
 */
/*!
 * \brief Wall-clock time spent in each phase of a step, in nanoseconds.
 *
 * Accumulated across steps and reset on demand. The cost is a handful of clock
 * reads per internal step, which is nothing beside the per-cell work they
 * bracket, and it is the only way to decide *what* to accelerate: a kernel that
 * looks expensive on the page and takes 2% of the run is not where the time is.
 * Phase 7 is meant to be driven by this rather than by intuition.
 */
struct FVQual_EXPORT StepProfile
{
    int64_t density    = 0; //!< Equation of state and vertical mixing.
    int64_t pressure   = 0; //!< Baroclinic pressure gradient.
    int64_t forcing    = 0; //!< Transfers, outlet groups, structures, sources.
    int64_t external   = 0; //!< The barotropic subcycle.
    int64_t momentum   = 0; //!< Internal momentum advance.
    int64_t geometry   = 0; //!< Surface update and mode reconciliation.
    int64_t fluxes     = 0; //!< Layer and vertical flux construction.
    int64_t transport  = 0; //!< Scalar advection, diffusion, heat flux.
    int64_t reactions  = 0; //!< Operator-split kinetics.
    int64_t total      = 0; //!< The whole step.

    void reset() noexcept { *this = StepProfile{}; }
};

struct FVQual_EXPORT InternalModeDiagnostics
{
    double timeStep              = 0.0; //!< Internal step taken (s).
    double maxHorizontalVelocity = 0.0; //!< Largest |u| over the cells (m/s).
    double maxVerticalFlux       = 0.0; //!< Largest interface flux (m/s).

    /*!
     * \brief Largest free-surface continuity residual (m/s).
     *
     * The vertical velocity is integrated up from a no-flux bed, so its value
     * at the free surface must vanish if the internal fluxes are consistent
     * with the external mode. A non-negligible value here means the modes
     * have drifted apart and scalar transport can no longer be trusted.
     */
    double surfaceFluxResidual = 0.0;

    double totalVolume    = 0.0; //!< Domain water volume (m³).
    double maxRichardson  = 0.0; //!< Largest gradient Richardson number.
    int64_t externalSteps = 0;   //!< External substeps taken.
    int64_t steps         = 0;   //!< Internal steps taken by advance().
};

/*!
 * \brief Three-dimensional hydrostatic baroclinic solver.
 */
class FVQual_EXPORT InternalModeSolver
{
public:
    InternalModeSolver();
    ~InternalModeSolver();

    InternalModeSolver(const InternalModeSolver &)            = delete;
    InternalModeSolver &operator=(const InternalModeSolver &) = delete;

    /*!
     * \brief Binds the solver to a mesh, a state store, and a species set.
     *
     * The mesh and state must outlive the solver; both are mutated by step().
     * Every column must be wet: wetting and drying of the three-dimensional
     * mode is deferred (F3.4).
     *
     * \param[in,out] mesh    Layered mesh; its layer geometry is updated each
     *                        step as the free surface moves.
     * \param[in,out] state   State store, already allocated for \p mesh.
     * \param[in]     species Species registry the state was allocated with.
     * \param[in]     options Physical and numerical settings.
     * \param[out]    message Failure description when returning false.
     */
    [[nodiscard]] bool initialize(Mesh::LayeredMesh &mesh, State::StateStore &state,
                                  const State::SpeciesRegistry &species,
                                  const InternalModeOptions &options,
                                  std::string &message);

    /*!
     * \brief Replaces the pressure closure.
     *
     * Defaults to HydrostaticPressureClosure. Passing a non-hydrostatic
     * closure is how the correction is introduced later without touching the
     * momentum assembly.
     */
    void setPressureClosure(std::unique_ptr<IPressureClosure> closure);

    [[nodiscard]] const IPressureClosure &pressureClosure() const noexcept
    {
        return *m_pressure;
    }

    /*!
     * \brief Largest stable internal time step (s), from the advective and
     *        internal-wave Courant limits.
     */
    [[nodiscard]] double stableTimeStep() const;

    /*!
     * \brief Advances the full three-dimensional state by \p dt seconds,
     *        subcycling the external mode within it.
     */
    [[nodiscard]] bool step(double dt, std::string &message);

    /*!
     * \brief Advances by \p duration seconds, substepping at the internal
     *        Courant limit and landing exactly on the target time.
     */
    [[nodiscard]] bool advance(double duration, std::string &message);

    // ── Access ────────────────────────────────────────────────────────────

    [[nodiscard]] BarotropicSolver       &external() noexcept { return m_external; }
    [[nodiscard]] const BarotropicSolver &external() const noexcept { return m_external; }

    [[nodiscard]] double time() const noexcept { return m_time; }

    [[nodiscard]] const InternalModeDiagnostics &diagnostics() const noexcept
    {
        return m_diagnostics;
    }

    /*!
     * \brief Upward volume flux through each layer interface, per unit plan
     *        area (m/s) — the grid-relative vertical velocity that transport
     *        uses. Indexed by Mesh::LayeredMesh::interfaceSlot().
     */
    [[nodiscard]] const std::vector<double> &verticalFlux() const noexcept
    {
        return m_verticalFlux;
    }

    /*!
     * \brief Vertical eddy viscosity and diffusivity at interfaces (m²/s).
     */
    [[nodiscard]] const std::vector<double> &eddyViscosity() const noexcept
    {
        return m_km;
    }
    [[nodiscard]] const std::vector<double> &eddyDiffusivity() const noexcept
    {
        return m_kz;
    }

    /*!
     * \brief Inflows and withdrawals acting on this model.
     *
     * Configure before or during a run; placement is recomputed every step,
     * so a plunging inflow follows the thermocline and a withdrawal envelope
     * responds to changing stratification.
     */
    [[nodiscard]] Forcing::SourceSet       &sources() noexcept { return m_sources; }
    [[nodiscard]] const Forcing::SourceSet &sources() const noexcept { return m_sources; }

    /*!
     * \brief Outlet works whose discharge follows from the water level.
     *
     * Evaluated once per step, before the sources are placed, so a rating
     * curve always sees the level the step begins from. Register a structure
     * with `structures().add(spec, sources(), message)` — the pairing between
     * a structure and the source it drives is made there so it cannot come
     * apart.
     */
    /*!
     * \brief Accumulated per-phase timings (ns).
     */
    [[nodiscard]] const StepProfile &profile() const noexcept { return m_profile; }

    void resetProfile() noexcept { m_profile.reset(); }

    [[nodiscard]] Forcing::TransferSet &transfers() noexcept { return m_transfers; }
    [[nodiscard]] const Forcing::TransferSet &transfers() const noexcept
    {
        return m_transfers;
    }

    [[nodiscard]] Forcing::OutletGroupSet &outletGroups() noexcept { return m_groups; }
    [[nodiscard]] const Forcing::OutletGroupSet &outletGroups() const noexcept
    {
        return m_groups;
    }

    [[nodiscard]] Forcing::StructureSet &structures() noexcept { return m_structures; }
    [[nodiscard]] const Forcing::StructureSet &structures() const noexcept
    {
        return m_structures;
    }

    /*!
     * \brief Rate laws applied to the transported scalars each step.
     *
     * Configure the network here, then call compileReactions() before
     * stepping. Empty by default, so hydrodynamics-only runs pay nothing.
     */
    [[nodiscard]] Kinetics::ReactionNetwork &reactions() noexcept { return m_reactions; }
    [[nodiscard]] const Kinetics::ReactionNetwork &reactions() const noexcept
    {
        return m_reactions;
    }

    /*!
     * \brief Resolves the reaction network's names against this solver's species.
     *
     * Exists rather than exposing the registry so the network cannot be bound
     * against a *different* registry than the state was allocated with — the
     * indices would resolve and the wrong species would react.
     *
     * \param[out] message Failure description when returning false.
     */
    [[nodiscard]] bool compileReactions(std::string &message);

    /*!
     * \brief Recomputes cell density from the current scalars.
     */
    void updateDensity();

    /*!
     * \brief Volume-weighted total of a transported species (concentration x m³).
     */
    [[nodiscard]] double totalSpeciesMass(int species) const;

private:
    void   updateVerticalMixing();
    void   advanceMomentum(double dt);
    void   applyModeAdjustment();
    void   computeLayerFluxes();
    void   diagnoseVerticalFlux(double dt);
    void   transportScalars(double dt);
    void   applyHorizontalDiffusion(double dt, std::span<double> values);
    void   applySurfaceHeatFlux(double dt);
    double bottomDragCoefficient(int64_t column) const;

    Mesh::LayeredMesh            *m_mesh    = nullptr;
    State::StateStore            *m_state   = nullptr;
    const State::SpeciesRegistry *m_species = nullptr;

    InternalModeOptions               m_options;
    InternalModeDiagnostics           m_diagnostics;
    BarotropicSolver                  m_external;
    Transport::EulerianTransport      m_transport;
    Forcing::SourceSet                m_sources;
    Forcing::StructureSet             m_structures;
    Forcing::OutletGroupSet           m_groups;
    Forcing::TransferSet              m_transfers;
    Kinetics::ReactionNetwork         m_reactions;
    StepProfile                       m_profile;
    std::unique_ptr<IPressureClosure> m_pressure;

    double m_time             = 0.0;
    int    m_temperatureIndex = -1;
    int    m_salinityIndex    = -1;

    // Working arrays, allocated once.
    std::vector<double> m_accelerationX;
    std::vector<double> m_accelerationY;
    std::vector<double> m_columnAccelerationX;
    std::vector<double> m_columnAccelerationY;
    std::vector<double> m_km;
    std::vector<double> m_kz;
    std::vector<double> m_edgeLayerFlux;
    std::vector<double> m_verticalFlux;
    std::vector<double> m_oldVolume;
    std::vector<double> m_workU;
    std::vector<double> m_workV;
    std::vector<double> m_scalarWork;
    std::vector<double> m_cellSource;
    std::vector<double> m_columnSource;

    /*!
     * \brief Volume each structure released over the current internal step
     *        (m³), accumulated across the barotropic subcycle.
     */
    std::vector<double> m_structureVolume;
    std::vector<double> m_sourceValue;

    // Tridiagonal scratch, sized to the layer count.
    std::vector<double> m_lower;
    std::vector<double> m_diagonal;
    std::vector<double> m_upper;
    std::vector<double> m_rhs;
    std::vector<double> m_solution;
    std::vector<double> m_sweep;
};

} // namespace FVQual::Hydro

#endif // FVQual_HYDRO_INTERNALMODE_H
