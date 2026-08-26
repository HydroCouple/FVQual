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
#include "fvqual/hydro/heatbudget.h"
#include "fvqual/hydro/turbulence.h"
#include "fvqual/hydro/pressureclosure.h"
#include "fvqual/kinetics/lightfield.h"
#include "fvqual/kinetics/reactionnetwork.h"
#include "fvqual/forcing/forcingfield.h"
#include "fvqual/forcing/sources.h"
#include "fvqual/forcing/outletgroups.h"
#include "fvqual/forcing/structures.h"
#include "fvqual/forcing/transfers.h"
#include "fvqual/transport/euleriantransport.h"
#include "fvqual/transport/settling.h"
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

    /*!
     * \brief Use the two-equation k-epsilon closure instead of Richardson
     *        damping (F5.4).
     *
     * Off by default: Richardson damping is what every gate through Phase 3 was
     * written against, and a closure that changes every mixing coefficient in
     * the model should be opted into and measured, not switched on under
     * existing results.
     *
     * Measured on DeGray, Richardson damping has no setting between off and
     * overwhelming — see `turbulence.h`. This is the fix for that, and it is
     * only worth defaulting on once it beats both.
     */
    bool useTurbulenceClosure = false;

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
    int64_t settling   = 0; //!< Particulate settling and the benthic handoff.
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

    /*!
     * \brief Which constraint set the last stableTimeStep(), and where.
     *
     * `limitingColumn` is -1 when the horizontal constraint bound. A collapsing
     * step is otherwise diagnosed by guessing a column and tracing it, which is
     * slow and frequently wrong — the solver already knows.
     */
    int64_t limitingColumn  = -1;
    int     limitingLayer   = -1;
    double  limitingSpeed   = 0.0;
    bool    limitIsVertical = false;
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
     * \brief Meteorological forcing driving the surface heat budget.
     *
     * Inert until enableMeteorology() is called, so an existing model keeps its
     * prescribed `InternalModeOptions::surfaceHeatFlux` and nothing changes
     * under it.
     *
     * Each variable is uniform or per column (F3.7). Assigning a `Meteorology`
     * makes everything uniform, so `solver.meteorology() = weather;` still does
     * what it always did.
     */
    [[nodiscard]] Forcing::MeteorologyField &meteorology() noexcept
    {
        return m_meteorology;
    }
    [[nodiscard]] const Forcing::MeteorologyField &meteorology() const noexcept
    {
        return m_meteorology;
    }

    /*!
     * \brief Switches the surface heat flux from a prescribed value to the
     *        meteorological budget.
     *
     * \param[in] attenuation Attenuation of sunlight by clear water (1/m),
     *                        governing the depth over which the penetrating
     *                        shortwave is absorbed — and so, more than any other
     *                        single number here, how the model stratifies.
     */
    void enableMeteorology(double attenuation = 0.45) noexcept
    {
        m_useMeteorology   = true;
        m_solarAttenuation = attenuation;
    }

    void disableMeteorology() noexcept { m_useMeteorology = false; }

    [[nodiscard]] bool meteorologyEnabled() const noexcept { return m_useMeteorology; }

    /*!
     * \brief The solar field the heat budget last marched.
     *
     * Separate from light(), which serves the kinetics: photosynthesis responds
     * to the photosynthetically active band and is shaded by algae, while
     * heating follows the whole spectrum through clear water. Sharing one field
     * would force one attenuation to stand for both.
     */
    [[nodiscard]] const Kinetics::LightField &solarField() const noexcept
    {
        return m_solar;
    }

    /*!
     * \brief Per-species settling velocities and their benthic destinations.
     *
     * Applied after transport and before the reactions, so a benthic rate law
     * acts on material that has already arrived this step rather than on last
     * step's sediment.
     */
    [[nodiscard]] Transport::SettlingSet &settling() noexcept { return m_settling; }
    [[nodiscard]] const Transport::SettlingSet &settling() const noexcept
    {
        return m_settling;
    }

    /*!
     * \brief Irradiance down each column, and its optical properties.
     *
     * Recomputed each step from the transported field, immediately before the
     * reactions, so self-shading reflects the constituents the reactions are
     * about to see. Only computed when a rate law references `__LIGHT__` —
     * otherwise the march is work with no consumer.
     */
    [[nodiscard]] Kinetics::LightField &light() noexcept { return m_light; }
    [[nodiscard]] const Kinetics::LightField &light() const noexcept { return m_light; }

    /*!
     * \brief Incident shortwave and optical properties (runtime-settable).
     */
    [[nodiscard]] Kinetics::LightOptions &lightOptions() noexcept
    {
        return m_lightOptions;
    }
    [[nodiscard]] const Kinetics::LightOptions &lightOptions() const noexcept
    {
        return m_lightOptions;
    }

    /*!
     * \brief Binds the external mode's open boundaries into the 3-D transport.
     *
     * Call once, after the boundaries have been added through external(), and
     * before stepping. Each (boundary, adjacent column) pair becomes a
     * Forcing::Source whose discharge is refreshed from the external mode every
     * step — so boundary water enters through the same placement, composition
     * and conservation machinery a tributary inflow already uses.
     *
     * **Why a source and not an edge flux.** The external mode solves a
     * depth-averaged problem, so all it knows about a boundary is the
     * depth-integrated volume that crossed it. How that volume distributes over
     * the water column is a modelling choice in either formulation, and a
     * column source makes the choice explicit and selectable instead of burying
     * it in a flux reconstruction. It also keeps the two modes consistent for
     * free: the vertical velocity is diagnosed from continuity *including*
     * sources, so the surface-flux residual that guards mode consistency
     * remains valid with a boundary open.
     *
     * \param[in]  placement How inflow spreads over the column. DensityMatched
     *                       by default, because a tributary entering a
     *                       stratified reservoir plunges to its own density
     *                       rather than mixing across the surface.
     * \param[out] message   Failure description when returning false.
     */
    [[nodiscard]] bool bindOpenBoundaries(
        Forcing::Placement placement, std::string &message);

    //! As above, with the density-matched default.
    [[nodiscard]] bool bindOpenBoundaries(std::string &message)
    {
        return bindOpenBoundaries(Forcing::Placement::DensityMatched, message);
    }

    /*!
     * \brief Prescribes what an open boundary's inflow carries.
     *
     * Unset, a boundary carries **ambient** values: inflow arrives at the
     * receiving column's own concentration, so an unconfigured open boundary is
     * transparent to constituents rather than a source of distilled water. That
     * is the neutral choice, not a silent one — but it is also not a tributary,
     * and a model with a real inflow should say what is in it.
     *
     * \param[in]  boundary  Index returned by external().addBoundary().
     * \param[in]  speciesId Constituent, or the reserved temperature id.
     * \param[in]  value     Concentration, or temperature in deg C.
     * \param[out] message   Failure description when returning false.
     */
    [[nodiscard]] bool setBoundaryValue(int boundary, const std::string &speciesId,
                                        double value, std::string &message);

    //! Returns an open boundary to carrying ambient values.
    [[nodiscard]] bool clearBoundaryValue(int boundary, const std::string &speciesId,
                                          std::string &message);

    /*!
     * \brief Volume that crossed \p boundary into the domain last step (m³).
     *
     * Negative for net export. Reported so a constituent budget can be closed
     * against the boundary rather than assumed.
     */
    [[nodiscard]] double boundaryVolume(int boundary) const;

    /*!
     * \brief Index of \p id in this solver's registry, or -1.
     *
     * Returns -1 before initialize(), so a caller that needs the registry can
     * check this first and never reach species() with nothing behind it.
     */
    [[nodiscard]] int speciesIndex(const std::string &id) const noexcept;

    /*!
     * \brief The species registry this solver was initialized with.
     *
     * Read-only, and deliberately so: the reaction network is bound through
     * compileReactions() rather than by handing the registry to a caller, so
     * that a network cannot be resolved against a different registry than the
     * state was allocated with.
     *
     * \pre initialize() has succeeded.
     */
    [[nodiscard]] const State::SpeciesRegistry &species() const noexcept
    {
        return *m_species;
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
    void   updateVerticalMixing(double dt);
    void   advanceMomentum(double dt);
    void   applyModeAdjustment();
    void   computeLayerFluxes();
    void   diagnoseVerticalFlux(double dt);
    void   transportScalars(double dt);
    void   updateCellDepth();
    void   refreshBoundarySources(double dt);
    void   applyBoundaryComposition(int species);

    /*!
     * \brief One open-boundary edge group: a barotropic boundary meeting a
     *        column, and the Forcing::Source that carries its water.
     */
    struct BoundaryLink
    {
        int                  boundary = -1; //!< Barotropic boundary index.
        int64_t              column   = 0;  //!< Column it feeds.
        int                  source   = -1; //!< Source index in m_sources.
        Forcing::Placement   inflowPlacement = Forcing::Placement::Uniform;
        std::vector<int64_t> edges;         //!< Edges of this pair.
        double               volume   = 0.0; //!< Volume in last step (m³).
    };

    std::vector<BoundaryLink> m_boundaryLinks;

    //! Prescribed inflow values, by (boundary, species); NaN means ambient.
    std::vector<double> m_boundaryValue;
    void   applyHorizontalDiffusion(double dt, std::span<double> values);
    void   applySurfaceHeatFlux(double dt);
    double bottomDragCoefficient(int64_t column) const;

    Mesh::LayeredMesh            *m_mesh    = nullptr;
    State::StateStore            *m_state   = nullptr;
    const State::SpeciesRegistry *m_species = nullptr;

    InternalModeOptions               m_options;
    mutable InternalModeDiagnostics m_diagnostics;
    BarotropicSolver                  m_external;
    Transport::EulerianTransport      m_transport;
    Forcing::SourceSet                m_sources;
    Forcing::StructureSet             m_structures;
    Forcing::OutletGroupSet           m_groups;
    Forcing::TransferSet              m_transfers;
    Transport::SettlingSet            m_settling;
    Kinetics::ReactionNetwork         m_reactions;
    Kinetics::LightField              m_light;
    Kinetics::LightField              m_solar;
    TurbulenceClosure                 m_turbulence;
    Forcing::MeteorologyField         m_meteorology;
    bool                              m_useMeteorology   = false;
    double                            m_solarAttenuation = 0.45;

    /*!
     * \brief Per-column net shortwave handed to the solar LightField (W/m²).
     *
     * A member rather than a local so the per-step fill does not allocate.
     */
    std::vector<double> m_incidentShortwave;

    /*!
     * \brief Per-column surface wind stress (Pa), refreshed each step.
     *
     * Holds `InternalModeOptions::windStressX/Y` when meteorology is off and
     * the derived `windStress()` when it is on.
     */
    std::vector<double> m_windStressX;
    std::vector<double> m_windStressY;

    //! Per-thread scratch for the closure's per-interface inputs.
    std::vector<double> m_shearSquared;
    std::vector<double> m_buoyancyFrequency;
    Kinetics::LightOptions            m_lightOptions;
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
    std::vector<double> m_cellDepth; //!< Cell-centre depth below the surface (m).

    /*!
     * \brief Kinetics::CellPosition flags for every cell.
     *
     * Fixed for the life of the solver: the layer count does not change, so the
     * top and bottom cells of each column are always the top and bottom cells.
     * Built once rather than recomputed, and handed to the reaction network so
     * it can apply scoped rate laws without knowing the mesh's indexing
     * convention.
     */
    std::vector<unsigned char> m_cellPosition;
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
