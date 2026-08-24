// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   barotropicsolver.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  BarotropicSolver — the explicit external mode: depth-averaged
 *         shallow-water flow on the tri-quad mesh.
 *
 * This is the fast mode of the barotropic-baroclinic splitting. Surface
 * gravity waves travel at sqrt(gH) — 22 m/s in 50 m of water — so resolving
 * them explicitly costs a time step of order a second. Rather than drag the
 * whole 3-D state along at that rate, mode splitting confines the surface
 * wave to this cheap 2-D solver and subcycles it inside each internal step
 * Everything here is written to be that inner loop: flat face loops,
 * no linear solve, no allocation per step.
 *
 * Numerics:
 * - Cell-centred finite volume in conserved variables (h, hu, hv), one cell
 *   per water column, on the unique-edge list so each face is visited once
 *   and its flux scattered equal-and-opposite to both sides. Mass is
 *   conserved to machine precision by construction.
 * - **HLLC** approximate Riemann solver at faces, evaluated in the
 *   face-normal frame. The contact wave matters here: an HLL solver smears
 *   tangential momentum across shear, which a lake with a wind-driven gyre
 *   feels immediately.
 * - **Audusse hydrostatic reconstruction** for well-balancedness. A naive
 *   scheme with a bed source term produces spurious currents over sloping
 *   bathymetry — the "lake at rest" that will not stay at rest. Here the
 *   interface depths are reconstructed against the higher of the two beds
 *   and the bed force is applied as a matching interface correction, so
 *   still water over arbitrary bathymetry is preserved to machine precision.
 *   That property is a hard test, not an aspiration: see the Phase 2 gate.
 * - **Semi-implicit Manning friction**, which is unconditionally stable and
 *   cannot reverse the flow no matter how thin the film — the failure mode
 *   that plagues explicit friction in drying cells.
 * - **Positivity-preserving wet/dry**: depths are clamped at zero, momentum
 *   is dropped below a dry threshold, and the Riemann wave speeds switch to
 *   the dry-front estimates when either side is dry.
 *
 * Boundaries are reflective walls in Phase 2; inflows, withdrawals, and open
 * boundaries arrive with the internal mode in Phase 3 (F3.4).
 *
 * \threadsafety One solver per thread; step() mutates the owned state.
 */

#ifndef FVQual_HYDRO_BAROTROPICSOLVER_H
#define FVQual_HYDRO_BAROTROPICSOLVER_H

#include "fvqual/fvqual_export.h"
#include "fvqual/mesh/horizontalmesh.h"

#include <cstdint>
#include <functional>
#include <span>
#include <string>
#include <vector>

namespace FVQual::Hydro {

/*!
 * \brief Physical and numerical settings of the external mode.
 *
 * SI units throughout.
 */
struct FVQual_EXPORT BarotropicOptions
{
    double gravity      = 9.80665; //!< Gravitational acceleration (m/s²).
    double waterDensity = 998.2;   //!< Reference water density (kg/m³).

    /*!
     * \brief Uniform Manning roughness (s/m^(1/3)); override per column with
     *        setManning().
     */
    double manningN = 0.03;

    /*!
     * \brief Depth below which a column is treated as dry (m). Momentum is
     *        zeroed there and the cell contributes no wave speed.
     */
    double dryDepth = 1.0e-3;

    /*!
     * \brief Courant number for the adaptive time step. Explicit first-order
     *        finite volume is stable to 1; 0.45 leaves headroom for the
     *        wetting front, where the linearized limit is optimistic.
     */
    double cflNumber = 0.45;

    double minTimeStep = 1.0e-6; //!< Floor on the adaptive step (s).
    double maxTimeStep = 60.0;   //!< Ceiling on the adaptive step (s).

    double windStressX = 0.0; //!< Surface wind stress, x (Pa).
    double windStressY = 0.0; //!< Surface wind stress, y (Pa).

    bool frictionEnabled = true; //!< Apply bed friction.
};

/*!
 * \brief Per-step diagnostics.
 */
struct FVQual_EXPORT BarotropicDiagnostics
{
    double  timeStep    = 0.0; //!< Step actually taken (s).
    double  maxCourant  = 0.0; //!< Largest cell Courant number.
    double  totalVolume = 0.0; //!< Domain water volume (m³).
    int64_t wetColumns  = 0;   //!< Columns above the dry threshold.
    int64_t steps       = 0;   //!< Steps taken by the last advance().
};

/*!
 * \brief Explicit depth-averaged shallow-water solver.
 */
class FVQual_EXPORT BarotropicSolver
{
public:
    BarotropicSolver() = default;

    /*!
     * \brief Binds the solver to a mesh and allocates its state.
     *
     * The mesh must outlive the solver. All columns start dry; call
     * setUniformSurface() or setSurface() to establish an initial condition.
     *
     * \param[in]  mesh    Horizontal mesh with face bed elevations set.
     * \param[in]  options Physical and numerical settings.
     * \param[out] message Failure description when returning false.
     */
    [[nodiscard]] bool initialize(const Mesh::HorizontalMesh &mesh,
                                  const BarotropicOptions &options,
                                  std::string &message);

    /*!
     * \brief Sets a spatially uniform water-surface elevation at rest.
     *        Columns whose bed lies above it become dry.
     */
    void setUniformSurface(double elevation);

    /*!
     * \brief Sets the water-surface elevation per column, at rest.
     * \param[in]  elevation One value per column (m).
     * \param[out] message   Failure description when returning false.
     */
    [[nodiscard]] bool setSurface(std::span<const double> elevation, std::string &message);

    /*!
     * \brief Sets the depth-averaged velocity per column.
     */
    [[nodiscard]] bool setVelocity(std::span<const double> u, std::span<const double> v,
                                   std::string &message);

    /*!
     * \brief Overrides the uniform Manning roughness with per-column values.
     */
    [[nodiscard]] bool setManning(std::span<const double> manning, std::string &message);

    /*!
     * \brief Sets a depth-distributed body acceleration per column (m/s²).
     *
     * This is how the internal mode hands the external mode its share of the
     * baroclinic pressure gradient. Without it the two modes are decoupled:
     * the external mode would keep a flat surface while density differences
     * drive the interior, and the subsequent mode adjustment — which forces
     * the internal depth-mean onto the external velocity — would then cancel
     * the barotropic response that should have developed. A lock exchange
     * simply refuses to propagate.
     *
     * Unlike wind stress, which acts on the surface, this force acts through
     * the water column and so scales with depth.
     */
    [[nodiscard]] bool setBodyAcceleration(std::span<const double> ax,
                                           std::span<const double> ay,
                                           std::string &message);

    /*!
     * \brief Clears any body acceleration.
     */
    void clearBodyAcceleration() noexcept;

    /*!
     * \brief Sets a volumetric source per column (m³/s), positive adding water.
     *
     * This is how inflows and withdrawals reach the free surface. The internal
     * mode must hand this solver exactly the total it distributes over the
     * cells, or the two modes disagree about how much water there is and the
     * free-surface continuity residual stops vanishing.
     *
     * Momentum treatment differs by sign, and the asymmetry is physical:
     * inflow arrives carrying no momentum of its own, so it dilutes the local
     * velocity, while outflow leaves at the local velocity and therefore
     * changes the depth without changing the speed.
     */
    /*!
     * \brief What an open edge does with the water that reaches it.
     */
    enum class BoundaryType
    {
        //! Reflective wall: the mirror state. The default for every edge.
        Wall,

        //! Prescribed inflow, stated as discharge per unit edge length.
        Inflow,

        //! Prescribed water-surface elevation.
        Stage,

        //! Zero-gradient: whatever arrives leaves, with no reflection.
        Transmissive
    };

    /*!
     * \brief One open boundary.
     *
     * How many quantities may be imposed is not a modelling preference, it is
     * counted by the characteristics. In the outward-normal direction the
     * shallow-water equations carry information at speeds `u ± c`. When the
     * flow is subcritical exactly one of those leaves the domain, so exactly
     * one condition may be imposed and the other is supplied by the outgoing
     * Riemann invariant `R = u + 2c` evaluated inside. Impose two and the
     * problem is over-determined; the surplus reflects back into the domain as
     * spurious waves. Impose none and it is under-determined.
     *
     * The solver therefore reads the local Froude number and adapts:
     * subcritical inflow takes the discharge and gets its depth from the
     * invariant; subcritical outflow takes the stage and gets its velocity the
     * same way; supercritical outflow takes nothing, because nothing outside
     * can reach back in; supercritical inflow needs both, which is why
     * \c depth exists.
     */
    struct BoundaryCondition
    {
        std::string  id;
        BoundaryType type = BoundaryType::Wall;

        /*!
         * \brief Inflow per unit edge length (m²/s), positive into the domain.
         */
        double specificDischarge = 0.0;

        /*!
         * \brief Water-surface elevation (m) for a Stage boundary.
         */
        double stage = 0.0;

        /*!
         * \brief Inflow depth (m). Only needed for a *supercritical* inflow,
         *        where both characteristics enter and the invariant supplies
         *        nothing. Zero leaves the depth to the invariant, which is
         *        correct whenever the inflow is subcritical.
         */
        double depth = 0.0;
    };

    /*!
     * \brief Registers a boundary condition over a set of boundary edges.
     *
     * \param[in]  condition The condition.
     * \param[in]  edges     Edge indices; each must be an unpaired mesh edge.
     * \param[out] message   Failure description when returning -1.
     * \return Its index, or -1 when the specification is invalid.
     */
    int addBoundary(const BoundaryCondition &condition, std::span<const int64_t> edges,
                    std::string &message);

    [[nodiscard]] int boundaryCount() const noexcept
    {
        return static_cast<int>(m_boundaries.size());
    }

    /*!
     * \brief Sets an inflow rate mid-run (m²/s per unit edge length).
     */
    [[nodiscard]] bool setBoundaryDischarge(int index, double specificDischarge,
                                            std::string &message);

    /*!
     * \brief Sets a boundary stage mid-run (m).
     */
    [[nodiscard]] bool setBoundaryStage(int index, double stage, std::string &message);

    /*!
     * \brief Volume flux through a boundary on the last step (m³/s), signed
     *        positive into the domain.
     *
     * The quantity a mass budget needs, and the one that says whether a
     * boundary delivered what it was asked for.
     */
    [[nodiscard]] double boundaryFlux(int index) const noexcept
    {
        return m_boundaryFlux[static_cast<size_t>(index)];
    }

    /*!
     * \brief Whether any edge carries a condition other than a wall.
     */
    [[nodiscard]] bool hasOpenBoundaries() const noexcept
    {
        return !m_boundaries.empty();
    }

    /*!
     * \brief Boundary index for each edge, or -1 for an interior or wall edge.
     *
     * Exposed so the internal mode can work out which column each open-boundary
     * edge feeds, and with what composition, without duplicating the edge-to-
     * boundary assignment made here.
     */
    [[nodiscard]] const std::vector<int> &edgeBoundary() const noexcept
    {
        return m_edgeBoundary;
    }

    /*!
     * \brief Callback that refreshes the column sources before each substep.
     *
     * \param surfaceElevation Current water-surface elevation per column (m).
     * \param dt               The substep about to be taken (s).
     * \param columnSource     Per-column discharge to overwrite (m³/s).
     * \param message          Failure description when returning false.
     */
    using ColumnSourceUpdater =
        std::function<bool(std::span<const double> surfaceElevation, double dt,
                           std::span<double> columnSource, std::string &message)>;

    /*!
     * \brief Installs a hook called at the top of every substep.
     *
     * Sources whose discharge depends on the water level cannot be held fixed
     * across a subcycle: the level is precisely what the subcycle is changing.
     * A spillway rated once for the whole interval releases at its starting
     * head long after the head has gone, and empties the reservoir past its own
     * crest. This hook lets such a source be re-evaluated against the surface
     * as it actually evolves. Prescribed sources need no hook and should not
     * install one.
     */
    void setColumnSourceUpdater(ColumnSourceUpdater updater)
    {
        m_columnSourceUpdater = std::move(updater);
    }

    [[nodiscard]] bool setColumnSource(std::span<const double> discharge,
                                       std::string &message);

    /*!
     * \brief Clears any column source.
     */
    void clearColumnSource() noexcept;

    /*!
     * \brief Largest stable time step for the current state (s).
     *
     * The unstructured CFL estimate: for each wet cell, area divided by the
     * sum over its edges of length times (|u.n| + sqrt(gh)), scaled by the
     * Courant number. Returns maxTimeStep when the domain is entirely dry or
     * motionless.
     */
    [[nodiscard]] double stableTimeStep() const;

    /*!
     * \brief Advances one explicit step of exactly \p dt seconds.
     * \param[in]  dt      Time step (s, > 0).
     * \param[out] message Failure description when returning false.
     */
    [[nodiscard]] bool step(double dt, std::string &message);

    /*!
     * \brief Advances by \p duration seconds, substepping at the CFL limit.
     *
     * The final substep is trimmed so the solver lands exactly on the target
     * time — the contract the internal mode will rely on when it subcycles
     * this solver.
     *
     * \param[in]  duration Interval to advance (s, >= 0).
     * \param[out] message  Failure description when returning false.
     */
    [[nodiscard]] bool advance(double duration, std::string &message);

    // ── State ─────────────────────────────────────────────────────────────

    [[nodiscard]] int64_t columnCount() const noexcept { return m_columnCount; }

    /*!
     * \brief Water depth per column (m).
     */
    [[nodiscard]] const std::vector<double> &depth() const noexcept { return m_h; }

    /*!
     * \brief Water-surface elevation per column (m) — bed plus depth.
     */
    [[nodiscard]] const std::vector<double> &surfaceElevation() const noexcept
    {
        return m_eta;
    }

    /*!
     * \brief Depth-averaged velocity per column (m/s).
     */
    [[nodiscard]] const std::vector<double> &velocityX() const noexcept { return m_u; }
    [[nodiscard]] const std::vector<double> &velocityY() const noexcept { return m_v; }

    /*!
     * \brief Whether a column is wet.
     */
    [[nodiscard]] bool isWet(int64_t column) const noexcept
    {
        return m_h[static_cast<size_t>(column)] > m_options.dryDepth;
    }

    /*!
     * \brief Total water volume (m³) — the quantity the mass-closure gate
     *        tracks.
     */
    [[nodiscard]] double totalVolume() const noexcept;

    /*!
     * \brief Elapsed simulated time (s).
     */
    [[nodiscard]] double time() const noexcept { return m_time; }

    // ── Flux accumulation for mode splitting ──────────────────────────────
    //
    // The internal mode must transport scalars on fluxes consistent with the
    // free surface the external mode actually produced. Instantaneous fluxes
    // will not do: the external mode takes many substeps per internal step,
    // and the surface wave it resolves largely cancels over that interval.
    // What the internal mode needs is the *time-averaged* transport, which is
    // what these accumulate.

    /*!
     * \brief Clears the accumulated edge volumes and restarts the interval.
     */
    void resetFluxAccumulation() noexcept;

    /*!
     * \brief Volume that has crossed each edge since the last reset (m³),
     *        signed positive from the edge's left face toward its right.
     */
    [[nodiscard]] const std::vector<double> &accumulatedEdgeVolume() const noexcept
    {
        return m_edgeVolumeAccumulated;
    }

    /*!
     * \brief Duration covered by the accumulation (s).
     */
    [[nodiscard]] double accumulationInterval() const noexcept
    {
        return m_accumulationInterval;
    }

    [[nodiscard]] const BarotropicDiagnostics &diagnostics() const noexcept
    {
        return m_diagnostics;
    }

    [[nodiscard]] const BarotropicOptions &options() const noexcept { return m_options; }

private:
    /*!
     * \brief Accumulates face fluxes and bed source terms into m_dh/m_dhu/m_dhv.
     */
    void accumulateFluxes();

    /*!
     * \brief Applies wind stress and semi-implicit Manning friction.
     */
    void applySourceTerms(double dt);

    /*!
     * \brief Recomputes depth, surface elevation, and velocity from the
     *        conserved state, enforcing positivity and the dry threshold.
     */
    void reconstructPrimitives();

    const Mesh::HorizontalMesh *m_mesh = nullptr;
    BarotropicOptions           m_options;
    BarotropicDiagnostics       m_diagnostics;

    int64_t m_columnCount = 0;
    double  m_time        = 0.0;

    // Conserved state.
    std::vector<double> m_h;
    std::vector<double> m_hu;
    std::vector<double> m_hv;

    // Primitive state, refreshed from the conserved state each step.
    std::vector<double> m_eta;
    std::vector<double> m_u;
    std::vector<double> m_v;

    // Flux accumulators, allocated once.
    std::vector<double> m_dh;
    std::vector<double> m_dhu;
    std::vector<double> m_dhv;

    // Per-edge volume flux of the current step (m³/s), and its time integral
    // since the last reset (m³).
    std::vector<double> m_edgeVolumeFlux;
    std::vector<double> m_edgeVolumeAccumulated;
    double              m_accumulationInterval = 0.0;

    std::vector<double> m_manning;

    // Depth-distributed body acceleration (m/s²), empty when unused.
    std::vector<double> m_bodyAccelerationX;
    std::vector<double> m_bodyAccelerationY;

    // Volumetric source per column (m³/s), empty when unused.
    std::vector<double> m_columnSource;
    ColumnSourceUpdater m_columnSourceUpdater;

    // Per-edge flux storage, so the edge sweep does not scatter into cells
    // and both passes can be threaded (F7.1a).
    std::vector<double> m_edgeFluxMass;
    std::vector<double> m_edgeFluxX;
    std::vector<double> m_edgeFluxY;
    std::vector<double> m_edgeCorrectionL;
    std::vector<double> m_edgeCorrectionR;

    std::vector<BoundaryCondition> m_boundaries;

    //! Edge -> boundary index, or -1 for a reflective wall.
    std::vector<int>    m_edgeBoundary;
    std::vector<double> m_boundaryFlux;

    /*!
     * \brief Builds the exterior state an open edge should see.
     */
    void ghostState(const BoundaryCondition &condition, double hL, double unL,
                    double utL, double zL, double &hR, double &unR,
                    double &utR) const;
};

} // namespace FVQual::Hydro

#endif // FVQual_HYDRO_BAROTROPICSOLVER_H
