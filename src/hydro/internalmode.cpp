// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/hydro/internalmode.h"

#include "fvqual/hydro/equationofstate.h"

#include "fvqual/parallel.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>

namespace FVQual::Hydro {

namespace {

/*!
 * \brief Thomas algorithm for a tridiagonal system of \p n equations.
 *
 * The vertical column solves are the only implicit part of the model, and
 * they are tridiagonal by construction of the (column x layer) indexing, so
 * this is O(n) with no allocation and no pivoting — the systems are
 * diagonally dominant because the diffusion operator makes them so.
 */
void solveTridiagonal(int n, const double *lower, const double *diagonal,
                      const double *upper, const double *rhs, double *solution,
                      double *sweep)
{
    if (n <= 0)
        return;

    double beta = diagonal[0];
    solution[0] = rhs[0] / beta;

    for (int i = 1; i < n; ++i)
    {
        sweep[i]    = upper[i - 1] / beta;
        beta        = diagonal[i] - lower[i] * sweep[i];
        solution[i] = (rhs[i] - lower[i] * solution[i - 1]) / beta;
    }

    for (int i = n - 2; i >= 0; --i)
        solution[i] -= sweep[i + 1] * solution[i + 1];
}

/*!
 * \brief Adds its lifetime to a counter, in nanoseconds.
 */
class ScopedTimer
{
public:
    explicit ScopedTimer(int64_t &sink)
        : m_sink(sink), m_start(std::chrono::steady_clock::now())
    {
    }

    ~ScopedTimer()
    {
        m_sink += std::chrono::duration_cast<std::chrono::nanoseconds>(
                      std::chrono::steady_clock::now() - m_start)
                      .count();
    }

    ScopedTimer(const ScopedTimer &)            = delete;
    ScopedTimer &operator=(const ScopedTimer &) = delete;

private:
    int64_t                                            &m_sink;
    std::chrono::time_point<std::chrono::steady_clock>  m_start;
};

} // namespace

InternalModeSolver::InternalModeSolver()
    : m_pressure(std::make_unique<HydrostaticPressureClosure>())
{
}

InternalModeSolver::~InternalModeSolver() = default;

void InternalModeSolver::setPressureClosure(std::unique_ptr<IPressureClosure> closure)
{
    if (closure)
        m_pressure = std::move(closure);
}

bool InternalModeSolver::initialize(Mesh::LayeredMesh &mesh, State::StateStore &state,
                                    const State::SpeciesRegistry &species,
                                    const InternalModeOptions &options,
                                    std::string &message)
{
    if (mesh.columnCount() <= 0 || mesh.layerCount() <= 0)
    {
        message = "internal mode needs a layered mesh with at least one layer";
        return false;
    }

    if (state.cellCount() != mesh.cellCount())
    {
        message = "state store was allocated for a different mesh";
        return false;
    }

    if (options.referenceDensity <= 0.0 || options.gravity <= 0.0)
    {
        message = "gravity and reference density must be positive";
        return false;
    }

    if (options.cflNumber <= 0.0 || options.cflNumber > 1.0)
    {
        message = "the Courant number must lie in (0, 1]";
        return false;
    }

    m_mesh    = &mesh;
    m_state   = &state;
    m_species = &species;
    m_options = options;
    m_time    = 0.0;

    m_temperatureIndex = species.temperatureIndex();
    m_salinityIndex    = options.salinitySpeciesId.empty()
                             ? -1
                             : species.indexOf(options.salinitySpeciesId);

    if (!options.salinitySpeciesId.empty() && m_salinityIndex < 0)
    {
        message = "salinity species '" + options.salinitySpeciesId + "' is not registered";
        return false;
    }

    // The three-dimensional mode has no wetting and drying yet (F3.4), so a
    // dry column would silently divide by a zero thickness. Refuse instead.
    for (int64_t c = 0; c < mesh.columnCount(); ++c)
    {
        if (mesh.columnDepth()[static_cast<size_t>(c)] <= 0.0)
        {
            message = "column " + std::to_string(c) +
                      " is dry; the internal mode requires every column to be wet "
                      "(three-dimensional wetting and drying is deferred to F3.4)";
            return false;
        }
    }

    BarotropicOptions externalOptions;
    externalOptions.gravity         = options.gravity;
    externalOptions.waterDensity    = options.referenceDensity;
    externalOptions.manningN        = options.manningN;
    externalOptions.windStressX     = options.windStressX;
    externalOptions.windStressY     = options.windStressY;
    externalOptions.frictionEnabled = true;

    if (!m_external.initialize(mesh.horizontal(), externalOptions, message))
        return false;

    if (!m_external.setSurface(mesh.surfaceElevation(), message))
        return false;

    // Both modes must agree about how much water a column holds, so the storage
    // curve goes to the external solver as well as the mesh. Without this the
    // external mode would raise the surface on a constant area while the mesh
    // computed volumes on a varying one, and the residual that guards mode
    // consistency would stop vanishing — which is exactly what the gate checks.
    if (mesh.hasHypsography())
    {
        for (int64_t c = 0; c < mesh.columnCount(); ++c)
        {
            const std::vector<double> &elevation = mesh.hypsographyElevation(c);
            const std::vector<double> &planArea  = mesh.hypsographyArea(c);

            if (elevation.empty())
                continue;

            if (!m_external.setHypsography(c, elevation, planArea, message))
                return false;
        }
    }

    const size_t cells      = static_cast<size_t>(mesh.cellCount());
    const size_t interfaces = static_cast<size_t>(mesh.interfaceSlotCount());
    const size_t edgeLayers = static_cast<size_t>(mesh.horizontal().edgeCount()) *
                              static_cast<size_t>(mesh.layerCount());

    m_accelerationX.assign(cells, 0.0);
    m_accelerationY.assign(cells, 0.0);
    m_km.assign(interfaces, options.backgroundViscosity);
    m_kz.assign(interfaces, options.backgroundDiffusivity);
    m_cellDepth.assign(cells, 0.0);

    // A one-layer column is both surface and bed, and the flags say so rather
    // than picking one — a shallow pond really does reaerate and exert sediment
    // demand on the same water.
    m_cellPosition.assign(cells, Kinetics::CellInterior);
    for (int64_t c = 0; c < mesh.columnCount(); ++c)
    {
        m_cellPosition[static_cast<size_t>(mesh.cell(c, 0))] |= Kinetics::CellAtSurface;
        m_cellPosition[static_cast<size_t>(mesh.cell(c, mesh.layerCount() - 1))] |=
            Kinetics::CellAtBed;
    }

    // Sized unconditionally, even though the march itself only runs when a rate
    // law asks for irradiance. The alternative is allocating it when the network
    // is compiled, which puts an ordering trap in the API: the natural thing is
    // to set self-shading coefficients before compiling, and that would then
    // write into a field that does not exist yet.
    m_turbulence.initialize(mesh.columnCount(), mesh.layerCount());

    {
        const size_t interfaces = static_cast<size_t>(mesh.layerCount() + 1);
        const size_t threads    = static_cast<size_t>(maxThreads());

        m_shearSquared.assign(interfaces * threads, 0.0);
        m_buoyancyFrequency.assign(interfaces * threads, 0.0);
    }

    m_light.initialize(mesh, species);
    m_solar.initialize(mesh, species);
    m_settling.initialize(species);
    m_edgeLayerFlux.assign(edgeLayers, 0.0);
    m_verticalFlux.assign(interfaces, 0.0);
    m_oldVolume.assign(cells, 0.0);
    m_workU.assign(cells, 0.0);
    m_workV.assign(cells, 0.0);
    m_scalarWork.assign(cells, 0.0);
    m_cellSource.assign(cells, 0.0);
    m_columnSource.assign(static_cast<size_t>(mesh.columnCount()), 0.0);

    // Vertical-solve workspace, one slice per thread. A single shared buffer
    // is a data race the moment the column loop is threaded, and one that
    // would corrupt results quietly rather than crash — every column would
    // solve against whatever tridiagonal its neighbours were mid-way through
    // assembling.
    const size_t layers  = static_cast<size_t>(mesh.layerCount());
    const size_t scratch = layers * static_cast<size_t>(maxThreads());

    m_lower.assign(scratch, 0.0);
    m_diagonal.assign(scratch, 0.0);
    m_upper.assign(scratch, 0.0);
    m_rhs.assign(scratch, 0.0);
    m_solution.assign(scratch, 0.0);
    m_sweep.assign(scratch, 0.0);

    m_transport.initialize(mesh);

    // Head-driven structures follow the free surface inside the subcycle. Only
    // the column totals are refreshed here: the external mode consumes nothing
    // else, and repeating the vertical placement — envelopes, density, layer
    // weights — a hundred times per internal step would cost far more than the
    // placement moves in that time. Transport gets the properly placed sources
    // once the subcycle is done (step 4b).
    m_external.setColumnSourceUpdater(
        [this](std::span<const double> surfaceElevation, double substep,
               std::span<double> columnSource, std::string &error) -> bool {
            if (m_structures.count() == 0)
                return true;

            if (!m_structures.update(surfaceElevation, m_sources, error))
                return false;

            m_sources.columnTotals(columnSource);

            for (int s = 0; s < m_structures.count(); ++s)
            {
                m_structureVolume[static_cast<size_t>(s)] +=
                    m_structures.discharge(s) * substep;
            }

            return true;
        });

    updateDensity();
    updateVerticalMixing(0.0);

    m_diagnostics             = InternalModeDiagnostics{};
    m_diagnostics.totalVolume = mesh.totalVolume();

    message.clear();
    return true;
}

void InternalModeSolver::updateDensity()
{
    State::StateStore &state = *m_state;

    const int64_t cells = state.cellCount();

    FVQual_PARALLEL_FOR
    for (int64_t i = 0; i < cells; ++i)
    {
        const double temperature =
            m_temperatureIndex >= 0 ? state.speciesAt(m_temperatureIndex, i) : 4.0;

        const double salinity =
            m_salinityIndex >= 0 ? std::max(0.0, state.speciesAt(m_salinityIndex, i)) : 0.0;

        state.density[static_cast<size_t>(i)] =
            EquationOfState::density(temperature, salinity);
    }
}

void InternalModeSolver::updateVerticalMixing(double dt)
{
    const Mesh::LayeredMesh  &mesh  = *m_mesh;
    const State::StateStore  &state = *m_state;

    const int64_t columns = mesh.columnCount();
    const int     layers  = mesh.layerCount();

    // ── k-epsilon closure (F5.4) ──────────────────────────────────────────
    //
    // `dt > 0` also selects the Richardson path during initialize(), which has
    // no step to advance over: the closure needs a seed, and Ri damping is a
    // reasonable one.
    if (m_options.useTurbulenceClosure && dt > 0.0)
    {
        const int interfaces = layers + 1;

        const std::vector<double> &thickness = mesh.thickness();

        FVQual_PARALLEL_FOR
        for (int64_t c = 0; c < columns; ++c)
        {
            const size_t slot = static_cast<size_t>(threadIndex() * interfaces);

            double *shearSquared = m_shearSquared.data() + slot;
            double *buoyancy     = m_buoyancyFrequency.data() + slot;

            // Surface and bed interfaces carry no interior gradient; the walls
            // enter through the friction velocities instead.
            shearSquared[0] = 0.0;
            buoyancy[0]     = 0.0;
            shearSquared[layers] = 0.0;
            buoyancy[layers]     = 0.0;

            for (int k = 1; k < layers; ++k)
            {
                const size_t above = static_cast<size_t>(mesh.cell(c, k - 1));
                const size_t below = static_cast<size_t>(mesh.cell(c, k));

                const double spacing =
                    0.5 * (thickness[above] + thickness[below]);

                if (spacing <= 0.0)
                {
                    shearSquared[k] = 0.0;
                    buoyancy[k]     = 0.0;
                    continue;
                }

                buoyancy[k] = m_options.gravity / m_options.referenceDensity *
                              (state.density[below] - state.density[above]) / spacing;

                const double shearU = (state.u[above] - state.u[below]) / spacing;
                const double shearV = (state.v[above] - state.v[below]) / spacing;

                shearSquared[k] = shearU * shearU + shearV * shearV;
            }

            // Friction velocities close the closure at both walls. The surface
            // one is the wind stress F3.7 derives; the bed one is the quadratic
            // drag the momentum solve already applies, so the two schemes cannot
            // disagree about how rough the bed is.
            double surfaceStress = 0.0;

            if (!m_windStressX.empty())
            {
                surfaceStress = std::hypot(m_windStressX[static_cast<size_t>(c)],
                                           m_windStressY[static_cast<size_t>(c)]);
            }

            const double frictionSurface =
                std::sqrt(surfaceStress / m_options.referenceDensity);

            const size_t bottom = static_cast<size_t>(mesh.cell(c, layers - 1));

            const double frictionBed =
                std::sqrt(bottomDragCoefficient(c)) *
                std::hypot(state.u[bottom], state.v[bottom]);

            const size_t first = static_cast<size_t>(mesh.interfaceSlot(c, 0));

            m_turbulence.advanceColumn(
                c, dt, &thickness[static_cast<size_t>(mesh.cell(c, 0))],
                shearSquared, buoyancy, frictionSurface, frictionBed,
                &m_km[first], &m_kz[first]);
        }

        // The closure carries its own memory, so there is no Richardson number
        // to report. Left at zero rather than stale.
        m_diagnostics.maxRichardson = 0.0;

        return;
    }

    // Accumulated into a local so the threaded loop below does not
    // read-modify-write a shared field. It is only a diagnostic, but a racing
    // diagnostic reports a number nobody can reproduce, which is worse than no
    // number at all.
    double maxRichardson = 0.0;

    FVQual_PARALLEL_FOR_REDUCE(max : maxRichardson)
    for (int64_t c = 0; c < columns; ++c)
    {
        // Surface and bed interfaces carry no interior mixing.
        m_km[static_cast<size_t>(mesh.interfaceSlot(c, 0))] = 0.0;
        m_kz[static_cast<size_t>(mesh.interfaceSlot(c, 0))] = 0.0;
        m_km[static_cast<size_t>(mesh.interfaceSlot(c, layers))] = 0.0;
        m_kz[static_cast<size_t>(mesh.interfaceSlot(c, layers))] = 0.0;

        for (int k = 1; k < layers; ++k)
        {
            const size_t above = static_cast<size_t>(mesh.cell(c, k - 1));
            const size_t below = static_cast<size_t>(mesh.cell(c, k));
            const size_t slot  = static_cast<size_t>(mesh.interfaceSlot(c, k));

            const double spacing =
                0.5 * (mesh.thickness()[above] + mesh.thickness()[below]);

            if (spacing <= 0.0)
            {
                m_km[slot] = m_options.backgroundViscosity;
                m_kz[slot] = m_options.backgroundDiffusivity;
                continue;
            }

            // Buoyancy frequency: positive when density increases downward.
            const double buoyancy = m_options.gravity / m_options.referenceDensity *
                                    (state.density[below] - state.density[above]) / spacing;

            const double shearU = (state.u[above] - state.u[below]) / spacing;
            const double shearV = (state.v[above] - state.v[below]) / spacing;
            const double shear  = shearU * shearU + shearV * shearV;

            if (buoyancy < 0.0)
            {
                // Statically unstable: convective adjustment mixes hard until
                // the column restratifies.
                m_km[slot] = m_options.convectiveMixing;
                m_kz[slot] = m_options.convectiveMixing;
                continue;
            }

            const double richardson = buoyancy / std::max(shear, 1.0e-12);

            maxRichardson = std::max(maxRichardson, richardson);

            // Pacanowski-Philander: mixing collapses as the square of the
            // stratification, which is what lets a thermocline survive.
            const double damping = 1.0 + 5.0 * richardson;

            m_km[slot] = m_options.neutralViscosity / (damping * damping) +
                         m_options.backgroundViscosity;

            m_kz[slot] = m_km[slot] / damping + m_options.backgroundDiffusivity;
        }
    }

    m_diagnostics.maxRichardson = maxRichardson;
}

double InternalModeSolver::bottomDragCoefficient(int64_t column) const
{
    const double depth = m_mesh->columnDepth()[static_cast<size_t>(column)];

    if (depth <= 0.0 || m_options.manningN <= 0.0)
        return 0.0;

    return m_options.gravity * m_options.manningN * m_options.manningN /
           std::pow(depth, 1.0 / 3.0);
}

double InternalModeSolver::stableTimeStep() const
{
    const Mesh::LayeredMesh     &mesh       = *m_mesh;
    const Mesh::HorizontalMesh  &horizontal = mesh.horizontal();
    const State::StateStore     &state      = *m_state;

    const int64_t columns = mesh.columnCount();
    const int     layers  = mesh.layerCount();

    const double spacing = horizontal.minimumCentroidDistance();

    if (spacing <= 0.0)
        return m_options.maxTimeStep;

    double fastest = 0.0;

    // The density contrast is taken over the whole domain, not within each
    // column. A column-local range is zero for a horizontal density front —
    // a lock exchange starts with every column internally uniform — and the
    // estimate would then see no baroclinic wave at all and hand back the
    // maximum step, which is unconditionally unstable from the first step.
    double lightest  = std::numeric_limits<double>::infinity();
    double heaviest  = 0.0;
    double deepest   = 0.0;

    for (int64_t c = 0; c < columns; ++c)
    {
        deepest = std::max(deepest, mesh.columnDepth()[static_cast<size_t>(c)]);

        for (int k = 0; k < layers; ++k)
        {
            const size_t cell = static_cast<size_t>(mesh.cell(c, k));

            lightest = std::min(lightest, state.density[cell]);
            heaviest = std::max(heaviest, state.density[cell]);

            fastest = std::max(fastest, std::hypot(state.u[cell], state.v[cell]));
        }
    }

    const double reduced =
        m_options.gravity * (heaviest - lightest) / m_options.referenceDensity;

    if (reduced > 0.0)
        fastest = std::max(fastest, std::sqrt(reduced * deepest));

    double limit = fastest > 0.0 ? m_options.cflNumber * spacing / fastest
                                 : m_options.maxTimeStep;

    // Vertical advection has its own Courant limit, and it is usually the
    // binding one: sigma layers are metres thick where the horizontal cells
    // are tens of metres wide, so a vertical velocity a hundred times smaller
    // than the horizontal one still crosses its cell first. Omitting this
    // limit lets the horizontal constraint authorise a step that the vertical
    // direction cannot survive.
    for (int64_t c = 0; c < columns; ++c)
    {
        for (int k = 0; k < layers; ++k)
        {
            const size_t cell      = static_cast<size_t>(mesh.cell(c, k));
            const double thickness = mesh.thickness()[cell];

            if (thickness <= 0.0)
                continue;

            const double top =
                std::abs(m_verticalFlux[static_cast<size_t>(mesh.interfaceSlot(c, k))]);
            const double bottom = std::abs(
                m_verticalFlux[static_cast<size_t>(mesh.interfaceSlot(c, k + 1))]);

            const double speed = std::max(top, bottom);

            if (speed > 0.0)
                limit = std::min(limit, m_options.cflNumber * thickness / speed);
        }
    }

    return std::clamp(limit, m_options.minTimeStep, m_options.maxTimeStep);
}

void InternalModeSolver::advanceMomentum(double dt)
{
    const Mesh::LayeredMesh    &mesh       = *m_mesh;
    const Mesh::HorizontalMesh &horizontal = mesh.horizontal();
    State::StateStore          &state      = *m_state;

    const int64_t columns = mesh.columnCount();
    const int     layers  = mesh.layerCount();

    const std::vector<double>  &area            = horizontal.faceArea();
    const std::vector<double>  &edgeLength      = horizontal.edgeLength();
    const std::vector<double>  &edgeDistance    = horizontal.edgeCentroidDistance();
    const std::vector<int64_t> &edgeLeft        = horizontal.edgeLeftFace();
    const std::vector<int64_t> &edgeRight       = horizontal.edgeRightFace();
    const std::vector<int64_t> &faceEdgeOffsets = horizontal.faceEdgeOffsets();
    const std::vector<int64_t> &faceEdges       = horizontal.faceEdges();

    m_workU = state.u;
    m_workV = state.v;
    // Safe to thread: this loop reads state.u/v and writes only m_workU/m_workV,
    // so no iteration can observe another's output.

    // ── Explicit terms: advection, baroclinic pressure, horizontal
    //    diffusion, Coriolis ─────────────────────────────────────────────
    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < columns; ++c)
    {
        const int64_t begin = faceEdgeOffsets[static_cast<size_t>(c)];
        const int64_t end   = faceEdgeOffsets[static_cast<size_t>(c) + 1];

        for (int k = 0; k < layers; ++k)
        {
            const size_t cell   = static_cast<size_t>(mesh.cell(c, k));
            const double volume = mesh.cellVolume()[cell];

            if (volume <= 0.0)
                continue;

            double advectionU = 0.0;
            double advectionV = 0.0;
            double diffusionU = 0.0;
            double diffusionV = 0.0;

            for (int64_t index = begin; index < end; ++index)
            {
                const int64_t e  = faceEdges[static_cast<size_t>(index)];
                const size_t  ei = static_cast<size_t>(e);

                const bool    owned = edgeLeft[ei] == c;
                const int64_t other = owned ? edgeRight[ei] : edgeLeft[ei];

                if (other == Mesh::kNoFace)
                    continue;

                // Outflow from this cell, positive leaving.
                const double flux =
                    (owned ? 1.0 : -1.0) *
                    m_edgeLayerFlux[ei * static_cast<size_t>(layers) +
                                    static_cast<size_t>(k)];

                const size_t neighbour = static_cast<size_t>(mesh.cell(other, k));

                // Advective form: a uniform velocity field is untouched by
                // advection regardless of the divergence, which matters
                // because the layer geometry is breathing under it.
                const double upwindU = flux >= 0.0 ? state.u[cell] : state.u[neighbour];
                const double upwindV = flux >= 0.0 ? state.v[cell] : state.v[neighbour];

                advectionU -= flux * (upwindU - state.u[cell]);
                advectionV -= flux * (upwindV - state.v[cell]);

                if (m_options.horizontalViscosity > 0.0 && edgeDistance[ei] > 0.0)
                {
                    const double faceThickness =
                        0.5 * (mesh.thickness()[cell] + mesh.thickness()[neighbour]);

                    const double weight = m_options.horizontalViscosity *
                                          edgeLength[ei] * faceThickness /
                                          edgeDistance[ei];

                    diffusionU += weight * (state.u[neighbour] - state.u[cell]);
                    diffusionV += weight * (state.v[neighbour] - state.v[cell]);
                }
            }

            // Vertical advection, upwinded across the interfaces.
            const double topFlux =
                m_verticalFlux[static_cast<size_t>(mesh.interfaceSlot(c, k))] *
                area[static_cast<size_t>(c)];
            const double bottomFlux =
                m_verticalFlux[static_cast<size_t>(mesh.interfaceSlot(c, k + 1))] *
                area[static_cast<size_t>(c)];

            if (k > 0)
            {
                const size_t above = static_cast<size_t>(mesh.cell(c, k - 1));
                // Positive topFlux leaves this cell upward.
                const double upwindU = topFlux >= 0.0 ? state.u[cell] : state.u[above];
                const double upwindV = topFlux >= 0.0 ? state.v[cell] : state.v[above];

                advectionU -= topFlux * (upwindU - state.u[cell]);
                advectionV -= topFlux * (upwindV - state.v[cell]);
            }

            if (k + 1 < layers)
            {
                const size_t below = static_cast<size_t>(mesh.cell(c, k + 1));
                // Positive bottomFlux enters this cell from below.
                const double upwindU = bottomFlux >= 0.0 ? state.u[below] : state.u[cell];
                const double upwindV = bottomFlux >= 0.0 ? state.v[below] : state.v[cell];

                advectionU += bottomFlux * (upwindU - state.u[cell]);
                advectionV += bottomFlux * (upwindV - state.v[cell]);
            }

            const double inverseVolume = 1.0 / volume;

            double rateU = (advectionU + diffusionU) * inverseVolume +
                           m_accelerationX[cell];
            double rateV = (advectionV + diffusionV) * inverseVolume +
                           m_accelerationY[cell];

            if (m_options.coriolis != 0.0)
            {
                rateU += m_options.coriolis * state.v[cell];
                rateV -= m_options.coriolis * state.u[cell];
            }

            m_workU[cell] = state.u[cell] + dt * rateU;
            m_workV[cell] = state.v[cell] + dt * rateV;
        }
    }

    // ── Implicit vertical viscosity, one tridiagonal solve per column ─────

    // Surface wind stress per column was filled at the top of the step, so this
    // and the external mode's body acceleration read the same values.

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < columns; ++c)
    {
        // Each thread assembles its own tridiagonal in its own slice.
        const size_t slot        = static_cast<size_t>(threadIndex() * layers);
        double      *lowerBuf    = m_lower.data() + slot;
        double      *diagonalBuf = m_diagonal.data() + slot;
        double      *upperBuf    = m_upper.data() + slot;
        double      *rhsBuf      = m_rhs.data() + slot;
        double      *solutionBuf = m_solution.data() + slot;
        double      *sweepBuf    = m_sweep.data() + slot;

        const double drag = bottomDragCoefficient(c);

        for (int component = 0; component < 2; ++component)
        {
            std::vector<double> &field = component == 0 ? m_workU : m_workV;

            const double surfaceAcceleration =
                (component == 0 ? m_windStressX[static_cast<size_t>(c)]
                                : m_windStressY[static_cast<size_t>(c)]) /
                m_options.referenceDensity;

            for (int k = 0; k < layers; ++k)
            {
                const size_t cell      = static_cast<size_t>(mesh.cell(c, k));
                const double thickness = std::max(mesh.thickness()[cell], 1.0e-12);

                double lower = 0.0;
                double upper = 0.0;

                if (k > 0)
                {
                    const size_t above = static_cast<size_t>(mesh.cell(c, k - 1));
                    const double spacing =
                        0.5 * (mesh.thickness()[above] + mesh.thickness()[cell]);
                    lower = -dt * m_km[static_cast<size_t>(mesh.interfaceSlot(c, k))] /
                            (thickness * std::max(spacing, 1.0e-12));
                }

                if (k + 1 < layers)
                {
                    const size_t below = static_cast<size_t>(mesh.cell(c, k + 1));
                    const double spacing =
                        0.5 * (mesh.thickness()[cell] + mesh.thickness()[below]);
                    upper = -dt * m_km[static_cast<size_t>(mesh.interfaceSlot(c, k + 1))] /
                            (thickness * std::max(spacing, 1.0e-12));
                }

                lowerBuf[k]    = lower;
                upperBuf[k]    = upper;
                diagonalBuf[k] = 1.0 - lower - upper;
                rhsBuf[k]      = field[cell];

                if (k == 0)
                {
                    // Wind stress enters through the surface interface.
                    rhsBuf[0] += dt * surfaceAcceleration / thickness;
                }

                if (k == layers - 1 && drag > 0.0)
                {
                    // Semi-implicit quadratic bed drag on the bottom layer.
                    const double speed = std::hypot(m_workU[cell], m_workV[cell]);
                    diagonalBuf[k] += dt * drag * speed / thickness;
                }
            }

            solveTridiagonal(layers, lowerBuf, diagonalBuf, upperBuf, rhsBuf,
                             solutionBuf, sweepBuf);

            for (int k = 0; k < layers; ++k)
                field[static_cast<size_t>(mesh.cell(c, k))] = solutionBuf[k];
        }
    }

    state.u = m_workU;
    state.v = m_workV;
}

void InternalModeSolver::applyModeAdjustment()
{
    const Mesh::LayeredMesh &mesh  = *m_mesh;
    State::StateStore       &state = *m_state;

    const int64_t columns = mesh.columnCount();
    const int     layers  = mesh.layerCount();

    const std::vector<double> &externalU = m_external.velocityX();
    const std::vector<double> &externalV = m_external.velocityY();

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < columns; ++c)
    {
        double depth = 0.0;
        double meanU = 0.0;
        double meanV = 0.0;

        for (int k = 0; k < layers; ++k)
        {
            const size_t cell      = static_cast<size_t>(mesh.cell(c, k));
            const double thickness = mesh.thickness()[cell];

            depth += thickness;
            meanU += thickness * state.u[cell];
            meanV += thickness * state.v[cell];
        }

        if (depth <= 0.0)
            continue;

        meanU /= depth;
        meanV /= depth;

        // The external mode owns the depth-averaged transport; the internal
        // mode owns only the shear about it.
        const double correctionU = externalU[static_cast<size_t>(c)] - meanU;
        const double correctionV = externalV[static_cast<size_t>(c)] - meanV;

        for (int k = 0; k < layers; ++k)
        {
            const size_t cell = static_cast<size_t>(mesh.cell(c, k));
            state.u[cell] += correctionU;
            state.v[cell] += correctionV;
        }
    }

    state.eta = m_external.surfaceElevation();
    state.ubar = externalU;
    state.vbar = externalV;
}

void InternalModeSolver::computeLayerFluxes()
{
    const Mesh::LayeredMesh    &mesh       = *m_mesh;
    const Mesh::HorizontalMesh &horizontal = mesh.horizontal();
    const State::StateStore    &state      = *m_state;

    const int layers = mesh.layerCount();

    const std::vector<double>  &edgeLength = horizontal.edgeLength();
    const std::vector<double>  &normalX    = horizontal.edgeNormalX();
    const std::vector<double>  &normalY    = horizontal.edgeNormalY();
    const std::vector<int64_t> &edgeLeft   = horizontal.edgeLeftFace();
    const std::vector<int64_t> &edgeRight  = horizontal.edgeRightFace();

    const std::vector<double> &accumulated = m_external.accumulatedEdgeVolume();
    const double               interval    = m_external.accumulationInterval();

    for (int64_t e = 0; e < horizontal.edgeCount(); ++e)
    {
        const size_t  ei    = static_cast<size_t>(e);
        const size_t  base  = ei * static_cast<size_t>(layers);
        const int64_t left  = edgeLeft[ei];
        const int64_t right = edgeRight[ei];

        if (right == Mesh::kNoFace)
        {
            for (int k = 0; k < layers; ++k)
                m_edgeLayerFlux[base + static_cast<size_t>(k)] = 0.0;
            continue;
        }

        double rawTotal       = 0.0;
        double thicknessTotal = 0.0;

        for (int k = 0; k < layers; ++k)
        {
            const size_t leftCell  = static_cast<size_t>(mesh.cell(left, k));
            const size_t rightCell = static_cast<size_t>(mesh.cell(right, k));

            const double faceU = 0.5 * (state.u[leftCell] + state.u[rightCell]);
            const double faceV = 0.5 * (state.v[leftCell] + state.v[rightCell]);

            const double thickness =
                0.5 * (mesh.thickness()[leftCell] + mesh.thickness()[rightCell]);

            const double normal = faceU * normalX[ei] + faceV * normalY[ei];
            const double flux   = normal * thickness * edgeLength[ei];

            m_edgeLayerFlux[base + static_cast<size_t>(k)] = flux;

            rawTotal += flux;
            thicknessTotal += thickness;
        }

        // Flux adjustment: force the layer fluxes to sum to the transport the
        // external mode actually delivered, distributed by layer thickness.
        // Without this the two modes disagree and scalar transport stops
        // conserving mass — the classic mode-splitting failure.
        const double externalFlux =
            interval > 0.0 ? accumulated[ei] / interval : 0.0;

        if (thicknessTotal > 0.0)
        {
            const double deficit = externalFlux - rawTotal;

            for (int k = 0; k < layers; ++k)
            {
                const size_t leftCell  = static_cast<size_t>(mesh.cell(left, k));
                const size_t rightCell = static_cast<size_t>(mesh.cell(right, k));

                const double thickness =
                    0.5 * (mesh.thickness()[leftCell] + mesh.thickness()[rightCell]);

                m_edgeLayerFlux[base + static_cast<size_t>(k)] +=
                    deficit * thickness / thicknessTotal;
            }
        }
    }
}

void InternalModeSolver::diagnoseVerticalFlux(double dt)
{
    const Mesh::LayeredMesh    &mesh       = *m_mesh;
    const Mesh::HorizontalMesh &horizontal = mesh.horizontal();

    const int64_t columns = mesh.columnCount();
    const int     layers  = mesh.layerCount();

    const std::vector<double>  &area            = horizontal.faceArea();
    const std::vector<int64_t> &edgeLeft        = horizontal.edgeLeftFace();
    const std::vector<int64_t> &faceEdgeOffsets = horizontal.faceEdgeOffsets();
    const std::vector<int64_t> &faceEdges       = horizontal.faceEdges();

    m_diagnostics.surfaceFluxResidual = 0.0;
    m_diagnostics.maxVerticalFlux     = 0.0;

    for (int64_t c = 0; c < columns; ++c)
    {
        const int64_t begin       = faceEdgeOffsets[static_cast<size_t>(c)];
        const int64_t end         = faceEdgeOffsets[static_cast<size_t>(c) + 1];
        const double  inverseArea = 1.0 / area[static_cast<size_t>(c)];

        // No flux through the bed: the integration starts there.
        m_verticalFlux[static_cast<size_t>(mesh.interfaceSlot(c, layers))] = 0.0;

        for (int k = layers - 1; k >= 0; --k)
        {
            const size_t cell = static_cast<size_t>(mesh.cell(c, k));

            double horizontalOutflow = 0.0;

            for (int64_t index = begin; index < end; ++index)
            {
                const int64_t e  = faceEdges[static_cast<size_t>(index)];
                const size_t  ei = static_cast<size_t>(e);

                horizontalOutflow +=
                    (edgeLeft[ei] == c ? 1.0 : -1.0) *
                    m_edgeLayerFlux[ei * static_cast<size_t>(layers) +
                                    static_cast<size_t>(k)];
            }

            const double volumeRate = (mesh.cellVolume()[cell] - m_oldVolume[cell]) / dt;

            // Continuity for the layer: what does not leave sideways, does not
            // accumulate, and did not arrive through a source must pass
            // vertically. Because the column totals handed to the external
            // mode are the same numbers summed, the residual at the free
            // surface still cancels exactly.
            const double below =
                m_verticalFlux[static_cast<size_t>(mesh.interfaceSlot(c, k + 1))] *
                area[static_cast<size_t>(c)];

            const double top =
                below - horizontalOutflow - volumeRate + m_cellSource[cell];

            m_verticalFlux[static_cast<size_t>(mesh.interfaceSlot(c, k))] =
                top * inverseArea;
        }

        const double residual =
            m_verticalFlux[static_cast<size_t>(mesh.interfaceSlot(c, 0))];

        m_diagnostics.surfaceFluxResidual =
            std::max(m_diagnostics.surfaceFluxResidual, std::abs(residual));

        for (int k = 0; k <= layers; ++k)
        {
            m_diagnostics.maxVerticalFlux =
                std::max(m_diagnostics.maxVerticalFlux,
                         std::abs(m_verticalFlux[static_cast<size_t>(
                             mesh.interfaceSlot(c, k))]));
        }

        // The free surface is a material surface: force the residual to zero
        // so transport sees an exactly closed budget.
        m_verticalFlux[static_cast<size_t>(mesh.interfaceSlot(c, 0))] = 0.0;
    }
}

void InternalModeSolver::transportScalars(double dt)
{
    const Mesh::LayeredMesh    &mesh       = *m_mesh;
    const Mesh::HorizontalMesh &horizontal = mesh.horizontal();
    State::StateStore          &state      = *m_state;

    const int64_t columns = mesh.columnCount();
    const int     layers  = mesh.layerCount();

    const std::vector<double>  &area            = horizontal.faceArea();
    const std::vector<double>  &edgeLength      = horizontal.edgeLength();
    const std::vector<double>  &edgeDistance    = horizontal.edgeCentroidDistance();
    const std::vector<int64_t> &edgeLeft        = horizontal.edgeLeftFace();
    const std::vector<int64_t> &edgeRight       = horizontal.edgeRightFace();
    const std::vector<int64_t> &faceEdgeOffsets = horizontal.faceEdgeOffsets();
    const std::vector<int64_t> &faceEdges       = horizontal.faceEdges();

    for (int s = 0; s < state.speciesCount(); ++s)
    {
        if (!m_species->definition(s).isTransported())
            continue;

        std::span<double> values = state.speciesBlock(s);

        // An internal transfer's outfall delivers whatever its intake is
        // drawing, so the receiving inflow's composition is set here rather
        // than stated by the caller. It must happen before the values are
        // gathered and before this species is advected: the withdrawal removes
        // mass at the current, pre-advection concentrations, and the arriving
        // water has to be that same water for the transfer to conserve.
        if (m_transfers.count() > 0)
            m_transfers.carry(mesh, s, m_temperatureIndex, values, m_sources);

        // Advection is delegated to the Eulerian engine (MUSCL + Zalesak
        // flux correction). It receives fluxes rather than velocities
        // because the fluxes are what mode splitting already reconciled.
        // An open boundary's ambient composition depends on the placement
        // weights, which distribute() has now computed for this step.
        applyBoundaryComposition(s);

        m_sources.inflowConcentration(mesh, s, m_temperatureIndex, m_sourceValue);

        m_transport.advect(mesh, m_oldVolume, m_edgeLayerFlux, m_verticalFlux,
                           m_cellSource, m_sourceValue, dt, m_options.advectionScheme,
                           values);

        if (m_options.horizontalDiffusivity > 0.0)
            applyHorizontalDiffusion(dt, values);

        // ── Implicit vertical diffusion, one tridiagonal solve per column ─
        FVQual_PARALLEL_FOR
        for (int64_t c = 0; c < columns; ++c)
        {
            // Each thread assembles its own tridiagonal in its own slice.
            const size_t slot            = static_cast<size_t>(threadIndex() * layers);
            double      *lowerBuf    = m_lower.data() + slot;
            double      *diagonalBuf = m_diagonal.data() + slot;
            double      *upperBuf    = m_upper.data() + slot;
            double      *rhsBuf      = m_rhs.data() + slot;
            double      *solutionBuf = m_solution.data() + slot;
            double      *sweepBuf    = m_sweep.data() + slot;

            for (int k = 0; k < layers; ++k)
            {
                const size_t cell      = static_cast<size_t>(mesh.cell(c, k));
                const double thickness = std::max(mesh.thickness()[cell], 1.0e-12);

                double lower = 0.0;
                double upper = 0.0;

                if (k > 0)
                {
                    const size_t above = static_cast<size_t>(mesh.cell(c, k - 1));
                    const double spacing =
                        0.5 * (mesh.thickness()[above] + mesh.thickness()[cell]);
                    lower = -dt * m_kz[static_cast<size_t>(mesh.interfaceSlot(c, k))] /
                            (thickness * std::max(spacing, 1.0e-12));
                }

                if (k + 1 < layers)
                {
                    const size_t below = static_cast<size_t>(mesh.cell(c, k + 1));
                    const double spacing =
                        0.5 * (mesh.thickness()[cell] + mesh.thickness()[below]);
                    upper = -dt * m_kz[static_cast<size_t>(mesh.interfaceSlot(c, k + 1))] /
                            (thickness * std::max(spacing, 1.0e-12));
                }

                lowerBuf[k]    = lower;
                upperBuf[k]    = upper;
                diagonalBuf[k] = 1.0 - lower - upper;
                rhsBuf[k]      = values[cell];
            }

            solveTridiagonal(layers, lowerBuf, diagonalBuf, upperBuf, rhsBuf,
                             solutionBuf, sweepBuf);

            for (int k = 0; k < layers; ++k)
                values[static_cast<size_t>(mesh.cell(c, k))] = solutionBuf[k];
        }
    }
}

void InternalModeSolver::applyHorizontalDiffusion(double dt, std::span<double> values)
{
    const Mesh::LayeredMesh    &mesh       = *m_mesh;
    const Mesh::HorizontalMesh &horizontal = mesh.horizontal();

    const int     layers  = mesh.layerCount();
    const int64_t columns = mesh.columnCount();

    const std::vector<double>  &edgeLength      = horizontal.edgeLength();
    const std::vector<double>  &edgeDistance    = horizontal.edgeCentroidDistance();
    const std::vector<int64_t> &edgeLeft        = horizontal.edgeLeftFace();
    const std::vector<int64_t> &edgeRight       = horizontal.edgeRightFace();
    const std::vector<int64_t> &faceEdgeOffsets = horizontal.faceEdgeOffsets();
    const std::vector<int64_t> &faceEdges       = horizontal.faceEdges();

    for (int64_t i = 0; i < mesh.cellCount(); ++i)
        m_scalarWork[static_cast<size_t>(i)] = values[static_cast<size_t>(i)];

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < columns; ++c)
    {
        const int64_t begin = faceEdgeOffsets[static_cast<size_t>(c)];
        const int64_t end   = faceEdgeOffsets[static_cast<size_t>(c) + 1];

        for (int k = 0; k < layers; ++k)
        {
            const size_t cell   = static_cast<size_t>(mesh.cell(c, k));
            const double volume = mesh.cellVolume()[cell];

            if (volume <= 0.0)
                continue;

            double inflow = 0.0;

            for (int64_t index = begin; index < end; ++index)
            {
                const size_t  ei    = static_cast<size_t>(faceEdges[static_cast<size_t>(index)]);
                const bool    owned = edgeLeft[ei] == c;
                const int64_t other = owned ? edgeRight[ei] : edgeLeft[ei];

                if (other == Mesh::kNoFace || edgeDistance[ei] <= 0.0)
                    continue;

                const size_t neighbour = static_cast<size_t>(mesh.cell(other, k));

                const double faceThickness =
                    0.5 * (mesh.thickness()[cell] + mesh.thickness()[neighbour]);

                inflow += m_options.horizontalDiffusivity * edgeLength[ei] *
                          faceThickness *
                          (m_scalarWork[neighbour] - m_scalarWork[cell]) /
                          edgeDistance[ei];
            }

            values[cell] += dt * inflow / volume;
        }
    }
}

void InternalModeSolver::applySurfaceHeatFlux(double dt)
{
    if (m_temperatureIndex < 0)
        return;

    const Mesh::LayeredMesh &mesh  = *m_mesh;
    State::StateStore       &state = *m_state;

    const double capacity = m_options.referenceDensity * m_options.heatCapacity;
    const int    layers   = mesh.layerCount();

    if (!m_useMeteorology)
    {
        if (m_options.surfaceHeatFlux == 0.0)
            return;

        FVQual_PARALLEL_FOR
        for (int64_t c = 0; c < mesh.columnCount(); ++c)
        {
            const size_t surface   = static_cast<size_t>(mesh.cell(c, 0));
            const double thickness = mesh.thickness()[surface];

            if (thickness <= 0.0)
                continue;

            state.speciesAt(m_temperatureIndex, static_cast<int64_t>(surface)) +=
                dt * m_options.surfaceHeatFlux / (capacity * thickness);
        }

        return;
    }

    // The solar field is marched from the *penetrating* part of the net
    // shortwave, so the irradiance it reports is energy still in transit. The
    // absorbed energy in a cell is then the difference across it, which is what
    // heats that cell — and is the identity LightField's own energy gate pins.
    // Net shortwave per column, because albedo and incident radiation may both
    // vary along the water body (F3.7). Filled into a member rather than a
    // local so the per-step work does not allocate.
    m_incidentShortwave.resize(static_cast<size_t>(mesh.columnCount()));

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < mesh.columnCount(); ++c)
    {
        m_incidentShortwave[static_cast<size_t>(c)] =
            (1.0 - m_meteorology.albedo.at(c)) * m_meteorology.shortwave.at(c);
    }

    Kinetics::LightOptions solar;
    solar.incidentPerColumn     = &m_incidentShortwave;
    solar.surfaceAbsorption     = m_meteorology.base.surfaceAbsorption;
    solar.backgroundAttenuation = m_solarAttenuation;

    m_solar.update(mesh, state, solar);

    const std::vector<double> &interface = m_solar.interfaceIrradiance();

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < mesh.columnCount(); ++c)
    {
        const size_t surface = static_cast<size_t>(mesh.cell(c, 0));

        // The exchange depends on the surface temperature it is trying to set,
        // which is the feedback that gives a water body an equilibrium at all.
        const HeatFluxTerms terms =
            surfaceHeatFlux(m_meteorology.at(c),
                            state.speciesAt(m_temperatureIndex,
                                            static_cast<int64_t>(surface)));

        for (int k = 0; k < layers; ++k)
        {
            const size_t cell      = static_cast<size_t>(mesh.cell(c, k));
            const double thickness = mesh.thickness()[cell];

            if (thickness <= 0.0)
                continue;

            // Absorbed sunlight: what entered the top of this layer and did not
            // leave the bottom.
            double flux =
                interface[static_cast<size_t>(mesh.interfaceSlot(c, k))] -
                interface[static_cast<size_t>(mesh.interfaceSlot(c, k + 1))];

            // Everything else acts at the surface: longwave both ways, latent,
            // sensible, and the shortwave stopped in the first millimetres.
            if (k == 0)
                flux += terms.surfaceNet();

            state.speciesAt(m_temperatureIndex, static_cast<int64_t>(cell)) +=
                dt * flux / (capacity * thickness);
        }

        // Sunlight reaching the bed is absorbed there rather than lost. Adding
        // it to the bottom cell keeps the energy budget closed; the alternative
        // is a model that quietly loses heat through the floor of a shallow,
        // clear reservoir.
        const double toBed = interface[static_cast<size_t>(mesh.interfaceSlot(c, layers))];

        const size_t bed          = static_cast<size_t>(mesh.cell(c, layers - 1));
        const double bedThickness = mesh.thickness()[bed];

        if (toBed > 0.0 && bedThickness > 0.0)
        {
            state.speciesAt(m_temperatureIndex, static_cast<int64_t>(bed)) +=
                dt * toBed / (capacity * bedThickness);
        }
    }
}

void InternalModeSolver::updateCellDepth()
{
    const Mesh::LayeredMesh &mesh   = *m_mesh;
    const int64_t            columns = mesh.columnCount();
    const int                layers  = mesh.layerCount();

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < columns; ++c)
    {
        const double surface = mesh.z(c, 0);

        for (int k = 0; k < layers; ++k)
        {
            // Depth of the cell centre below the free surface, so it is
            // positive downward and moves with the surface rather than being
            // measured from a fixed datum.
            m_cellDepth[static_cast<size_t>(mesh.cell(c, k))] =
                surface - 0.5 * (mesh.z(c, k) + mesh.z(c, k + 1));
        }
    }
}

bool InternalModeSolver::bindOpenBoundaries(Forcing::Placement placement,
                                            std::string &message)
{
    if (m_mesh == nullptr)
    {
        message = "internal mode has not been initialized";
        return false;
    }

    if (!m_boundaryLinks.empty())
    {
        message = "open boundaries are already bound";
        return false;
    }

    const std::vector<int> &edgeBoundary = m_external.edgeBoundary();

    if (edgeBoundary.empty() || !m_external.hasOpenBoundaries())
    {
        message = "no open boundaries have been added to the external mode";
        return false;
    }

    const Mesh::HorizontalMesh &horizontal = m_mesh->horizontal();
    const std::vector<int64_t> &edgeLeft   = horizontal.edgeLeftFace();
    const std::vector<int64_t> &edgeRight  = horizontal.edgeRightFace();

    // One source per (boundary, column). A boundary may span several columns
    // and a column may sit on more than one boundary, so neither alone is the
    // right key.
    for (int64_t e = 0; e < horizontal.edgeCount(); ++e)
    {
        const int boundary = edgeBoundary[static_cast<size_t>(e)];

        if (boundary < 0)
            continue;

        // A boundary edge has no right-hand face by construction; the interior
        // column is always the left one.
        const int64_t column =
            edgeRight[static_cast<size_t>(e)] == Mesh::kNoFace
                ? edgeLeft[static_cast<size_t>(e)]
                : edgeRight[static_cast<size_t>(e)];

        auto link = std::find_if(m_boundaryLinks.begin(), m_boundaryLinks.end(),
                                 [&](const BoundaryLink &l) {
                                     return l.boundary == boundary && l.column == column;
                                 });

        if (link == m_boundaryLinks.end())
        {
            Forcing::Source source;
            source.column              = column;
            source.discharge           = 0.0;
            source.placement           = placement;
            source.externallyAccounted = true;

            const int index = m_sources.add(source, message);

            if (index < 0)
                return false;

            BoundaryLink created;
            created.boundary        = boundary;
            created.column          = column;
            created.source          = index;
            created.inflowPlacement = placement;

            m_boundaryLinks.push_back(std::move(created));
            link = std::prev(m_boundaryLinks.end());
        }

        link->edges.push_back(e);
    }

    if (m_boundaryLinks.empty())
    {
        message = "the external mode reports open boundaries but no edge carries one";
        return false;
    }

    // NaN marks "ambient", which is the default: an unconfigured boundary
    // carries the receiving column's own values and so is transparent to
    // constituents rather than a source of distilled water.
    m_boundaryValue.assign(static_cast<size_t>(m_external.boundaryCount()) *
                               static_cast<size_t>(m_species->count()),
                           std::numeric_limits<double>::quiet_NaN());

    message.clear();
    return true;
}

bool InternalModeSolver::setBoundaryValue(int boundary, const std::string &speciesId,
                                          double value, std::string &message)
{
    if (m_boundaryLinks.empty())
    {
        message = "open boundaries have not been bound; call bindOpenBoundaries()";
        return false;
    }

    if (boundary < 0 || boundary >= m_external.boundaryCount())
    {
        message = "boundary index " + std::to_string(boundary) + " is out of range";
        return false;
    }

    const int species = speciesIndex(speciesId);

    if (species < 0)
    {
        message = "species '" + speciesId + "' is not registered";
        return false;
    }

    m_boundaryValue[static_cast<size_t>(boundary) *
                        static_cast<size_t>(m_species->count()) +
                    static_cast<size_t>(species)] = value;

    message.clear();
    return true;
}

bool InternalModeSolver::clearBoundaryValue(int boundary, const std::string &speciesId,
                                            std::string &message)
{
    return setBoundaryValue(boundary, speciesId,
                            std::numeric_limits<double>::quiet_NaN(), message);
}

double InternalModeSolver::boundaryVolume(int boundary) const
{
    double total = 0.0;

    for (const BoundaryLink &link : m_boundaryLinks)
    {
        if (link.boundary == boundary)
            total += link.volume;
    }

    return total;
}

void InternalModeSolver::refreshBoundarySources(double dt)
{
    if (m_boundaryLinks.empty())
        return;

    const std::vector<double> &accumulated = m_external.accumulatedEdgeVolume();

    std::string message;

    for (BoundaryLink &link : m_boundaryLinks)
    {
        double volume = 0.0;

        for (int64_t e : link.edges)
        {
            // The edge volume is signed left-to-right and a boundary edge's
            // normal points out of the domain, so entering water is negative
            // there and positive here.
            volume -= accumulated[static_cast<size_t>(e)];
        }

        link.volume = volume;

        // The placement follows the sign of the flow, and has to: the
        // configured placement describes where *arriving* water belongs, and a
        // density-matched search has no meaning for water leaving. Measured
        // rather than reasoned — a density placement applied to the outflow
        // left the mode-consistency residual at 1.2e-2 where uniform gives
        // 1e-15, because the withdrawal was being placed somewhere the volume
        // change was not.
        //
        // Outflow leaves uniformly over the column. A boundary has no outlet
        // elevation or width, so the selective-withdrawal treatment that a real
        // port gets is not available to it; if that matters for a particular
        // opening, model it as a structure rather than a boundary.
        (void)m_sources.setPlacement(link.source,
                                     volume > 0.0 ? link.inflowPlacement
                                                  : Forcing::Placement::Uniform,
                                     message);

        // The external mode moved this volume over the interval it accumulated
        // across, which is the internal step; the source is a rate.
        (void)m_sources.setDischarge(link.source, volume / dt, message);
    }
}

void InternalModeSolver::applyBoundaryComposition(int species)
{
    if (m_boundaryLinks.empty())
        return;

    const int count = m_species->count();

    std::string message;

    for (const BoundaryLink &link : m_boundaryLinks)
    {
        if (link.volume <= 0.0)
            continue; // A withdrawal carries ambient values by construction.

        const double prescribed =
            m_boundaryValue[static_cast<size_t>(link.boundary) *
                                static_cast<size_t>(count) +
                            static_cast<size_t>(species)];

        // Ambient: the inflow arrives carrying exactly what the layers it is
        // being placed into already hold, so it changes the volume and not the
        // chemistry. That has to be the *placement-weighted* mean rather than
        // the column mean — an inflow spread by density into the hypolimnion
        // and given the whole column's average would import epilimnion water to
        // the bed, which is an artefact, not a neutral boundary.
        //
        // It also has to be computed here rather than when the discharge is
        // set: the weights come from distribute(), which has not run for this
        // step's discharges at that point. Reading them early returned zero
        // weights and so a zero ambient value — an open boundary quietly
        // importing distilled water, which is exactly the failure the ambient
        // default exists to avoid.
        const double value =
            std::isnan(prescribed)
                ? m_sources.weightedMean(*m_mesh, link.source,
                                         m_state->speciesBlock(species))
                : prescribed;

        (void)m_sources.setInflowValue(link.source, species, m_temperatureIndex, value,
                                       message);
    }
}

int InternalModeSolver::speciesIndex(const std::string &id) const noexcept
{
    if (m_species == nullptr)
        return -1;

    return m_species->indexOf(id);
}

bool InternalModeSolver::compileReactions(std::string &message)
{
    if (m_species == nullptr)
    {
        message = "internal mode has not been initialized";
        return false;
    }

    return m_reactions.compile(*m_species, message);
}

bool InternalModeSolver::step(double dt, std::string &message)
{
    if (m_mesh == nullptr)
    {
        message = "internal mode has not been initialized";
        return false;
    }

    if (!(dt > 0.0))
    {
        message = "internal time step must be positive";
        return false;
    }

    // Checked here rather than at the point of use, so a misconfigured network
    // refuses the step instead of failing halfway through one and leaving the
    // state transported but unreacted.
    if (!m_reactions.empty() && !m_reactions.compiled())
    {
        message = "reaction network has rate laws but has not been compiled; "
                  "call compileReactions()";
        return false;
    }

    Mesh::LayeredMesh &mesh = *m_mesh;

    ScopedTimer stepTimer(m_profile.total);

    // 0. Surface wind stress, per column, before anything reads it.
    //
    // Filled here rather than where it is first used because two consumers need
    // it at different points in the step — the external mode's body acceleration
    // and the internal mode's surface boundary condition — and they must see the
    // *same* stress. Computing it inside the later consumer would hand the
    // earlier one the previous step's wind, and a force lagged across a
    // subcycling interval is how Phase 3's bug 1 turned a millimetre of surface
    // tilt into several metres.
    //
    // It depends only on meteorology and geometry, never on the evolving state,
    // so there is nothing to gain by deferring it.
    {
        const int64_t columns = mesh.columnCount();

        m_windStressX.assign(static_cast<size_t>(columns), m_options.windStressX);
        m_windStressY.assign(static_cast<size_t>(columns), m_options.windStressY);

        if (m_useMeteorology)
        {
            FVQual_PARALLEL_FOR
            for (int64_t c = 0; c < columns; ++c)
            {
                const WindStress stress = windStress(m_meteorology.at(c));

                m_windStressX[static_cast<size_t>(c)] = stress.x;
                m_windStressY[static_cast<size_t>(c)] = stress.y;
            }
        }
    }

    // 1. Density, then the mixing it controls, then the pressure it exerts.
    {
        ScopedTimer timer(m_profile.density);
        updateDensity();
        updateVerticalMixing(dt);
    }

    {
        ScopedTimer timer(m_profile.pressure);
        m_pressure->baroclinicAcceleration(mesh, m_state->density,
                                           m_options.referenceDensity, m_options.gravity,
                                           m_accelerationX, m_accelerationY);
    }

    // 2. Geometry before the surface moves — transport needs both volumes.
    m_oldVolume = mesh.cellVolume();

    // 3. Hand the external mode its share of the baroclinic forcing, depth
    //    averaged. Skipping this decouples the modes: the surface would stay
    //    flat while the interior accelerates, and the mode adjustment in
    //    step 6 would then erase the barotropic response instead of
    //    reconciling it.
    const int layers = mesh.layerCount();

    m_columnAccelerationX.assign(static_cast<size_t>(mesh.columnCount()), 0.0);
    m_columnAccelerationY.assign(static_cast<size_t>(mesh.columnCount()), 0.0);

    for (int64_t c = 0; c < mesh.columnCount(); ++c)
    {
        double depth = 0.0;
        double sumX  = 0.0;
        double sumY  = 0.0;

        for (int k = 0; k < layers; ++k)
        {
            const size_t cell      = static_cast<size_t>(mesh.cell(c, k));
            const double thickness = mesh.thickness()[cell];

            depth += thickness;
            sumX += thickness * m_accelerationX[cell];
            sumY += thickness * m_accelerationY[cell];
        }

        if (depth > 0.0)
        {
            m_columnAccelerationX[static_cast<size_t>(c)] = sumX / depth;
            m_columnAccelerationY[static_cast<size_t>(c)] = sumY / depth;

        }
    }

    if (!m_external.setBodyAcceleration(m_columnAccelerationX, m_columnAccelerationY,
                                        message))
    {
        return false;
    }

    // Wind reaches the external mode as a surface stress, not folded into the
    // body force. See BarotropicSolver::setWindStress for why the two are not
    // interchangeable.
    if (!m_windStressX.empty() &&
        !m_external.setWindStress(m_windStressX, m_windStressY, message))
    {
        return false;
    }

    // Forcing: transfers, outlet groups, structures, and source placement. In
    // its own scope so the timer records even on the early-return paths below.
    {
        ScopedTimer timer(m_profile.forcing);

    // 3a--. Internal transfers move water between columns; their two halves are
    //       written together so they cannot drift apart into a net source.
    if (m_transfers.count() > 0 && !m_transfers.updateDischarges(m_sources, message))
        return false;

    // 3a-. Operate the multi-port outlets first: the group decides how the
    //      scheduled release is split between ports, then the ports rate
    //      themselves against that schedule. The order matters — a port's
    //      capacity limit has to be applied to a demand that already exists.
    if (m_groups.count() > 0 && m_temperatureIndex >= 0 &&
        !m_groups.update(mesh, m_sources, m_state->density,
                         m_state->speciesBlock(m_temperatureIndex), m_options.gravity,
                         m_options.referenceDensity, m_structures, message))
    {
        return false;
    }

    // 3a. Let the outlet works read the level they are about to release
    //     against. A weir whose reservoir has drawn down passes less than it
    //     did last step, and that feedback is the difference between a model
    //     that drains to its crest and stops and one that empties itself.
    //
    //     This is the opening rating only. The feedback itself lives in the
    //     external subcycle, where the head actually moves — see the updater
    //     installed in initialize() and StructureSet::update().
    if (m_structures.count() > 0 &&
        !m_structures.update(m_external.surfaceElevation(), m_sources, message))
    {
        return false;
    }

    // 3b. Place the inflows and withdrawals against the *current* density
    //     field, then give the external mode exactly the column totals that
    //     were distributed over the cells. The two must agree to the bit: the
    //     free-surface continuity residual is the difference between what the
    //     surface thinks happened and what the layers think happened, and a
    //     source counted once here and differently there would show up there
    //     first.
    m_sources.distribute(mesh, m_state->density, m_options.gravity,
                         m_options.referenceDensity, m_cellSource, m_columnSource);

    if (!m_external.setColumnSource(m_columnSource, message))
        return false;
    }

    // 4. External mode, accumulating the transport it delivers.
    m_structureVolume.assign(static_cast<size_t>(m_structures.count()), 0.0);

    m_external.resetFluxAccumulation();

    {
        ScopedTimer timer(m_profile.external);

        if (!m_external.advance(dt, message))
            return false;
    }

    m_diagnostics.externalSteps = m_external.diagnostics().steps;

    // 4b. The subcycle re-rated every structure as the level fell, so the
    //     release that actually happened is not the rating at either end of the
    //     step but the mean over it. Transport has to be handed that mean:
    //     giving it the end-of-step rating would move a different volume of
    //     water than the free surface just did, and the two modes would
    //     disagree about a spillway by exactly the amount the head changed.
    const bool restated = (m_structures.count() > 0 && dt > 0.0) ||
                          !m_boundaryLinks.empty();

    if (m_structures.count() > 0 && dt > 0.0)
    {
        for (int s = 0; s < m_structures.count(); ++s)
        {
            const double mean = m_structureVolume[static_cast<size_t>(s)] / dt;

            if (!m_sources.setDischarge(m_structures.sourceIndex(s), mean, message))
                return false;
        }
    }

    // 4c. An open boundary is known only after the fact, for the same reason a
    //     structure is: how much water crossed it is what the subcycle worked
    //     out, not something stated in advance. So its source is set here from
    //     the volume the external mode actually passed, and redistributed with
    //     the rest.
    if (!m_boundaryLinks.empty() && dt > 0.0)
        refreshBoundarySources(dt);

    if (restated)
    {
        m_sources.distribute(mesh, m_state->density, m_options.gravity,
                             m_options.referenceDensity, m_cellSource, m_columnSource);
    }

    // 5. Internal momentum on the old geometry (advection uses the previous
    //    step's fluxes, which is what makes the scheme explicit).
    {
        ScopedTimer timer(m_profile.momentum);
        advanceMomentum(dt);
    }

    // 6. New layer geometry from the new free surface.
    {
        ScopedTimer timer(m_profile.geometry);

        if (!mesh.updateSurface(m_external.surfaceElevation(), message))
            return false;

        // 7. Reconcile the two modes.
        applyModeAdjustment();
    }

    // 7b. Build fluxes consistent with both.
    {
        ScopedTimer timer(m_profile.fluxes);
        computeLayerFluxes();
        diagnoseVerticalFlux(dt);
    }

    // 8. Scalars ride the reconciled fluxes.
    {
        ScopedTimer timer(m_profile.transport);
        transportScalars(dt);
        applySurfaceHeatFlux(dt);
    }

    // 8b. Particulates sink. Before the reactions rather than after, so a
    // benthic rate law acts on material that arrived this step.
    if (m_settling.any())
    {
        ScopedTimer timer(m_profile.settling);
        m_settling.apply(mesh, *m_state, dt);
    }

    // 9. Reactions act on the transported field: first-order Lie splitting.
    //
    // After transport rather than before, and before the density update rather
    // than after, so a reaction that changes temperature is reflected in the
    // density the *next* step's pressure gradient is built from.
    if (!m_reactions.empty())
    {
        ScopedTimer timer(m_profile.reactions);

        updateCellDepth();

        // Light is marched from the just-transported field, so self-shading
        // reflects the constituents the reactions are about to act on. Skipped
        // entirely when nothing reads it: the march costs a pass over every
        // cell with a species loop inside, unlike the depth field above, which
        // is one subtraction and not worth branching on.
        if (m_reactions.usesLight())
            m_light.update(mesh, *m_state, m_lightOptions);

        Kinetics::DerivedFields fields;
        fields.light     = m_light.meanIrradiance();
        fields.depth     = m_cellDepth;
        fields.thickness = mesh.thickness();

        if (!m_reactions.integrate(*m_species, *m_state, m_state->cellCount(), dt,
                                   message, fields, m_cellPosition))
        {
            return false;
        }
    }

    updateDensity();

    m_time += dt;

    m_diagnostics.timeStep    = dt;
    m_diagnostics.totalVolume = mesh.totalVolume();

    m_diagnostics.maxHorizontalVelocity = 0.0;
    for (int64_t i = 0; i < m_state->cellCount(); ++i)
    {
        m_diagnostics.maxHorizontalVelocity =
            std::max(m_diagnostics.maxHorizontalVelocity,
                     std::hypot(m_state->u[static_cast<size_t>(i)],
                                m_state->v[static_cast<size_t>(i)]));
    }

    message.clear();
    return true;
}

bool InternalModeSolver::advance(double duration, std::string &message)
{
    if (m_mesh == nullptr)
    {
        message = "internal mode has not been initialized";
        return false;
    }

    if (duration < 0.0)
    {
        message = "advance duration must be non-negative";
        return false;
    }

    const double target = m_time + duration;
    int64_t      steps  = 0;

    while (m_time < target - 1.0e-12)
    {
        const double dt = std::min(stableTimeStep(), target - m_time);

        if (dt <= 0.0)
            break;

        if (!step(dt, message))
            return false;

        ++steps;
    }

    m_time              = target;
    m_diagnostics.steps = steps;

    message.clear();
    return true;
}

double InternalModeSolver::totalSpeciesMass(int species) const
{
    return m_state->totalMass(species, *m_mesh);
}

} // namespace FVQual::Hydro
