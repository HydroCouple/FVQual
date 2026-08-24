// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/hydro/barotropicsolver.h"

#include "fvqual/parallel.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace FVQual::Hydro {

namespace {

/*!
 * \brief Conserved state in the face-normal frame.
 */
struct NormalState
{
    double h  = 0.0; //!< Depth.
    double un = 0.0; //!< Normal velocity.
    double ut = 0.0; //!< Tangential velocity.
};

/*!
 * \brief Flux in the face-normal frame.
 */
struct NormalFlux
{
    double mass       = 0.0;
    double normal     = 0.0;
    double tangential = 0.0;
};

/*!
 * \brief Physical shallow-water flux of a state in the normal frame.
 */
NormalFlux physicalFlux(const NormalState &s, double gravity)
{
    NormalFlux f;
    f.mass       = s.h * s.un;
    f.normal     = s.h * s.un * s.un + 0.5 * gravity * s.h * s.h;
    f.tangential = s.h * s.un * s.ut;
    return f;
}

/*!
 * \brief HLLC approximate Riemann solver for the shallow-water equations in
 *        the face-normal frame.
 *
 * Wave speeds follow the two-rarefaction estimate of Toro, with the dry-front
 * speeds substituted when either side is dry: a wetting front advances at
 * u ± 2c, not u ± c, and using the wet estimate there is the classic way to
 * lose positivity at the shoreline.
 *
 * The contact wave carries tangential momentum: it is simply upwinded across
 * the contact, which is what distinguishes HLLC from HLL and what keeps shear
 * layers — wind-driven gyres, jet-like inflows — from smearing.
 */
NormalFlux hllc(const NormalState &left, const NormalState &right, double gravity,
                double dryDepth)
{
    const bool leftDry  = left.h <= dryDepth;
    const bool rightDry = right.h <= dryDepth;

    if (leftDry && rightDry)
        return NormalFlux{};

    const double cl = std::sqrt(gravity * std::max(left.h, 0.0));
    const double cr = std::sqrt(gravity * std::max(right.h, 0.0));

    double sl = 0.0;
    double sr = 0.0;

    if (leftDry)
    {
        // Dry-left: the front runs into the dry bed at the left.
        sl = right.un - 2.0 * cr;
        sr = right.un + cr;
    }
    else if (rightDry)
    {
        sl = left.un - cl;
        sr = left.un + 2.0 * cl;
    }
    else
    {
        // Two-rarefaction estimate of the star depth, then Einfeldt-style
        // speeds corrected for a shock on either side.
        const double hStar =
            std::pow(0.5 * (cl + cr) + 0.25 * (left.un - right.un), 2.0) / gravity;

        const double ql = hStar > left.h
                              ? std::sqrt(0.5 * (hStar + left.h) * hStar / (left.h * left.h))
                              : 1.0;
        const double qr =
            hStar > right.h
                ? std::sqrt(0.5 * (hStar + right.h) * hStar / (right.h * right.h))
                : 1.0;

        sl = left.un - cl * ql;
        sr = right.un + cr * qr;
    }

    if (sl >= 0.0)
        return physicalFlux(left, gravity);

    if (sr <= 0.0)
        return physicalFlux(right, gravity);

    const NormalFlux fl = physicalFlux(left, gravity);
    const NormalFlux fr = physicalFlux(right, gravity);

    const double denominator = sr - sl;

    NormalFlux flux;
    flux.mass = (sr * fl.mass - sl * fr.mass + sl * sr * (right.h - left.h)) / denominator;
    flux.normal = (sr * fl.normal - sl * fr.normal +
                   sl * sr * (right.h * right.un - left.h * left.un)) /
                  denominator;

    // Contact speed; tangential momentum is upwinded across it.
    const double contact =
        (sl * right.h * (right.un - sr) - sr * left.h * (left.un - sl)) /
        (right.h * (right.un - sr) - left.h * (left.un - sl));

    flux.tangential = flux.mass * (contact >= 0.0 ? left.ut : right.ut);

    return flux;
}

} // namespace

bool BarotropicSolver::initialize(const Mesh::HorizontalMesh &mesh,
                                  const BarotropicOptions &options, std::string &message)
{
    if (mesh.faceCount() <= 0)
    {
        message = "barotropic solver needs a mesh with at least one face";
        return false;
    }

    if (options.gravity <= 0.0 || options.waterDensity <= 0.0)
    {
        message = "gravity and water density must be positive";
        return false;
    }

    if (options.cflNumber <= 0.0 || options.cflNumber > 1.0)
    {
        message = "the Courant number must lie in (0, 1]";
        return false;
    }

    if (options.dryDepth < 0.0)
    {
        message = "dryDepth must be non-negative";
        return false;
    }

    if (options.minTimeStep <= 0.0 || options.maxTimeStep < options.minTimeStep)
    {
        message = "time-step bounds must satisfy 0 < minTimeStep <= maxTimeStep";
        return false;
    }

    m_mesh        = &mesh;
    m_options     = options;
    m_columnCount = mesh.faceCount();
    m_time        = 0.0;
    m_diagnostics = BarotropicDiagnostics{};

    const size_t columns = static_cast<size_t>(m_columnCount);

    m_h.assign(columns, 0.0);
    m_hu.assign(columns, 0.0);
    m_hv.assign(columns, 0.0);

    m_eta.assign(columns, 0.0);
    m_u.assign(columns, 0.0);
    m_v.assign(columns, 0.0);

    m_dh.assign(columns, 0.0);
    m_dhu.assign(columns, 0.0);
    m_dhv.assign(columns, 0.0);

    m_manning.assign(columns, options.manningN);

    const size_t edges = static_cast<size_t>(mesh.edgeCount());

    m_edgeFluxMass.assign(edges, 0.0);
    m_edgeFluxX.assign(edges, 0.0);
    m_edgeFluxY.assign(edges, 0.0);
    m_edgeCorrectionL.assign(edges, 0.0);
    m_edgeCorrectionR.assign(edges, 0.0);

    m_edgeVolumeFlux.assign(static_cast<size_t>(mesh.edgeCount()), 0.0);
    m_edgeVolumeAccumulated.assign(static_cast<size_t>(mesh.edgeCount()), 0.0);
    m_accumulationInterval = 0.0;

    reconstructPrimitives();

    message.clear();
    return true;
}

void BarotropicSolver::setUniformSurface(double elevation)
{
    const std::vector<double> &bed = m_mesh->faceBed();

    for (int64_t c = 0; c < m_columnCount; ++c)
    {
        const size_t i = static_cast<size_t>(c);
        m_h[i]         = std::max(0.0, elevation - bed[i]);
        m_hu[i]        = 0.0;
        m_hv[i]        = 0.0;
    }

    reconstructPrimitives();
}

bool BarotropicSolver::setSurface(std::span<const double> elevation, std::string &message)
{
    if (static_cast<int64_t>(elevation.size()) != m_columnCount)
    {
        message = "surface elevation array has " + std::to_string(elevation.size()) +
                  " entries; expected " + std::to_string(m_columnCount);
        return false;
    }

    const std::vector<double> &bed = m_mesh->faceBed();

    for (int64_t c = 0; c < m_columnCount; ++c)
    {
        const size_t i = static_cast<size_t>(c);
        m_h[i]         = std::max(0.0, elevation[i] - bed[i]);
        m_hu[i]        = 0.0;
        m_hv[i]        = 0.0;
    }

    reconstructPrimitives();

    message.clear();
    return true;
}

bool BarotropicSolver::setVelocity(std::span<const double> u, std::span<const double> v,
                                   std::string &message)
{
    if (static_cast<int64_t>(u.size()) != m_columnCount ||
        static_cast<int64_t>(v.size()) != m_columnCount)
    {
        message = "velocity arrays must have one entry per column";
        return false;
    }

    for (int64_t c = 0; c < m_columnCount; ++c)
    {
        const size_t i = static_cast<size_t>(c);
        const bool   wet = m_h[i] > m_options.dryDepth;

        m_hu[i] = wet ? m_h[i] * u[i] : 0.0;
        m_hv[i] = wet ? m_h[i] * v[i] : 0.0;
    }

    reconstructPrimitives();

    message.clear();
    return true;
}

bool BarotropicSolver::setManning(std::span<const double> manning, std::string &message)
{
    if (static_cast<int64_t>(manning.size()) != m_columnCount)
    {
        message = "manning array must have one entry per column";
        return false;
    }

    for (double value : manning)
    {
        if (value < 0.0)
        {
            message = "manning roughness must be non-negative";
            return false;
        }
    }

    m_manning.assign(manning.begin(), manning.end());

    message.clear();
    return true;
}

bool BarotropicSolver::setBodyAcceleration(std::span<const double> ax,
                                           std::span<const double> ay,
                                           std::string &message)
{
    if (static_cast<int64_t>(ax.size()) != m_columnCount ||
        static_cast<int64_t>(ay.size()) != m_columnCount)
    {
        message = "body acceleration arrays must have one entry per column";
        return false;
    }

    m_bodyAccelerationX.assign(ax.begin(), ax.end());
    m_bodyAccelerationY.assign(ay.begin(), ay.end());

    message.clear();
    return true;
}

void BarotropicSolver::clearBodyAcceleration() noexcept
{
    m_bodyAccelerationX.clear();
    m_bodyAccelerationY.clear();
}

int BarotropicSolver::addBoundary(const BoundaryCondition &condition,
                                  std::span<const int64_t> edges, std::string &message)
{
    if (m_mesh == nullptr)
    {
        message = "barotropic solver has not been initialized";
        return -1;
    }

    if (edges.empty())
    {
        message = "boundary '" + condition.id + "' has no edges";
        return -1;
    }

    if (condition.type == BoundaryType::Inflow && condition.specificDischarge < 0.0)
    {
        message = "an inflow's specific discharge is positive into the domain; "
                  "use a Stage or Transmissive boundary to let water leave";
        return -1;
    }

    const int64_t edgeCount = m_mesh->edgeCount();

    if (m_edgeBoundary.empty())
        m_edgeBoundary.assign(static_cast<size_t>(edgeCount), -1);

    for (int64_t e : edges)
    {
        if (e < 0 || e >= edgeCount)
        {
            message = "boundary '" + condition.id + "' references edge " +
                      std::to_string(e) + ", which this mesh does not have";
            return -1;
        }

        // Only an unpaired edge is on the boundary. Opening an interior edge
        // would tear a hole in the middle of the domain.
        if (m_mesh->edgeRightFace()[static_cast<size_t>(e)] != Mesh::kNoFace)
        {
            message = "boundary '" + condition.id + "' references edge " +
                      std::to_string(e) +
                      ", which is interior; only unpaired edges lie on the boundary";
            return -1;
        }

        if (m_edgeBoundary[static_cast<size_t>(e)] >= 0)
        {
            message = "edge " + std::to_string(e) +
                      " already carries a boundary condition";
            return -1;
        }
    }

    m_boundaries.push_back(condition);
    m_boundaryFlux.push_back(0.0);

    const int index = static_cast<int>(m_boundaries.size()) - 1;

    for (int64_t e : edges)
        m_edgeBoundary[static_cast<size_t>(e)] = index;

    message.clear();
    return index;
}

bool BarotropicSolver::setBoundaryDischarge(int index, double specificDischarge,
                                            std::string &message)
{
    if (index < 0 || index >= boundaryCount())
    {
        message = "boundary index is out of range";
        return false;
    }

    if (specificDischarge < 0.0)
    {
        message = "an inflow's specific discharge is positive into the domain";
        return false;
    }

    m_boundaries[static_cast<size_t>(index)].specificDischarge = specificDischarge;

    message.clear();
    return true;
}

bool BarotropicSolver::setBoundaryStage(int index, double stage, std::string &message)
{
    if (index < 0 || index >= boundaryCount())
    {
        message = "boundary index is out of range";
        return false;
    }

    m_boundaries[static_cast<size_t>(index)].stage = stage;

    message.clear();
    return true;
}

void BarotropicSolver::ghostState(const BoundaryCondition &condition, double hL,
                                  double unL, double utL, double zL, double &hR,
                                  double &unR, double &utR) const
{
    const double gravity = m_options.gravity;
    const double celerity = std::sqrt(gravity * std::max(hL, 0.0));

    // The outgoing Riemann invariant, carried to the boundary from inside along
    // the characteristic travelling at u + c. This is the piece the exterior
    // does not get to choose.
    const double invariant = unL + 2.0 * celerity;

    utR = utL;

    switch (condition.type)
    {
    case BoundaryType::Wall:
        hR  = hL;
        unR = -unL;
        return;

    case BoundaryType::Transmissive:
        hR  = hL;
        unR = unL;
        return;

    case BoundaryType::Stage:
    {
        // Supercritical outflow: both characteristics leave, so nothing
        // outside can reach in and the stage is simply not a boundary
        // condition here. Imposing it anyway is the classic way to
        // manufacture a spurious upstream-propagating shock.
        if (unL >= celerity)
        {
            hR  = hL;
            unR = unL;
            return;
        }

        const double depth = condition.stage - zL;

        if (depth <= m_options.dryDepth)
        {
            hR  = hL;
            unR = -unL;
            return;
        }

        hR  = depth;
        unR = invariant - 2.0 * std::sqrt(gravity * hR);
        return;
    }

    case BoundaryType::Inflow:
    {
        const double q = condition.specificDischarge;

        if (q <= 0.0)
        {
            hR  = hL;
            unR = -unL;
            return;
        }

        // A stated depth means the caller is describing a supercritical
        // inflow, where both characteristics enter and the invariant has
        // nothing to contribute.
        if (condition.depth > 0.0)
        {
            hR  = condition.depth;
            unR = -q / condition.depth;
            return;
        }

        // Subcritical: the depth is whatever makes the imposed discharge
        // consistent with the invariant the interior is sending out,
        //
        //     -q/h + 2 sqrt(g h) = R .
        //
        // The left side increases monotonically in h from -inf to +inf, so the
        // root is unique and bisection cannot fail. Newton could, near h -> 0.
        const auto residual = [&](double h) {
            return -q / h + 2.0 * std::sqrt(gravity * h) - invariant;
        };

        double low  = 1.0e-6;
        double high = std::max(hL, 1.0);

        while (residual(high) < 0.0 && high < 1.0e6)
            high *= 2.0;

        for (int iteration = 0; iteration < 60; ++iteration)
        {
            const double mid = 0.5 * (low + high);

            if (residual(mid) < 0.0)
                low = mid;
            else
                high = mid;
        }

        hR  = 0.5 * (low + high);
        unR = -q / hR;
        return;
    }
    }
}

bool BarotropicSolver::setColumnSource(std::span<const double> discharge,
                                       std::string &message)
{
    if (static_cast<int64_t>(discharge.size()) != m_columnCount)
    {
        message = "column source array must have one entry per column";
        return false;
    }

    m_columnSource.assign(discharge.begin(), discharge.end());

    message.clear();
    return true;
}

void BarotropicSolver::clearColumnSource() noexcept
{
    m_columnSource.clear();
}

void BarotropicSolver::reconstructPrimitives()
{
    const std::vector<double> &bed = m_mesh->faceBed();

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < m_columnCount; ++c)
    {
        const size_t i = static_cast<size_t>(c);

        // Positivity: the Riemann solver can only guarantee non-negative
        // depth up to round-off, so clamp rather than let a -1e-18 depth
        // propagate into a square root.
        if (m_h[i] < 0.0)
            m_h[i] = 0.0;

        if (m_h[i] > m_options.dryDepth)
        {
            m_u[i] = m_hu[i] / m_h[i];
            m_v[i] = m_hv[i] / m_h[i];
        }
        else
        {
            // A dry or near-dry column carries no momentum: dividing a tiny
            // discharge by a tiny depth is how thin films generate absurd
            // velocities and destroy the time step.
            m_hu[i] = 0.0;
            m_hv[i] = 0.0;
            m_u[i]  = 0.0;
            m_v[i]  = 0.0;
        }

        m_eta[i] = bed[i] + m_h[i];
    }
}

double BarotropicSolver::stableTimeStep() const
{
    const Mesh::HorizontalMesh &mesh = *m_mesh;

    const std::vector<double> &area   = mesh.faceArea();
    const std::vector<double> &length = mesh.edgeLength();
    const std::vector<double> &nx     = mesh.edgeNormalX();
    const std::vector<double> &ny     = mesh.edgeNormalY();

    std::vector<double> speedSum(static_cast<size_t>(m_columnCount), 0.0);

    for (int64_t e = 0; e < mesh.edgeCount(); ++e)
    {
        const size_t ei    = static_cast<size_t>(e);
        const size_t left  = static_cast<size_t>(mesh.edgeLeftFace()[ei]);
        const int64_t right = mesh.edgeRightFace()[ei];

        const double normalX = nx[ei];
        const double normalY = ny[ei];

        if (m_h[left] > m_options.dryDepth)
        {
            const double un = m_u[left] * normalX + m_v[left] * normalY;
            speedSum[left] +=
                length[ei] * (std::abs(un) + std::sqrt(m_options.gravity * m_h[left]));
        }

        if (right != Mesh::kNoFace)
        {
            const size_t r = static_cast<size_t>(right);

            if (m_h[r] > m_options.dryDepth)
            {
                const double un = m_u[r] * normalX + m_v[r] * normalY;
                speedSum[r] +=
                    length[ei] * (std::abs(un) + std::sqrt(m_options.gravity * m_h[r]));
            }
        }
    }

    double smallest = std::numeric_limits<double>::infinity();

    for (int64_t c = 0; c < m_columnCount; ++c)
    {
        const size_t i = static_cast<size_t>(c);

        if (speedSum[i] <= 0.0)
            continue;

        smallest = std::min(smallest, area[i] / speedSum[i]);
    }

    if (std::isinf(smallest))
        return m_options.maxTimeStep;

    return std::clamp(m_options.cflNumber * smallest, m_options.minTimeStep,
                      m_options.maxTimeStep);
}

void BarotropicSolver::accumulateFluxes()
{
    const Mesh::HorizontalMesh &mesh = *m_mesh;

    // Only the boundary tally needs clearing. The derivative arrays are
    // assigned outright by the gather below, not accumulated into, so zeroing
    // them was dead work left over from the scattering version.
    std::fill(m_boundaryFlux.begin(), m_boundaryFlux.end(), 0.0);

    const std::vector<double> &bed    = mesh.faceBed();
    const std::vector<double> &area   = mesh.faceArea();
    const std::vector<double> &length = mesh.edgeLength();
    const std::vector<double> &nxs    = mesh.edgeNormalX();
    const std::vector<double> &nys    = mesh.edgeNormalY();

    const double gravity = m_options.gravity;

    const int64_t edgeCount = mesh.edgeCount();

    FVQual_PARALLEL_FOR
    for (int64_t e = 0; e < edgeCount; ++e)
    {
        const size_t  ei    = static_cast<size_t>(e);
        const size_t  left  = static_cast<size_t>(mesh.edgeLeftFace()[ei]);
        const int64_t right = mesh.edgeRightFace()[ei];

        const double nx = nxs[ei];
        const double ny = nys[ei];

        const double hL   = m_h[left];
        const double etaL = m_eta[left];
        const double zL   = bed[left];

        const double unL = m_u[left] * nx + m_v[left] * ny;
        const double utL = -m_u[left] * ny + m_v[left] * nx;

        // An unpaired edge is on the boundary. A wall is the mirror state —
        // same depth and bed, tangential velocity preserved, normal velocity
        // reversed — and an open boundary is whatever the characteristics
        // allow the exterior to impose (see ghostState).
        const bool   wall = right == Mesh::kNoFace;
        const size_t r    = wall ? left : static_cast<size_t>(right);

        double hR;
        double unR;
        double utR;

        const int boundary =
            wall && !m_edgeBoundary.empty() ? m_edgeBoundary[ei] : -1;

        if (!wall)
        {
            hR  = m_h[r];
            unR = m_u[r] * nx + m_v[r] * ny;
            utR = -m_u[r] * ny + m_v[r] * nx;
        }
        else if (boundary >= 0)
        {
            ghostState(m_boundaries[static_cast<size_t>(boundary)], hL, unL, utL, zL,
                       hR, unR, utR);
        }
        else
        {
            hR  = hL;
            unR = -unL;
            utR = utL;
        }

        // The exterior bed continues the interior one; an open boundary is a
        // cut through the domain, not a step in the bathymetry.
        const double zR   = wall ? zL : bed[r];
        const double etaR = wall ? zR + hR : m_eta[r];

        // Audusse hydrostatic reconstruction against the higher bed.
        const double zFace = std::max(zL, zR);
        const double hLs   = std::max(0.0, etaL - zFace);
        const double hRs   = std::max(0.0, etaR - zFace);

        const NormalState stateL{hLs, unL, utL};
        const NormalState stateR{hRs, unR, utR};

        const NormalFlux flux = hllc(stateL, stateR, gravity, m_options.dryDepth);

        // Rotate the flux back into Cartesian components.
        const double fMass = flux.mass;

        // Edge volume flux, signed left-to-right, retained so the internal
        // mode can transport scalars on a consistent depth-integrated
        // transport.
        m_edgeVolumeFlux[ei] = fMass * length[ei];

        if (boundary >= 0)
        {
            // Reported positive into the domain, which is the sign a mass
            // budget wants; the edge normal points out. This is the one
            // remaining scatter in the edge pass. It is over a handful of
            // boundary groups rather than every cell, and only boundary edges
            // reach it, so an atomic costs nothing measurable here where one
            // over the cell accumulators would have dominated.
            FVQual_ATOMIC
            m_boundaryFlux[static_cast<size_t>(boundary)] -= fMass * length[ei];
        }
        const double fMomX = flux.normal * nx - flux.tangential * ny;
        const double fMomY = flux.normal * ny + flux.tangential * nx;

        // Interface bed correction. Pairing g/2(h² - h*²) with the
        // reconstructed flux is what makes still water stay still: at rest
        // the two terms collapse to g/2 h² times the outward normal, whose
        // sum around a closed cell is identically zero.
        //
        // Stored per edge rather than scattered into the two adjacent cells.
        // The scatter is what made this loop unparallelisable: two edges
        // sharing a cell race on the same accumulator, and that is a property
        // of every unstructured finite-volume sweep, not of this one. Storing
        // here and gathering below costs three doubles per edge and makes both
        // passes race-free by construction — no colouring, no atomics, and the
        // shape a GPU wants (F7.1a).
        m_edgeFluxMass[ei] = fMass;
        m_edgeFluxX[ei]    = fMomX;
        m_edgeFluxY[ei]    = fMomY;

        m_edgeCorrectionL[ei] = 0.5 * gravity * (hL * hL - hLs * hLs);
        m_edgeCorrectionR[ei] = wall ? 0.0 : 0.5 * gravity * (hR * hR - hRs * hRs);
    }

    // Gather: each cell sums the edges it owns. Parallel over cells, and every
    // write below targets a cell only this iteration touches.
    const std::vector<int64_t> &edgeOffsets = mesh.faceEdgeOffsets();
    const std::vector<int64_t> &faceEdges   = mesh.faceEdges();

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < m_columnCount; ++c)
    {
        const size_t i = static_cast<size_t>(c);

        double dh  = 0.0;
        double dhu = 0.0;
        double dhv = 0.0;

        const int64_t begin = edgeOffsets[i];
        const int64_t end   = edgeOffsets[i + 1];

        for (int64_t slot = begin; slot < end; ++slot)
        {
            const size_t ei = static_cast<size_t>(faceEdges[static_cast<size_t>(slot)]);

            // The edge normal points out of its left face, so this cell takes
            // the flux with a minus sign when it is that left face.
            const bool   owns = mesh.edgeLeftFace()[ei] == c;
            const double sign = owns ? -1.0 : 1.0;

            const double correction = owns ? m_edgeCorrectionL[ei] : m_edgeCorrectionR[ei];
            const double scale      = sign * length[ei] / area[i];

            dh += scale * m_edgeFluxMass[ei];
            dhu += scale * (m_edgeFluxX[ei] + correction * nxs[ei]);
            dhv += scale * (m_edgeFluxY[ei] + correction * nys[ei]);
        }

        m_dh[i]  = dh;
        m_dhu[i] = dhu;
        m_dhv[i] = dhv;
    }
}

void BarotropicSolver::applySourceTerms(double dt)
{
    const double windAccelerationX = m_options.windStressX / m_options.waterDensity;
    const double windAccelerationY = m_options.windStressY / m_options.waterDensity;

    const bool wind = windAccelerationX != 0.0 || windAccelerationY != 0.0;
    const bool body = !m_bodyAccelerationX.empty();

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < m_columnCount; ++c)
    {
        const size_t i = static_cast<size_t>(c);

        if (m_h[i] <= m_options.dryDepth)
        {
            m_hu[i] = 0.0;
            m_hv[i] = 0.0;
            continue;
        }

        if (wind)
        {
            m_hu[i] += dt * windAccelerationX;
            m_hv[i] += dt * windAccelerationY;
        }

        if (body)
        {
            // A body force acts through the whole column, so it scales with
            // depth — unlike the surface wind stress above.
            m_hu[i] += dt * m_h[i] * m_bodyAccelerationX[i];
            m_hv[i] += dt * m_h[i] * m_bodyAccelerationY[i];
        }

        if (!m_options.frictionEnabled)
            continue;

        const double manning = m_manning[i];

        if (manning <= 0.0)
            continue;

        const double u     = m_hu[i] / m_h[i];
        const double v     = m_hv[i] / m_h[i];
        const double speed = std::hypot(u, v);

        if (speed <= 0.0)
            continue;

        // Semi-implicit Manning drag: linearizing about the current speed
        // gives a strictly positive denominator, so friction can slow the
        // flow to a halt but never reverse it — the overshoot an explicit
        // treatment produces in thin water.
        const double drag = m_options.gravity * manning * manning * speed /
                            std::pow(m_h[i], 4.0 / 3.0);

        const double factor = 1.0 / (1.0 + dt * drag);

        m_hu[i] *= factor;
        m_hv[i] *= factor;
    }
}

bool BarotropicSolver::step(double dt, std::string &message)
{
    if (m_mesh == nullptr)
    {
        message = "barotropic solver has not been initialized";
        return false;
    }

    if (!(dt > 0.0))
    {
        message = "time step must be positive";
        return false;
    }

    // Level-driven sources are re-rated against the surface this substep is
    // about to act on. m_eta was reconstructed at the end of the previous step
    // (or by initialize), so it is the current free surface here.
    if (m_columnSourceUpdater)
    {
        if (m_columnSource.empty())
            m_columnSource.assign(static_cast<size_t>(m_columnCount), 0.0);

        if (!m_columnSourceUpdater(m_eta, dt, m_columnSource, message))
            return false;
    }

    const std::vector<double> &area   = m_mesh->faceArea();
    const bool                 source = !m_columnSource.empty();

    accumulateFluxes();

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < m_columnCount; ++c)
    {
        const size_t i = static_cast<size_t>(c);

        m_h[i] += dt * m_dh[i];
        m_hu[i] += dt * m_dhu[i];
        m_hv[i] += dt * m_dhv[i];

        if (source && m_columnSource[i] != 0.0)
        {
            const double rate = m_columnSource[i] / area[i];

            m_h[i] += dt * rate;

            // Outflow carries the local velocity away with it, so the depth
            // falls but the speed does not rise. Inflow carries none, so it
            // dilutes — which is why only the negative branch touches
            // momentum.
            if (rate < 0.0)
            {
                m_hu[i] += dt * rate * m_u[i];
                m_hv[i] += dt * rate * m_v[i];
            }
        }

        if (m_h[i] < 0.0)
            m_h[i] = 0.0;
    }

    applySourceTerms(dt);
    reconstructPrimitives();

    for (size_t e = 0; e < m_edgeVolumeFlux.size(); ++e)
        m_edgeVolumeAccumulated[e] += m_edgeVolumeFlux[e] * dt;

    m_accumulationInterval += dt;
    m_time += dt;

    m_diagnostics.timeStep    = dt;
    m_diagnostics.totalVolume = totalVolume();
    m_diagnostics.wetColumns  = 0;
    m_diagnostics.maxCourant  = 0.0;

    const double stable = stableTimeStep();
    if (stable > 0.0)
        m_diagnostics.maxCourant = m_options.cflNumber * dt / stable;

    for (int64_t c = 0; c < m_columnCount; ++c)
        m_diagnostics.wetColumns += isWet(c) ? 1 : 0;

    message.clear();
    return true;
}

bool BarotropicSolver::advance(double duration, std::string &message)
{
    if (m_mesh == nullptr)
    {
        message = "barotropic solver has not been initialized";
        return false;
    }

    if (duration < 0.0)
    {
        message = "advance duration must be non-negative";
        return false;
    }

    const double target = m_time + duration;
    int64_t      steps  = 0;

    while (m_time < target)
    {
        // Trim the final substep so the solver lands exactly on the target:
        // the internal mode subcycles this solver and needs the two clocks
        // to agree to the bit.
        const double dt = std::min(stableTimeStep(), target - m_time);

        if (dt <= 0.0)
            break;

        if (!step(dt, message))
            return false;

        ++steps;
    }

    m_time                 = target;
    m_diagnostics.steps    = steps;

    message.clear();
    return true;
}

void BarotropicSolver::resetFluxAccumulation() noexcept
{
    std::fill(m_edgeVolumeAccumulated.begin(), m_edgeVolumeAccumulated.end(), 0.0);
    m_accumulationInterval = 0.0;
}

double BarotropicSolver::totalVolume() const noexcept
{
    const std::vector<double> &area = m_mesh->faceArea();

    double total = 0.0;

    for (int64_t c = 0; c < m_columnCount; ++c)
    {
        const size_t i = static_cast<size_t>(c);
        total += m_h[i] * area[i];
    }

    return total;
}

} // namespace FVQual::Hydro
