// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/hydro/pressureclosure.h"

#include "fvqual/parallel.h"

#include <algorithm>
#include <cmath>

namespace FVQual::Hydro {

namespace {

/*!
 * \brief Minmod limiter: the smaller slope when they agree in sign, zero
 *        otherwise. Prevents the reconstruction from inventing an extremum
 *        across a sharp thermocline.
 */
double minmod(double a, double b) noexcept
{
    if (a * b <= 0.0)
        return 0.0;

    return std::abs(a) < std::abs(b) ? a : b;
}

} // namespace

std::string HydrostaticPressureClosure::name() const
{
    return m_options.linearDensityReconstruction ? "hydrostatic (linear reconstruction)"
                                                 : "hydrostatic (piecewise constant)";
}

void HydrostaticPressureClosure::prepareColumns(const Mesh::LayeredMesh &mesh,
                                                std::span<const double> density,
                                                double referenceDensity)
{
    const int64_t columns    = mesh.columnCount();
    const int     layers     = mesh.layerCount();
    const int     interfaces = mesh.interfaceCount();

    m_cumulative.assign(static_cast<size_t>(columns) * static_cast<size_t>(interfaces),
                        0.0);
    m_slope.assign(static_cast<size_t>(mesh.cellCount()), 0.0);

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < columns; ++c)
    {
        const size_t base = static_cast<size_t>(c) * static_cast<size_t>(interfaces);

        // Cumulative integral of the density *anomaly* downward from the free
        // surface. Integrating the full density instead would be a serious
        // error rather than a scaling choice: since
        //
        //     integral(z..eta) rho  =  rho0 (eta - z) + integral(z..eta) rho'
        //
        // the full-density form carries the entire barotropic surface-slope
        // term rho0 d(eta)/dx inside what is supposed to be the baroclinic
        // part. The external mode already computes that term through its own
        // pressure flux, so including it here hands the external mode a
        // duplicate of its own restoring force — and hands it *lagged*, held
        // fixed across a whole subcycling interval. A restoring force that
        // keeps pushing after the surface has already responded is a
        // resonantly forced oscillator, and the free surface grows without
        // bound.
        //
        // Note this is invisible while the surface is flat, which is why the
        // rest-state tests pass either way and only a tilted surface exposes
        // it.
        for (int k = 0; k < layers; ++k)
        {
            const double thickness = mesh.z(c, k) - mesh.z(c, k + 1);
            const double anomaly =
                density[static_cast<size_t>(mesh.cell(c, k))] - referenceDensity;

            m_cumulative[base + static_cast<size_t>(k) + 1] =
                m_cumulative[base + static_cast<size_t>(k)] + anomaly * thickness;
        }

        if (!m_options.linearDensityReconstruction || layers < 2)
            continue;

        // One-sided boundary slopes retain an exactly linear density profile.
        // Setting either to zero introduces different surface/bed integrals
        // in columns with different sigma spacing, even in resting water.
        for (int k : {0,layers-1}) {
            const int neighbour=k==0 ? 1 : layers-2;
            const double centre=.5*(mesh.z(c,k)+mesh.z(c,k+1));
            const double adjacent=.5*(mesh.z(c,neighbour)+mesh.z(c,neighbour+1));
            const size_t cell=static_cast<size_t>(mesh.cell(c,k));
            const size_t other=static_cast<size_t>(mesh.cell(c,neighbour));
            m_slope[cell]=(density[cell]-density[other])/(centre-adjacent);
        }

        for (int k = 1; k + 1 < layers; ++k)
        {
            const size_t above  = static_cast<size_t>(mesh.cell(c, k - 1));
            const size_t centre = static_cast<size_t>(mesh.cell(c, k));
            const size_t below  = static_cast<size_t>(mesh.cell(c, k + 1));

            const double zAbove  = 0.5 * (mesh.z(c, k - 1) + mesh.z(c, k));
            const double zCentre = 0.5 * (mesh.z(c, k) + mesh.z(c, k + 1));
            const double zBelow  = 0.5 * (mesh.z(c, k + 1) + mesh.z(c, k + 2));

            const double upper = zAbove - zCentre;
            const double lower = zCentre - zBelow;

            if (upper <= 0.0 || lower <= 0.0)
                continue;

            m_slope[centre] = minmod((density[above] - density[centre]) / upper,
                                     (density[centre] - density[below]) / lower);
        }
    }
}

double HydrostaticPressureClosure::columnIntegral(const Mesh::LayeredMesh &mesh,
                                                  std::span<const double> density,
                                                  double referenceDensity, int64_t column,
                                                  double z) const
{
    const int    layers     = mesh.layerCount();
    const int    interfaces = mesh.interfaceCount();
    const size_t base = static_cast<size_t>(column) * static_cast<size_t>(interfaces);

    const double surface = mesh.z(column, 0);
    const double bed     = mesh.z(column, layers);

    if (z >= surface)
        return 0.0;

    if (z <= bed)
    {
        // Below the column's bed the deepest density is extended downward, so
        // a face between columns of unequal depth still compares like with
        // like rather than comparing water against nothing.
        const double deepest =
            density[static_cast<size_t>(mesh.cell(column, layers - 1))] - referenceDensity;
        return m_cumulative[base + static_cast<size_t>(layers)] + deepest * (bed - z);
    }

    // Locate the layer containing z: interface k is the top of layer k.
    int low  = 0;
    int high = layers;

    while (low + 1 < high)
    {
        const int middle = (low + high) / 2;

        if (mesh.z(column, middle) > z)
            low = middle;
        else
            high = middle;
    }

    const size_t cell    = static_cast<size_t>(mesh.cell(column, low));
    const double topZ    = mesh.z(column, low);
    const double centreZ = 0.5 * (mesh.z(column, low) + mesh.z(column, low + 1));
    const double rho     = density[cell] - referenceDensity;

    // Integral of rho_k + s_k (z' - z_centre) from z up to the layer top.
    const double span   = topZ - z;
    const double linear = 0.5 * m_slope[cell] *
                          ((topZ - centreZ) * (topZ - centreZ) -
                           (z - centreZ) * (z - centreZ));

    return m_cumulative[base + static_cast<size_t>(low)] + rho * span + linear;
}

void HydrostaticPressureClosure::baroclinicAcceleration(const Mesh::LayeredMesh &mesh,
                                                        std::span<const double> density,
                                                        double referenceDensity,
                                                        double gravity,
                                                        std::vector<double> &accelerationX,
                                                        std::vector<double> &accelerationY)
{
    const Mesh::HorizontalMesh &horizontal = mesh.horizontal();

    const int64_t columns = mesh.columnCount();
    const int     layers  = mesh.layerCount();

    accelerationX.assign(static_cast<size_t>(mesh.cellCount()), 0.0);
    accelerationY.assign(static_cast<size_t>(mesh.cellCount()), 0.0);

    prepareColumns(mesh, density, referenceDensity);

    const std::vector<double>  &area            = horizontal.faceArea();
    const std::vector<double>  &edgeLength      = horizontal.edgeLength();
    const std::vector<double>  &normalX         = horizontal.edgeNormalX();
    const std::vector<double>  &normalY         = horizontal.edgeNormalY();
    const std::vector<int64_t> &edgeLeft        = horizontal.edgeLeftFace();
    const std::vector<int64_t> &edgeRight       = horizontal.edgeRightFace();
    const std::vector<int64_t> &faceEdgeOffsets = horizontal.faceEdgeOffsets();
    const std::vector<int64_t> &faceEdges       = horizontal.faceEdges();

    const double scale = -gravity / referenceDensity;

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < columns; ++c)
    {
        const int64_t begin = faceEdgeOffsets[static_cast<size_t>(c)];
        const int64_t end   = faceEdgeOffsets[static_cast<size_t>(c) + 1];

        for (int k = 0; k < layers; ++k)
        {
            // Evaluate every face of this cell at the cell's own layer-centre
            // elevation. Doing it at each face's own elevation instead is
            // precisely the mistake that makes sigma coordinates leak
            // spurious currents.
            const double z = 0.5 * (mesh.z(c, k) + mesh.z(c, k + 1));

            const double here = columnIntegral(mesh, density, referenceDensity, c, z);

            double gradientX = 0.0;
            double gradientY = 0.0;

            for (int64_t index = begin; index < end; ++index)
            {
                const int64_t e  = faceEdges[static_cast<size_t>(index)];
                const size_t  ei = static_cast<size_t>(e);

                const bool    owned = edgeLeft[ei] == c;
                const int64_t other = owned ? edgeRight[ei] : edgeLeft[ei];
                const double  sign  = owned ? 1.0 : -1.0;

                const double facePressure =
                    // Below the neighbour's bed this face is a solid sidewall.
                    // Extrapolating the neighbour's water below that bed adds a
                    // fictitious pressure contrast over a steep depth step.
                    (other == Mesh::kNoFace || z <= mesh.z(other,layers))
                        ? here
                        : 0.5 * (here + columnIntegral(mesh, density, referenceDensity, other, z));

                const double weight = edgeLength[ei] * sign;

                gradientX += weight * facePressure * normalX[ei];
                gradientY += weight * facePressure * normalY[ei];
            }

            const double inverseArea = 1.0 / area[static_cast<size_t>(c)];
            const size_t cell        = static_cast<size_t>(mesh.cell(c, k));

            accelerationX[cell] = scale * gradientX * inverseArea;
            accelerationY[cell] = scale * gradientY * inverseArea;
        }
    }
}

} // namespace FVQual::Hydro
