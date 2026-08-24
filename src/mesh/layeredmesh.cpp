// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/mesh/layeredmesh.h"

#include "hydrocouplesdk/tools/sigmagrid.h"

#include <algorithm>
#include <cmath>

namespace FVQual::Mesh {

namespace Tools = HydroCouple::SDK::Tools;

std::string HaneyReport::summary() const
{
    std::string text = "Haney steepness max r = " + std::to_string(maxRatio) +
                       " (limit " + std::to_string(limit) + ")";

    if (exceedances == 0)
        return text + " — pass";

    return text + " — " + std::to_string(exceedances) +
           " interior edge(s) exceed the limit, worst at edge " +
           std::to_string(worstEdge) +
           "; consider bathymetry smoothing or a hybrid z-sigma coordinate";
}

bool LayeredMesh::build(const HorizontalMesh &horizontal, const SigmaSpec &spec,
                        std::span<const double> surfaceElevation, std::string &message)
{
    if (horizontal.faceCount() <= 0)
    {
        message = "horizontal mesh has no faces";
        return false;
    }

    std::vector<double> sigma;
    if (!spec.resolve(sigma, message))
        return false;

    m_horizontal        = &horizontal;
    m_sigma             = std::move(sigma);
    m_layerCount        = static_cast<int>(m_sigma.size()) - 1;
    m_columnCount       = horizontal.faceCount();
    m_minLayerThickness = spec.minLayerThickness;

    return regenerate(surfaceElevation, message);
}

bool LayeredMesh::build(const HorizontalMesh &horizontal, const SigmaSpec &spec,
                        double uniformSurfaceElevation, std::string &message)
{
    const std::vector<double> surface(static_cast<size_t>(horizontal.faceCount()),
                                      uniformSurfaceElevation);

    return build(horizontal, spec, surface, message);
}

bool LayeredMesh::updateSurface(std::span<const double> surfaceElevation,
                                std::string &message)
{
    if (m_horizontal == nullptr)
    {
        message = "layered mesh has not been built";
        return false;
    }

    return regenerate(surfaceElevation, message);
}

bool LayeredMesh::regenerate(std::span<const double> surfaceElevation, std::string &message)
{
    if (static_cast<int64_t>(surfaceElevation.size()) != m_columnCount)
    {
        message = "surface elevation array has " +
                  std::to_string(surfaceElevation.size()) + " entries; expected " +
                  std::to_string(m_columnCount);
        return false;
    }

    const std::vector<double> &bed = m_horizontal->faceBed();

    for (int64_t c = 0; c < m_columnCount; ++c)
    {
        if (surfaceElevation[static_cast<size_t>(c)] <= bed[static_cast<size_t>(c)])
        {
            message = "column " + std::to_string(c) +
                      " has a free surface at or below its bed (" +
                      std::to_string(surfaceElevation[static_cast<size_t>(c)]) + " <= " +
                      std::to_string(bed[static_cast<size_t>(c)]) +
                      "); dry columns are handled from Phase 3 onward";
            return false;
        }
    }

    Tools::SigmaGrid grid;
    if (!Tools::SigmaGridGenerator::build(bed, surfaceElevation, m_sigma,
                                          m_minLayerThickness, grid, message))
    {
        return false;
    }

    m_interfaceZ = std::move(grid.interfaceZ);

    m_surfaceElevation.assign(surfaceElevation.begin(), surfaceElevation.end());
    m_columnDepth.assign(static_cast<size_t>(m_columnCount), 0.0);
    m_thickness.assign(static_cast<size_t>(cellCount()), 0.0);
    m_cellVolume.assign(static_cast<size_t>(cellCount()), 0.0);

    const std::vector<double> &area = m_horizontal->faceArea();

    for (int64_t c = 0; c < m_columnCount; ++c)
    {
        m_columnDepth[static_cast<size_t>(c)] =
            surfaceElevation[static_cast<size_t>(c)] - bed[static_cast<size_t>(c)];

        for (int k = 0; k < m_layerCount; ++k)
        {
            const double dz = z(c, k) - z(c, k + 1);
            const size_t i  = static_cast<size_t>(cell(c, k));

            m_thickness[i]  = dz;
            m_cellVolume[i] = dz * area[static_cast<size_t>(c)];
        }
    }

    message.clear();
    return true;
}

double LayeredMesh::totalVolume() const noexcept
{
    double total = 0.0;

    for (double volume : m_cellVolume)
        total += volume;

    return total;
}

HaneyReport LayeredMesh::haneyReport(double limit) const
{
    HaneyReport report;
    report.limit = limit;

    if (m_horizontal == nullptr)
        return report;

    const HorizontalMesh &mesh = *m_horizontal;

    for (int64_t e = 0; e < mesh.edgeCount(); ++e)
    {
        if (mesh.isBoundaryEdge(e))
            continue;

        const size_t left  = static_cast<size_t>(mesh.edgeLeftFace()[static_cast<size_t>(e)]);
        const size_t right = static_cast<size_t>(mesh.edgeRightFace()[static_cast<size_t>(e)]);

        const double hl = m_columnDepth[left];
        const double hr = m_columnDepth[right];
        const double sum = hl + hr;

        if (sum <= 0.0)
            continue;

        const double ratio = std::abs(hl - hr) / sum;

        if (ratio > report.maxRatio)
        {
            report.maxRatio  = ratio;
            report.worstEdge = e;
        }

        if (ratio > limit)
            ++report.exceedances;
    }

    return report;
}

} // namespace FVQual::Mesh
