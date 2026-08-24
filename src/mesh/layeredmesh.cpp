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
    m_interfaceArea.assign(static_cast<size_t>(interfaceSlotCount()), 0.0);

    const std::vector<double> &area = m_horizontal->faceArea();

    for (int64_t c = 0; c < m_columnCount; ++c)
    {
        m_columnDepth[static_cast<size_t>(c)] =
            surfaceElevation[static_cast<size_t>(c)] - bed[static_cast<size_t>(c)];

        for (int k = 0; k <= m_layerCount; ++k)
        {
            m_interfaceArea[static_cast<size_t>(interfaceSlot(c, k))] =
                planArea(c, z(c, k));
        }

        for (int k = 0; k < m_layerCount; ++k)
        {
            const double dz = z(c, k) - z(c, k + 1);
            const size_t i  = static_cast<size_t>(cell(c, k));

            m_thickness[i] = dz;

            // The exact integral of the area curve over the cell, which reduces
            // to `dz * area` when the column has no curve — so a model without
            // hypsography gets bit-for-bit what it always did.
            m_cellVolume[i] = storageBetween(c, z(c, k + 1), z(c, k));
        }
    }

    message.clear();
    return true;
}

bool LayeredMesh::setHypsography(int64_t column, std::span<const double> elevation,
                                 std::span<const double> area, std::string &message)
{
    if (m_horizontal == nullptr)
    {
        // Sized from the horizontal mesh, which build() has not seen yet — so
        // the curve is stashed against a column index that may not exist. Refuse
        // rather than resize later and hope the indices still line up.
        message = "attach the horizontal mesh with build() before setting hypsography";
        return false;
    }

    if (column < 0 || column >= m_columnCount)
    {
        message = "column " + std::to_string(column) + " is out of range";
        return false;
    }

    if (elevation.size() != area.size())
    {
        message = "hypsography needs one area per elevation; got " +
                  std::to_string(elevation.size()) + " and " +
                  std::to_string(area.size());
        return false;
    }

    if (elevation.size() < 2)
    {
        message = "hypsography needs at least two samples to interpolate between";
        return false;
    }

    for (size_t i = 0; i < elevation.size(); ++i)
    {
        if (!(area[i] >= 0.0))
        {
            message = "hypsographic area must be non-negative";
            return false;
        }

        if (i == 0)
            continue;

        if (!(elevation[i] > elevation[i - 1]))
        {
            message = "hypsographic elevations must increase strictly";
            return false;
        }

        // A water body that widens downward is a cave, not a reservoir, and it
        // would make the storage curve non-monotone — so a free surface would
        // no longer follow uniquely from a volume.
        if (area[i] < area[i - 1])
        {
            message = "hypsographic area must not decrease with elevation";
            return false;
        }
    }

    m_hypsography.resize(static_cast<size_t>(m_columnCount));

    Hypsography &curve = m_hypsography[static_cast<size_t>(column)];

    curve.elevation.assign(elevation.begin(), elevation.end());
    curve.area.assign(area.begin(), area.end());

    m_hasHypsography = true;

    message.clear();
    return true;
}

namespace {

//! Shared empty table, so the accessors can return a reference either way.
const std::vector<double> kNoSamples;

} // namespace

const std::vector<double> &LayeredMesh::hypsographyElevation(int64_t column) const
{
    if (!m_hasHypsography || m_hypsography.empty())
        return kNoSamples;

    return m_hypsography[static_cast<size_t>(column)].elevation;
}

const std::vector<double> &LayeredMesh::hypsographyArea(int64_t column) const
{
    if (!m_hasHypsography || m_hypsography.empty())
        return kNoSamples;

    return m_hypsography[static_cast<size_t>(column)].area;
}

double LayeredMesh::planArea(int64_t column, double zLevel) const
{
    if (!m_hasHypsography || m_hypsography.empty())
        return m_horizontal->faceArea()[static_cast<size_t>(column)];

    const Hypsography &curve = m_hypsography[static_cast<size_t>(column)];

    if (curve.elevation.empty())
        return m_horizontal->faceArea()[static_cast<size_t>(column)];

    // Clamped rather than extrapolated: below the surveyed bed a linear
    // extrapolation runs negative, and above the highest sample it invents
    // storage the survey never saw.
    if (zLevel <= curve.elevation.front())
        return curve.area.front();

    if (zLevel >= curve.elevation.back())
        return curve.area.back();

    const auto upper =
        std::upper_bound(curve.elevation.begin(), curve.elevation.end(), zLevel);

    const size_t high = static_cast<size_t>(upper - curve.elevation.begin());
    const size_t low  = high - 1;

    const double span = curve.elevation[high] - curve.elevation[low];

    if (span <= 0.0)
        return curve.area[low];

    const double fraction = (zLevel - curve.elevation[low]) / span;

    return curve.area[low] + fraction * (curve.area[high] - curve.area[low]);
}

double LayeredMesh::storageBetween(int64_t column, double lower, double upper) const
{
    if (upper <= lower)
        return 0.0;

    if (!m_hasHypsography || m_hypsography.empty() ||
        m_hypsography[static_cast<size_t>(column)].elevation.empty())
    {
        return (upper - lower) *
               m_horizontal->faceArea()[static_cast<size_t>(column)];
    }

    const Hypsography &curve = m_hypsography[static_cast<size_t>(column)];

    // Integrate the piecewise-linear curve exactly by splitting at every sample
    // inside the range. A midpoint rule would be exact only where the range
    // contains no kink, and a cell straddling the thermocline break in a
    // surveyed curve is precisely where it would not.
    double total = 0.0;
    double at    = lower;

    while (at < upper)
    {
        double next = upper;

        // The next sample elevation strictly above `at`, if it lies inside.
        const auto above = std::upper_bound(curve.elevation.begin(),
                                            curve.elevation.end(), at);

        if (above != curve.elevation.end() && *above < next)
            next = *above;

        const double a0 = planArea(column, at);
        const double a1 = planArea(column, next);

        total += 0.5 * (a0 + a1) * (next - at);

        at = next;
    }

    return total;
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
