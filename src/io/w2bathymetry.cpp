// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/io/w2bathymetry.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

namespace FVQual::IO {

namespace {

//! Splits a CSV line into trimmed fields.
std::vector<std::string> splitCsv(const std::string &line)
{
    std::vector<std::string> out;
    std::string              field;
    std::istringstream       in(line);

    while (std::getline(in, field, ','))
    {
        const auto first = field.find_first_not_of(" \t\r\n");

        if (first == std::string::npos)
        {
            out.emplace_back();
            continue;
        }

        const auto last = field.find_last_not_of(" \t\r\n");

        out.push_back(field.substr(first, last - first + 1));
    }

    return out;
}

//! Reads the numeric fields after the label, stopping at the first blank.
bool readRow(const std::vector<std::string> &fields, size_t expected,
             std::vector<double> &out, std::string &why)
{
    out.clear();

    for (size_t i = 1; i < fields.size() && out.size() < expected; ++i)
    {
        if (fields[i].empty())
            continue;

        try
        {
            out.push_back(std::stod(fields[i]));
        }
        catch (const std::exception &)
        {
            why = "'" + fields[i] + "' is not a number";
            return false;
        }
    }

    if (out.size() != expected)
    {
        why = "expected " + std::to_string(expected) + " values, found " +
              std::to_string(out.size());
        return false;
    }

    return true;
}

std::string upper(std::string text)
{
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return text;
}

} // namespace

bool W2Bathymetry::read(const std::string &path, double bottomElevation,
                        std::string &message)
{
    std::ifstream input(path);

    if (!input)
    {
        message = "cannot open bathymetry file '" + path + "'";
        return false;
    }

    return parse(input, path, bottomElevation, message);
}

bool W2Bathymetry::parse(std::istream &input, const std::string &origin,
                         double bottomElevation, std::string &message)
{
    std::vector<double>              length;
    std::vector<double>              orientation;
    std::vector<double>              friction;
    std::vector<double>              surface;
    std::vector<double>              layerHeight;
    std::vector<std::vector<double>> width;

    size_t segments = 0;
    int    number   = 0;
    bool   inLayers = false;

    std::string line;
    std::string why;

    const auto fail = [&](const std::string &reason) {
        message = origin + ":" + std::to_string(number) + ": " + reason;
        return false;
    };

    while (std::getline(input, line))
    {
        ++number;

        const std::vector<std::string> fields = splitCsv(line);

        if (fields.empty())
            continue;

        const std::string label = upper(fields.front());

        if (label == "$")
            continue;

        if (!inLayers)
        {
            if (label.empty())
            {
                // The segment-index row. Its populated entries fix how many
                // segments every later row must carry.
                segments = 0;

                for (size_t i = 1; i < fields.size(); ++i)
                {
                    if (!fields[i].empty())
                        ++segments;
                }

                if (segments == 0)
                    return fail("the segment-index row lists no segments");

                continue;
            }

            if (segments == 0)
                return fail("'" + label + "' appears before the segment-index row");

            if (label == "DLX")
            {
                if (!readRow(fields, segments, length, why))
                    return fail("DLX: " + why);
            }
            else if (label == "ELWS")
            {
                if (!readRow(fields, segments, surface, why))
                    return fail("ELWS: " + why);
            }
            else if (label == "PHI0")
            {
                if (!readRow(fields, segments, orientation, why))
                    return fail("PHI0: " + why);
            }
            else if (label == "FRIC")
            {
                if (!readRow(fields, segments, friction, why))
                    return fail("FRIC: " + why);
            }
            else if (label == "LAYERH")
            {
                inLayers = true;
            }
            else
            {
                return fail("unknown row '" + fields.front() + "'");
            }

            continue;
        }

        // Layer rows: thickness, then one width per segment, then the layer
        // number. The trailing K is ignored — it is a label, and trusting the
        // row order is what every other reader of this format does.
        if (fields.front().empty())
            continue;

        double thickness = 0.0;

        try
        {
            thickness = std::stod(fields.front());
        }
        catch (const std::exception &)
        {
            return fail("layer thickness '" + fields.front() + "' is not a number");
        }

        if (!(thickness > 0.0))
            return fail("layer thickness must be positive");

        std::vector<double> row;

        if (!readRow(fields, segments, row, why))
            return fail("layer widths: " + why);

        for (double w : row)
        {
            if (w < 0.0)
                return fail("a layer width is negative");
        }

        layerHeight.push_back(thickness);
        width.push_back(std::move(row));
    }

    if (length.empty() || layerHeight.empty())
    {
        message = origin + ": no segments or no layers were read";
        return false;
    }

    if (orientation.empty())
        orientation.assign(segments, 0.0);

    if (friction.empty())
        friction.assign(segments, 0.03);

    if (surface.empty())
        surface.assign(segments, bottomElevation);

    m_length         = std::move(length);
    m_orientation    = std::move(orientation);
    m_friction       = std::move(friction);
    m_initialSurface = std::move(surface);
    m_layerHeight    = std::move(layerHeight);
    m_width          = std::move(width);
    m_bottom         = bottomElevation;

    // W2 pads a branch with zero-width segments at each end; they are grid
    // bookkeeping, not reservoir, and become no column.
    m_active.clear();

    for (int i = 0; i < segmentCount(); ++i)
    {
        double widest = 0.0;

        for (int k = 0; k < layerCount(); ++k)
            widest = std::max(widest, m_width[static_cast<size_t>(k)][static_cast<size_t>(i)]);

        if (widest > 0.0)
            m_active.push_back(i);
    }

    if (m_active.empty())
    {
        message = origin + ": every segment has zero width";
        return false;
    }

    message.clear();
    return true;
}

double W2Bathymetry::width(int layer, int segment) const
{
    if (layer < 0 || layer >= layerCount() || segment < 0 || segment >= segmentCount())
        return 0.0;

    return m_width[static_cast<size_t>(layer)][static_cast<size_t>(segment)];
}

double W2Bathymetry::layerBottom(int layer) const
{
    // Layer 0 is the top of the grid, so accumulate upward from the bed.
    double elevation = m_bottom;

    for (int k = layerCount() - 1; k > layer; --k)
        elevation += m_layerHeight[static_cast<size_t>(k)];

    return elevation;
}

double W2Bathymetry::volumeAt(double surfaceElevation) const
{
    double total = 0.0;

    for (int k = 0; k < layerCount(); ++k)
    {
        const double bottom = layerBottom(k);
        const double top    = bottom + m_layerHeight[static_cast<size_t>(k)];

        const double wetted = std::min(top, surfaceElevation) - bottom;

        if (wetted <= 0.0)
            continue;

        for (int i : m_active)
        {
            total += width(k, i) * m_length[static_cast<size_t>(i)] * wetted;
        }
    }

    return total;
}

double W2Bathymetry::surfaceAreaAt(double surfaceElevation) const
{
    double total = 0.0;

    for (int i : m_active)
    {
        for (int k = 0; k < layerCount(); ++k)
        {
            const double bottom = layerBottom(k);
            const double top    = bottom + m_layerHeight[static_cast<size_t>(k)];

            if (surfaceElevation > bottom && surfaceElevation <= top && width(k, i) > 0.0)
            {
                total += width(k, i) * m_length[static_cast<size_t>(i)];
                break;
            }
        }
    }

    return total;
}

double W2Bathymetry::surfaceWidth(int segment) const
{
    const double wse = m_initialSurface[static_cast<size_t>(segment)];

    for (int k = 0; k < layerCount(); ++k)
    {
        const double bottom = layerBottom(k);
        const double top    = bottom + m_layerHeight[static_cast<size_t>(k)];

        if (wse > bottom && wse <= top && width(k, segment) > 0.0)
            return width(k, segment);
    }

    // The surface may sit above every wetted layer; fall back to the widest.
    double widest = 0.0;

    for (int k = 0; k < layerCount(); ++k)
        widest = std::max(widest, width(k, segment));

    return widest;
}

bool W2Bathymetry::toMeshDefinition(HydroCouple::SDK::IO::MeshDefinition &definition,
                                    std::vector<double> &bed, std::string &message) const
{
    if (m_active.empty())
    {
        message = "no active segments to build a mesh from";
        return false;
    }

    definition = HydroCouple::SDK::IO::MeshDefinition{};
    definition.meshName = "mesh2d";

    bed.clear();

    // Walk the chain of segment centers, then hang a quad off each interval.
    // PHI0 is measured clockwise from north, so the along-channel unit vector is
    // (sin, cos) and the cross-channel one is its perpendicular.
    const size_t count = m_active.size();

    std::vector<double> centreX(count + 1, 0.0);
    std::vector<double> centreY(count + 1, 0.0);
    std::vector<double> normalX(count + 1, 0.0);
    std::vector<double> normalY(count + 1, 0.0);

    double x = 0.0;
    double y = 0.0;

    for (size_t j = 0; j <= count; ++j)
    {
        // At an interior boundary the direction is that of the segment just
        // upstream of it; at the ends, the single adjacent segment's.
        const size_t segment = std::min(j, count - 1);
        const int    index   = m_active[segment];

        const double phi = m_orientation[static_cast<size_t>(index)];

        const double dx = std::sin(phi);
        const double dy = std::cos(phi);

        centreX[j] = x;
        centreY[j] = y;

        // Perpendicular, so the quad spans the channel width.
        normalX[j] = -dy;
        normalY[j] = dx;

        if (j < count)
        {
            const double length = m_length[static_cast<size_t>(m_active[segment])];

            x += dx * length;
            y += dy * length;
        }
    }

    // Half-width at each boundary, from the widths of the adjacent segments at
    // the initial water surface: that is the plan geometry the model runs at.
    for (size_t j = 0; j <= count; ++j)
    {
        double half = 0.0;

        if (j == 0)
            half = 0.5 * surfaceWidth(m_active[0]);
        else if (j == count)
            half = 0.5 * surfaceWidth(m_active[count - 1]);
        else
            half = 0.25 * (surfaceWidth(m_active[j - 1]) + surfaceWidth(m_active[j]));

        if (!(half > 0.0))
        {
            message = "segment " + std::to_string(m_active[std::min(j, count - 1)]) +
                      " has no width at the initial water surface";
            return false;
        }

        definition.nodeX.push_back(centreX[j] - half * normalX[j]);
        definition.nodeY.push_back(centreY[j] - half * normalY[j]);
        definition.nodeX.push_back(centreX[j] + half * normalX[j]);
        definition.nodeY.push_back(centreY[j] + half * normalY[j]);
    }

    definition.faceNodeOffsets.push_back(0);

    for (size_t j = 0; j < count; ++j)
    {
        const int64_t left0  = static_cast<int64_t>(2 * j);
        const int64_t right0 = static_cast<int64_t>(2 * j + 1);
        const int64_t left1  = static_cast<int64_t>(2 * (j + 1));
        const int64_t right1 = static_cast<int64_t>(2 * (j + 1) + 1);

        definition.faceNodes.push_back(left0);
        definition.faceNodes.push_back(left1);
        definition.faceNodes.push_back(right1);
        definition.faceNodes.push_back(right0);

        definition.faceNodeOffsets.push_back(
            static_cast<int64_t>(definition.faceNodes.size()));

        // The bed is the bottom of the deepest layer that carries water.
        const int segment = m_active[j];

        double lowest = m_initialSurface[static_cast<size_t>(segment)];

        for (int k = 0; k < layerCount(); ++k)
        {
            if (width(k, segment) > 0.0)
                lowest = std::min(lowest, layerBottom(k));
        }

        bed.push_back(lowest);
    }

    message.clear();
    return true;
}

bool W2Bathymetry::applyHypsography(Mesh::LayeredMesh &mesh, std::string &message) const
{
    if (static_cast<size_t>(mesh.columnCount()) != m_active.size())
    {
        message = "the mesh has " + std::to_string(mesh.columnCount()) +
                  " columns but the bathymetry has " + std::to_string(m_active.size()) +
                  " active segments";
        return false;
    }

    for (size_t j = 0; j < m_active.size(); ++j)
    {
        const int segment = m_active[j];

        std::vector<double> elevation;
        std::vector<double> area;

        // Sampled at layer midpoints, deepest first. The integral of a linear
        // interpolant over a layer equals its midpoint value times the
        // thickness, so each layer's volume comes out at exactly B x DLX x H —
        // sampling at layer boundaries would average adjacent widths and lose
        // storage wherever the reservoir narrows.
        for (int k = layerCount() - 1; k >= 0; --k)
        {
            const double w = width(k, segment);

            if (w <= 0.0)
                continue;

            const double bottom = layerBottom(k);
            const double middle = bottom + 0.5 * m_layerHeight[static_cast<size_t>(k)];

            const double planArea = w * m_length[static_cast<size_t>(segment)];

            if (!elevation.empty() && planArea < area.back())
            {
                // W2 permits a segment to narrow upward; a storage curve cannot,
                // because it would make stage a non-unique function of volume.
                // Hold the wider value rather than refusing the whole file.
                elevation.push_back(middle);
                area.push_back(area.back());
                continue;
            }

            elevation.push_back(middle);
            area.push_back(planArea);
        }

        if (elevation.size() < 2)
        {
            message = "segment " + std::to_string(segment) +
                      " has fewer than two wetted layers, so it has no storage curve";
            return false;
        }

        if (!mesh.setHypsography(static_cast<int64_t>(j), elevation, area, message))
            return false;
    }

    message.clear();
    return true;
}


namespace {

/*!
 * \brief Overlap of lane [a, b] with the section [-halfWidth, +halfWidth] (m).
 *
 * The outermost lanes extend outward without bound, which is what makes the
 * shares sum to the full width `2 * halfWidth` at *every* stage rather than only
 * at the one the lane spans were laid out for. It matters: in DeGray all thirty
 * segments are wider below the top of the grid than at the initial surface — by
 * up to 4.1x — because the grid extends 14 m above the water. Lanes bounded at
 * the surface width would silently discard that section as the stage rose, and a
 * single-lane mesh would then not reproduce the laterally averaged model it is
 * supposed to be the control for.
 *
 * Physically it says the outer cells take the extra width when the reservoir
 * spreads above its initial stage, which is what a bank does.
 */
double lateralOverlap(double a, double b, double halfWidth, bool openInner,
                      bool openOuter)
{
    const double lower = openInner ? -halfWidth : std::max(a, -halfWidth);
    const double upper = openOuter ? halfWidth : std::min(b, halfWidth);

    return std::max(0.0, upper - lower);
}

} // namespace

bool W2Bathymetry::toLateralMeshDefinition(int lateralCells,
                                           HydroCouple::SDK::IO::MeshDefinition &definition,
                                           std::vector<double> &bed,
                                           std::string &message) const
{
    if (lateralCells < 1)
    {
        message = "lateralCells must be at least 1";
        return false;
    }

    if (m_active.empty())
    {
        message = "no active segments to build a mesh from";
        return false;
    }

    definition          = HydroCouple::SDK::IO::MeshDefinition{};
    definition.meshName = "mesh2d";

    bed.clear();

    const size_t count = m_active.size();
    const size_t lanes = static_cast<size_t>(lateralCells);

    // Along-channel walk, identical to toMeshDefinition(): PHI0 is clockwise
    // from north, so the along-channel unit vector is (sin, cos).
    std::vector<double> centreX(count + 1, 0.0);
    std::vector<double> centreY(count + 1, 0.0);
    std::vector<double> normalX(count + 1, 0.0);
    std::vector<double> normalY(count + 1, 0.0);

    double x = 0.0;
    double y = 0.0;

    for (size_t j = 0; j <= count; ++j)
    {
        const size_t segment = std::min(j, count - 1);
        const double phi     = m_orientation[static_cast<size_t>(m_active[segment])];

        const double dx = std::sin(phi);
        const double dy = std::cos(phi);

        centreX[j] = x;
        centreY[j] = y;
        normalX[j] = -dy;
        normalY[j] = dx;

        if (j < count)
        {
            x += dx * m_length[static_cast<size_t>(m_active[segment])];
            y += dy * m_length[static_cast<size_t>(m_active[segment])];
        }
    }

    for (size_t j = 0; j <= count; ++j)
    {
        double half = 0.0;

        if (j == 0)
            half = 0.5 * surfaceWidth(m_active[0]);
        else if (j == count)
            half = 0.5 * surfaceWidth(m_active[count - 1]);
        else
            half = 0.25 * (surfaceWidth(m_active[j - 1]) + surfaceWidth(m_active[j]));

        if (!(half > 0.0))
        {
            message = "segment " + std::to_string(m_active[std::min(j, count - 1)]) +
                      " has no width at the initial water surface";
            return false;
        }

        for (size_t l = 0; l <= lanes; ++l)
        {
            const double t      = static_cast<double>(l) / static_cast<double>(lanes);
            const double offset = -half + 2.0 * half * t;

            definition.nodeX.push_back(centreX[j] + offset * normalX[j]);
            definition.nodeY.push_back(centreY[j] + offset * normalY[j]);
        }
    }

    definition.faceNodeOffsets.push_back(0);

    const auto node = [&](size_t j, size_t l) {
        return static_cast<int64_t>(j * (lanes + 1) + l);
    };

    for (size_t j = 0; j < count; ++j)
    {
        const int    segment = m_active[j];
        const double halfTop = 0.5 * surfaceWidth(segment);

        for (size_t l = 0; l < lanes; ++l)
        {
            definition.faceNodes.push_back(node(j, l));
            definition.faceNodes.push_back(node(j + 1, l));
            definition.faceNodes.push_back(node(j + 1, l + 1));
            definition.faceNodes.push_back(node(j, l + 1));

            definition.faceNodeOffsets.push_back(
                static_cast<int64_t>(definition.faceNodes.size()));

            // This lane's span, in the segment's own coordinates.
            const double a = -halfTop + 2.0 * halfTop * static_cast<double>(l) /
                                            static_cast<double>(lanes);
            const double b = -halfTop + 2.0 * halfTop * static_cast<double>(l + 1) /
                                            static_cast<double>(lanes);

            // The bed is the lowest layer whose section still reaches this lane.
            // An outer lane therefore sits higher than a central one, which is
            // the cross-section shape a laterally averaged column cannot hold.
            double lowest = m_initialSurface[static_cast<size_t>(segment)];

            for (int k = 0; k < layerCount(); ++k)
            {
                if (lateralOverlap(a, b, 0.5 * width(k, segment), l == 0,
                                   l + 1 == lanes) > 0.0)
                    lowest = std::min(lowest, layerBottom(k));
            }

            bed.push_back(lowest);
        }
    }

    message.clear();
    return true;
}

bool W2Bathymetry::applyLateralHypsography(Mesh::LayeredMesh &mesh, int lateralCells,
                                           std::string &message) const
{
    if (lateralCells < 1)
    {
        message = "lateralCells must be at least 1";
        return false;
    }

    const size_t lanes    = static_cast<size_t>(lateralCells);
    const size_t expected = m_active.size() * lanes;

    if (static_cast<size_t>(mesh.columnCount()) != expected)
    {
        message = "the mesh has " + std::to_string(mesh.columnCount()) +
                  " columns but " + std::to_string(m_active.size()) +
                  " active segments at " + std::to_string(lateralCells) +
                  " cells wide needs " + std::to_string(expected);
        return false;
    }

    for (size_t j = 0; j < m_active.size(); ++j)
    {
        const int    segment = m_active[j];
        const double halfTop = 0.5 * surfaceWidth(segment);
        const double length  = m_length[static_cast<size_t>(segment)];

        for (size_t l = 0; l < lanes; ++l)
        {
            const double a = -halfTop + 2.0 * halfTop * static_cast<double>(l) /
                                            static_cast<double>(lanes);
            const double b = -halfTop + 2.0 * halfTop * static_cast<double>(l + 1) /
                                            static_cast<double>(lanes);

            std::vector<double> elevation;
            std::vector<double> area;

            // Layer midpoints, deepest first, exactly as applyHypsography does
            // and for the same reason: the integral of a linear interpolant over
            // a layer is its midpoint value times the thickness, so each layer's
            // volume comes out exact.
            for (int k = layerCount() - 1; k >= 0; --k)
            {
                const double overlap = lateralOverlap(a, b, 0.5 * width(k, segment),
                                                      l == 0, l + 1 == lanes);

                if (overlap <= 0.0)
                    continue;

                const double bottom = layerBottom(k);
                const double middle =
                    bottom + 0.5 * m_layerHeight[static_cast<size_t>(k)];

                const double planArea = overlap * length;

                if (!elevation.empty() && planArea < area.back())
                {
                    elevation.push_back(middle);
                    area.push_back(area.back());
                    continue;
                }

                elevation.push_back(middle);
                area.push_back(planArea);
            }

            if (elevation.size() < 2)
            {
                message = "segment " + std::to_string(segment) + " lane " +
                          std::to_string(l) +
                          " has fewer than two wetted layers, so it has no storage "
                          "curve; try fewer lateral cells";
                return false;
            }

            const int64_t column = static_cast<int64_t>(j * lanes + l);

            if (!mesh.setHypsography(column, elevation, area, message))
                return false;
        }
    }

    message.clear();
    return true;
}


W2Bathymetry::LateralPlan W2Bathymetry::planLateralCells(int maxLanes,
                                                         double minCellWidth) const
{
    LateralPlan plan;

    plan.lanes.reserve(m_active.size());
    plan.offset.reserve(m_active.size());

    const int ceiling = std::max(1, maxLanes);

    for (size_t j = 0; j < m_active.size(); ++j)
    {
        int lanes = ceiling;

        if (minCellWidth > 0.0)
        {
            const double w = surfaceWidth(m_active[j]);
            const int    fits =
                static_cast<int>(std::floor(w / minCellWidth));

            lanes = std::clamp(fits, 1, ceiling);
        }

        plan.offset.push_back(plan.columnCount);
        plan.lanes.push_back(lanes);
        plan.columnCount += lanes;
    }

    return plan;
}

namespace {

/*!
 * \brief Sorted union of two even partitions of [0, 1], deduplicated.
 *
 * These are the lane boundaries a shared segment boundary must carry so that a
 * coarse cell and its finer neighbor still meet edge to edge.
 */
std::vector<double> partitionUnion(int left, int right)
{
    std::vector<double> f;

    for (int i = 0; i <= left; ++i)
        f.push_back(static_cast<double>(i) / static_cast<double>(left));

    for (int i = 0; i <= right; ++i)
        f.push_back(static_cast<double>(i) / static_cast<double>(right));

    std::sort(f.begin(), f.end());

    f.erase(std::unique(f.begin(), f.end(),
                        [](double a, double b) { return std::abs(a - b) < 1.0e-12; }),
            f.end());

    return f;
}

} // namespace

bool W2Bathymetry::toPlannedMeshDefinition(const LateralPlan &plan,
                                           HydroCouple::SDK::IO::MeshDefinition &definition,
                                           std::vector<double> &bed,
                                           std::string &message) const
{
    if (m_active.empty())
    {
        message = "no active segments to build a mesh from";
        return false;
    }

    if (plan.lanes.size() != m_active.size())
    {
        message = "the lateral plan covers " + std::to_string(plan.lanes.size()) +
                  " segments but the bathymetry has " + std::to_string(m_active.size());
        return false;
    }

    definition          = HydroCouple::SDK::IO::MeshDefinition{};
    definition.meshName = "mesh2d";

    bed.clear();

    const size_t count = m_active.size();

    std::vector<double> centreX(count + 1, 0.0);
    std::vector<double> centreY(count + 1, 0.0);
    std::vector<double> normalX(count + 1, 0.0);
    std::vector<double> normalY(count + 1, 0.0);

    double x = 0.0;
    double y = 0.0;

    for (size_t j = 0; j <= count; ++j)
    {
        const size_t segment = std::min(j, count - 1);
        const double phi     = m_orientation[static_cast<size_t>(m_active[segment])];

        const double dx = std::sin(phi);
        const double dy = std::cos(phi);

        centreX[j] = x;
        centreY[j] = y;
        normalX[j] = -dy;
        normalY[j] = dx;

        if (j < count)
        {
            x += dx * m_length[static_cast<size_t>(m_active[segment])];
            y += dy * m_length[static_cast<size_t>(m_active[segment])];
        }
    }

    // Each boundary carries the union of the lane partitions on either side, so
    // a four-lane segment meeting a one-lane one still shares whole edges.
    std::vector<std::vector<double>> fractions(count + 1);
    std::vector<int64_t>             firstNode(count + 1, 0);

    for (size_t j = 0; j <= count; ++j)
    {
        const int left  = j == 0 ? plan.lanes.front() : plan.lanes[j - 1];
        const int right = j == count ? plan.lanes.back() : plan.lanes[j];

        fractions[j] = partitionUnion(left, right);

        double half = 0.0;

        if (j == 0)
            half = 0.5 * surfaceWidth(m_active[0]);
        else if (j == count)
            half = 0.5 * surfaceWidth(m_active[count - 1]);
        else
            half = 0.25 * (surfaceWidth(m_active[j - 1]) + surfaceWidth(m_active[j]));

        if (!(half > 0.0))
        {
            message = "segment " + std::to_string(m_active[std::min(j, count - 1)]) +
                      " has no width at the initial water surface";
            return false;
        }

        firstNode[j] = static_cast<int64_t>(definition.nodeX.size());

        for (double f : fractions[j])
        {
            const double offset = -half + 2.0 * half * f;

            definition.nodeX.push_back(centreX[j] + offset * normalX[j]);
            definition.nodeY.push_back(centreY[j] + offset * normalY[j]);
        }
    }

    definition.faceNodeOffsets.push_back(0);

    // Indices in `fractions[j]` lying within [lo, hi], inclusive of both ends.
    const auto span = [&](size_t boundary, double lo, double hi) {
        std::vector<int64_t> nodes;

        for (size_t i = 0; i < fractions[boundary].size(); ++i)
        {
            const double f = fractions[boundary][i];

            if (f >= lo - 1.0e-12 && f <= hi + 1.0e-12)
                nodes.push_back(firstNode[boundary] + static_cast<int64_t>(i));
        }

        return nodes;
    };

    for (size_t j = 0; j < count; ++j)
    {
        const int    segment = m_active[j];
        const int    lanes   = plan.lanes[j];
        const double halfTop = 0.5 * surfaceWidth(segment);

        for (int l = 0; l < lanes; ++l)
        {
            const double lo = static_cast<double>(l) / static_cast<double>(lanes);
            const double hi = static_cast<double>(l + 1) / static_cast<double>(lanes);

            const std::vector<int64_t> upstream   = span(j, lo, hi);
            const std::vector<int64_t> downstream = span(j + 1, lo, hi);

            if (upstream.size() < 2 || downstream.size() < 2)
            {
                message = "segment " + std::to_string(segment) + " lane " +
                          std::to_string(l) + " has no well-defined boundary span";
                return false;
            }

            // Walk the upstream side outward, then the downstream side back:
            // a simple polygon, with extra collinear vertices wherever the
            // neighbor is finer.
            for (int64_t n : upstream)
                definition.faceNodes.push_back(n);

            for (auto it = downstream.rbegin(); it != downstream.rend(); ++it)
                definition.faceNodes.push_back(*it);

            definition.faceNodeOffsets.push_back(
                static_cast<int64_t>(definition.faceNodes.size()));

            const double a = -halfTop + 2.0 * halfTop * lo;
            const double b = -halfTop + 2.0 * halfTop * hi;

            double lowest = m_initialSurface[static_cast<size_t>(segment)];

            for (int k = 0; k < layerCount(); ++k)
            {
                if (lateralOverlap(a, b, 0.5 * width(k, segment), l == 0,
                                   l + 1 == lanes) > 0.0)
                    lowest = std::min(lowest, layerBottom(k));
            }

            bed.push_back(lowest);
        }
    }

    message.clear();
    return true;
}

bool W2Bathymetry::applyPlannedHypsography(Mesh::LayeredMesh &mesh,
                                           const LateralPlan &plan,
                                           std::string &message) const
{
    if (mesh.columnCount() != plan.columnCount)
    {
        message = "the mesh has " + std::to_string(mesh.columnCount()) +
                  " columns but the lateral plan describes " +
                  std::to_string(plan.columnCount);
        return false;
    }

    for (size_t j = 0; j < m_active.size(); ++j)
    {
        const int    segment = m_active[j];
        const int    lanes   = plan.lanes[j];
        const double halfTop = 0.5 * surfaceWidth(segment);
        const double length  = m_length[static_cast<size_t>(segment)];

        for (int l = 0; l < lanes; ++l)
        {
            const double a =
                -halfTop + 2.0 * halfTop * static_cast<double>(l) / lanes;
            const double b =
                -halfTop + 2.0 * halfTop * static_cast<double>(l + 1) / lanes;

            std::vector<double> elevation;
            std::vector<double> area;

            for (int k = layerCount() - 1; k >= 0; --k)
            {
                const double overlap = lateralOverlap(a, b, 0.5 * width(k, segment),
                                                      l == 0, l + 1 == lanes);

                if (overlap <= 0.0)
                    continue;

                const double bottom = layerBottom(k);
                const double middle =
                    bottom + 0.5 * m_layerHeight[static_cast<size_t>(k)];

                const double planArea = overlap * length;

                if (!elevation.empty() && planArea < area.back())
                {
                    elevation.push_back(middle);
                    area.push_back(area.back());
                    continue;
                }

                elevation.push_back(middle);
                area.push_back(planArea);
            }

            if (elevation.size() < 2)
            {
                message = "segment " + std::to_string(segment) + " lane " +
                          std::to_string(l) + " has no storage curve";
                return false;
            }

            if (!mesh.setHypsography(plan.offset[j] + l, elevation, area, message))
                return false;
        }
    }

    message.clear();
    return true;
}

} // namespace FVQual::IO
