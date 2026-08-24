// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/mesh/horizontalmesh.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <utility>

namespace FVQual::Mesh {

namespace {

/*!
 * \brief Signed plan-view area of a polygon by the shoelace formula.
 *        Positive when wound counter-clockwise.
 */
double signedArea(const std::vector<double> &x, const std::vector<double> &y,
                  const int64_t *nodes, int64_t count)
{
    double twice = 0.0;

    for (int64_t i = 0; i < count; ++i)
    {
        const size_t a = static_cast<size_t>(nodes[i]);
        const size_t b = static_cast<size_t>(nodes[(i + 1) % count]);
        twice += x[a] * y[b] - x[b] * y[a];
    }

    return 0.5 * twice;
}

/*!
 * \brief Area-weighted centroid of a simple polygon. Falls back to the vertex
 *        mean for degenerate (zero-area) polygons, which validation rejects
 *        anyway but which must not divide by zero first.
 */
void polygonCentroid(const std::vector<double> &x, const std::vector<double> &y,
                     const int64_t *nodes, int64_t count, double area,
                     double &cx, double &cy)
{
    if (std::abs(area) < 1e-30)
    {
        cx = 0.0;
        cy = 0.0;

        for (int64_t i = 0; i < count; ++i)
        {
            cx += x[static_cast<size_t>(nodes[i])];
            cy += y[static_cast<size_t>(nodes[i])];
        }

        cx /= static_cast<double>(count);
        cy /= static_cast<double>(count);
        return;
    }

    double sx = 0.0;
    double sy = 0.0;

    for (int64_t i = 0; i < count; ++i)
    {
        const size_t a = static_cast<size_t>(nodes[i]);
        const size_t b = static_cast<size_t>(nodes[(i + 1) % count]);
        const double cross = x[a] * y[b] - x[b] * y[a];
        sx += (x[a] + x[b]) * cross;
        sy += (y[a] + y[b]) * cross;
    }

    cx = sx / (6.0 * area);
    cy = sy / (6.0 * area);
}

} // namespace

bool HorizontalMesh::build(const HydroCouple::SDK::IO::MeshDefinition &definition,
                           std::string &message)
{
    if (!definition.validate(message))
        return false;

    if (definition.faceCount() <= 0)
    {
        message = "horizontal mesh has no faces";
        return false;
    }

    m_nodeX = definition.nodeX;
    m_nodeY = definition.nodeY;

    m_nodeBed.assign(m_nodeX.size(), 0.0);
    if (definition.nodeZ.size() == m_nodeX.size())
        m_nodeBed = definition.nodeZ;

    m_faceNodeOffsets = definition.faceNodeOffsets;
    m_faceNodes       = definition.faceNodes;

    const int64_t faces = faceCount();

    m_faceArea.assign(static_cast<size_t>(faces), 0.0);
    m_faceCentroidX.assign(static_cast<size_t>(faces), 0.0);
    m_faceCentroidY.assign(static_cast<size_t>(faces), 0.0);
    m_faceBed.assign(static_cast<size_t>(faces), 0.0);

    // ── Face geometry; normalize winding to counter-clockwise ─────────────
    for (int64_t f = 0; f < faces; ++f)
    {
        const int64_t begin = m_faceNodeOffsets[static_cast<size_t>(f)];
        const int64_t count = faceNodeCount(f);

        if (count < 3)
        {
            message = "face " + std::to_string(f) + " has fewer than 3 nodes";
            return false;
        }

        int64_t *nodes = m_faceNodes.data() + begin;
        double   area  = signedArea(m_nodeX, m_nodeY, nodes, count);

        if (area < 0.0)
        {
            std::reverse(nodes, nodes + count);
            area = -area;
        }

        if (area <= 0.0)
        {
            message = "face " + std::to_string(f) +
                      " is degenerate (zero or ill-defined area)";
            return false;
        }

        m_faceArea[static_cast<size_t>(f)] = area;

        polygonCentroid(m_nodeX, m_nodeY, nodes, count, area,
                        m_faceCentroidX[static_cast<size_t>(f)],
                        m_faceCentroidY[static_cast<size_t>(f)]);

        double bed = 0.0;
        for (int64_t i = 0; i < count; ++i)
            bed += m_nodeBed[static_cast<size_t>(nodes[i])];

        m_faceBed[static_cast<size_t>(f)] = bed / static_cast<double>(count);
    }

    // ── Unique edges ──────────────────────────────────────────────────────
    //
    // Keyed on the sorted endpoint pair. The first face to claim an edge
    // becomes its left face and fixes the normal's sign; a second claimant
    // becomes the right face. A third is a non-manifold error.
    std::map<std::pair<int64_t, int64_t>, int64_t> edgeIndex;

    m_edgeNodes.clear();
    m_edgeLeftFace.clear();
    m_edgeRightFace.clear();

    m_faceEdgeOffsets.assign(static_cast<size_t>(faces) + 1, 0);
    m_faceEdges.clear();
    m_faceEdges.reserve(m_faceNodes.size());

    for (int64_t f = 0; f < faces; ++f)
    {
        const int64_t        begin = m_faceNodeOffsets[static_cast<size_t>(f)];
        const int64_t        count = faceNodeCount(f);
        const int64_t *const nodes = m_faceNodes.data() + begin;

        for (int64_t i = 0; i < count; ++i)
        {
            const int64_t a = nodes[i];
            const int64_t b = nodes[(i + 1) % count];
            const auto    key = std::minmax(a, b);

            auto it = edgeIndex.find({key.first, key.second});

            if (it == edgeIndex.end())
            {
                const int64_t edge = static_cast<int64_t>(m_edgeLeftFace.size());
                edgeIndex.emplace(std::pair{key.first, key.second}, edge);

                m_edgeNodes.push_back(a);
                m_edgeNodes.push_back(b);
                m_edgeLeftFace.push_back(f);
                m_edgeRightFace.push_back(kNoFace);
                m_faceEdges.push_back(edge);
            }
            else
            {
                const int64_t edge = it->second;

                if (m_edgeRightFace[static_cast<size_t>(edge)] != kNoFace)
                {
                    message = "edge (" + std::to_string(key.first) + ", " +
                              std::to_string(key.second) +
                              ") is shared by more than two faces (non-manifold mesh)";
                    return false;
                }

                m_edgeRightFace[static_cast<size_t>(edge)] = f;
                m_faceEdges.push_back(edge);
            }
        }

        m_faceEdgeOffsets[static_cast<size_t>(f) + 1] =
            static_cast<int64_t>(m_faceEdges.size());
    }

    // ── Edge geometry ─────────────────────────────────────────────────────
    const int64_t edges = static_cast<int64_t>(m_edgeLeftFace.size());

    m_edgeLength.assign(static_cast<size_t>(edges), 0.0);
    m_edgeNormalX.assign(static_cast<size_t>(edges), 0.0);
    m_edgeNormalY.assign(static_cast<size_t>(edges), 0.0);
    m_edgeMidX.assign(static_cast<size_t>(edges), 0.0);
    m_edgeMidY.assign(static_cast<size_t>(edges), 0.0);
    m_edgeCentroidDistance.assign(static_cast<size_t>(edges), 0.0);

    for (int64_t e = 0; e < edges; ++e)
    {
        const size_t a = static_cast<size_t>(m_edgeNodes[static_cast<size_t>(e) * 2]);
        const size_t b = static_cast<size_t>(m_edgeNodes[static_cast<size_t>(e) * 2 + 1]);

        const double dx     = m_nodeX[b] - m_nodeX[a];
        const double dy     = m_nodeY[b] - m_nodeY[a];
        const double length = std::hypot(dx, dy);

        if (length <= 0.0)
        {
            message = "edge " + std::to_string(e) + " has zero length";
            return false;
        }

        m_edgeLength[static_cast<size_t>(e)] = length;
        m_edgeMidX[static_cast<size_t>(e)]   = 0.5 * (m_nodeX[a] + m_nodeX[b]);
        m_edgeMidY[static_cast<size_t>(e)]   = 0.5 * (m_nodeY[a] + m_nodeY[b]);

        // With counter-clockwise winding the outward normal of the left face
        // is the edge tangent rotated -90 degrees.
        m_edgeNormalX[static_cast<size_t>(e)] = dy / length;
        m_edgeNormalY[static_cast<size_t>(e)] = -dx / length;

        const size_t left  = static_cast<size_t>(m_edgeLeftFace[static_cast<size_t>(e)]);
        const int64_t right = m_edgeRightFace[static_cast<size_t>(e)];

        if (right == kNoFace)
        {
            m_edgeCentroidDistance[static_cast<size_t>(e)] =
                std::hypot(m_edgeMidX[static_cast<size_t>(e)] - m_faceCentroidX[left],
                           m_edgeMidY[static_cast<size_t>(e)] - m_faceCentroidY[left]);
        }
        else
        {
            const size_t r = static_cast<size_t>(right);
            m_edgeCentroidDistance[static_cast<size_t>(e)] =
                std::hypot(m_faceCentroidX[r] - m_faceCentroidX[left],
                           m_faceCentroidY[r] - m_faceCentroidY[left]);
        }
    }

    message.clear();
    return true;
}

double HorizontalMesh::totalArea() const noexcept
{
    double total = 0.0;

    for (double area : m_faceArea)
        total += area;

    return total;
}

double HorizontalMesh::minimumCentroidDistance() const noexcept
{
    double smallest = std::numeric_limits<double>::infinity();

    for (int64_t e = 0; e < edgeCount(); ++e)
    {
        if (isBoundaryEdge(e))
            continue;

        smallest = std::min(smallest, m_edgeCentroidDistance[static_cast<size_t>(e)]);
    }

    return std::isinf(smallest) ? 0.0 : smallest;
}

} // namespace FVQual::Mesh
