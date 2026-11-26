// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/mesh/horizontalmesh.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <utility>

namespace FVQual::Mesh {

namespace {

struct Box {
    double x0, y0, x1, y1;
    bool intersects(const Box &b) const {
        return x0 <= b.x1 && x1 >= b.x0 && y0 <= b.y1 && y1 >= b.y0;
    }
};

// Median-split index bounds the work by local candidates instead of comparing
// every face pair. This matters for the 155,688-cell Bellinge mesh.
class FaceIndex {
    struct Branch { Box box; size_t begin, end; int left = -1, right = -1; };
    const std::vector<Box> &boxes;
    std::vector<size_t> order;
    std::vector<Branch> branches;
    int partition(size_t begin, size_t end) {
        Box bounds = boxes[order[begin]];
        for (size_t i = begin + 1; i < end; ++i) {
            const auto &b = boxes[order[i]];
            bounds.x0 = std::min(bounds.x0, b.x0); bounds.y0 = std::min(bounds.y0, b.y0);
            bounds.x1 = std::max(bounds.x1, b.x1); bounds.y1 = std::max(bounds.y1, b.y1);
        }
        const int slot = static_cast<int>(branches.size());
        branches.push_back({bounds, begin, end});
        if (end - begin > 8) {
            const bool x = bounds.x1 - bounds.x0 >= bounds.y1 - bounds.y0;
            const size_t middle = begin + (end - begin) / 2;
            std::nth_element(order.begin() + begin, order.begin() + middle, order.begin() + end,
                [&](size_t a, size_t b) {
                    const double ca = x ? boxes[a].x0 + boxes[a].x1 : boxes[a].y0 + boxes[a].y1;
                    const double cb = x ? boxes[b].x0 + boxes[b].x1 : boxes[b].y0 + boxes[b].y1;
                    return ca == cb ? a < b : ca < cb;
                });
            const int left = partition(begin, middle), right = partition(middle, end);
            branches[slot].left = left; branches[slot].right = right;
        }
        return slot;
    }
public:
    explicit FaceIndex(const std::vector<Box> &b) : boxes(b), order(b.size()) {
        std::iota(order.begin(), order.end(), 0);
        branches.reserve(2 * b.size());
        partition(0, b.size());
    }
    template<class Visitor> bool visit(const Box &box, Visitor visitor, int slot = 0) const {
        const auto &b = branches[slot];
        if (!box.intersects(b.box)) return true;
        if (b.left >= 0) return visit(box, visitor, b.left) && visit(box, visitor, b.right);
        for (size_t i = b.begin; i < b.end; ++i)
            if (box.intersects(boxes[order[i]]) && !visitor(order[i])) return false;
        return true;
    }
};

bool validateGlobalGeometry(const HorizontalMesh &mesh, std::string &message)
{
    const auto &x = mesh.nodeX(), &y = mesh.nodeY();
    const auto &offsets = mesh.faceNodeOffsets(), &nodes = mesh.faceNodes();
    std::set<int64_t> used(nodes.begin(), nodes.end());
    std::map<std::pair<double, double>, int64_t> coordinates;
    for (auto n : used) {
        auto [it, inserted] = coordinates.emplace(std::pair{x[n], y[n]}, n);
        if (!inserted) {
            message = "distinct used nodes " + std::to_string(it->second) + " and " +
                std::to_string(n) + " have duplicate coordinates; merge and remap explicitly";
            return false;
        }
    }
    auto sign = [&](int64_t a, int64_t b, int64_t c) {
        const long double dx = x[b] - x[a], dy = y[b] - y[a];
        const long double px = x[c] - x[a], py = y[c] - y[a];
        const long double cross = dx * py - dy * px;
        const long double tolerance = 16 * std::numeric_limits<double>::epsilon() *
            (std::abs(dx * py) + std::abs(dy * px));
        return cross > tolerance ? 1 : cross < -tolerance ? -1 : 0;
    };
    auto onSegment = [&](int64_t a, int64_t b, int64_t p) {
        return sign(a, b, p) == 0 && x[p] >= std::min(x[a], x[b]) &&
            x[p] <= std::max(x[a], x[b]) && y[p] >= std::min(y[a], y[b]) &&
            y[p] <= std::max(y[a], y[b]);
    };
    auto inside = [&](int64_t p, size_t f) {
        bool result = false;
        for (auto i = offsets[f]; i < offsets[f+1]; ++i) {
            const auto a = nodes[i], b = nodes[i + 1 == offsets[f+1] ? offsets[f] : i+1];
            if (onSegment(a, b, p)) return false;
            if ((y[a] > y[p]) != (y[b] > y[p]) &&
                static_cast<long double>(x[p]) - x[a] <
                (static_cast<long double>(y[p]) - y[a]) * (x[b] - x[a]) / (y[b] - y[a]))
                result = !result;
        }
        return result;
    };
    std::vector<Box> boxes;
    for (size_t f = 0; f + 1 < offsets.size(); ++f) {
        auto n = nodes[offsets[f]];
        Box box{x[n], y[n], x[n], y[n]};
        for (auto i = offsets[f] + 1; i < offsets[f+1]; ++i) {
            n = nodes[i]; box.x0 = std::min(box.x0, x[n]); box.x1 = std::max(box.x1, x[n]);
            box.y0 = std::min(box.y0, y[n]); box.y1 = std::max(box.y1, y[n]);
        }
        boxes.push_back(box);
    }
    FaceIndex index(boxes);
    for (size_t f = 0; f < boxes.size(); ++f) {
        if (!index.visit(boxes[f], [&](size_t g) {
            if (g <= f) return true;
            auto fail = [&](const char *reason) {
                message = "faces " + std::to_string(f) + " and " + std::to_string(g) + reason;
                return false;
            };
            for (auto i = offsets[f]; i < offsets[f+1]; ++i) {
                const auto a = nodes[i], b = nodes[i+1 == offsets[f+1] ? offsets[f] : i+1];
                for (auto j = offsets[g]; j < offsets[g+1]; ++j) {
                    const auto c = nodes[j], d = nodes[j+1 == offsets[g+1] ? offsets[g] : j+1];
                    if (sign(a,b,c) * sign(a,b,d) < 0 && sign(c,d,a) * sign(c,d,b) < 0)
                        return fail(" overlap (crossing edges)");
                    if ((c != a && c != b && onSegment(a,b,c)) ||
                        (d != a && d != b && onSegment(a,b,d)) ||
                        (a != c && a != d && onSegment(c,d,a)) ||
                        (b != c && b != d && onSegment(c,d,b)))
                        return fail(" have a T-junction or partially overlapping edges");
                }
                if (inside(a, g)) return fail(" overlap (contained vertex)");
            }
            for (auto j = offsets[g]; j < offsets[g+1]; ++j)
                if (inside(nodes[j], f)) return fail(" overlap (contained vertex)");
            return true;
        })) return false;
    }
    // A conforming boundary is a union of closed, non-branching loops. Holes
    // and disconnected domains remain valid and are reported separately.
    std::vector<int> degree(x.size(), 0);
    for (int64_t e = 0; e < mesh.edgeCount(); ++e) if (mesh.isBoundaryEdge(e)) {
        ++degree[mesh.edgeNodes()[2*e]]; ++degree[mesh.edgeNodes()[2*e+1]];
    }
    for (size_t n = 0; n < degree.size(); ++n) if (degree[n] != 0 && degree[n] != 2) {
        message = "boundary branches at node " + std::to_string(n) + " (non-manifold vertex)";
        return false;
    }
    return true;
}

/*!
 * \brief Signed plan-view area of a polygon by the shoelace formula.
 *        Positive when wound counter-clockwise.
 */
double signedArea(const std::vector<double> &x, const std::vector<double> &y,
                  const int64_t *nodes, int64_t count)
{
    double twice = 0.0;
    const double ox = x[static_cast<size_t>(nodes[0])];
    const double oy = y[static_cast<size_t>(nodes[0])];

    for (int64_t i = 0; i < count; ++i)
    {
        const size_t a = static_cast<size_t>(nodes[i]);
        const size_t b = static_cast<size_t>(nodes[(i + 1) % count]);
        twice += (x[a] - ox) * (y[b] - oy) - (x[b] - ox) * (y[a] - oy);
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

    const double ox = x[static_cast<size_t>(nodes[0])];
    const double oy = y[static_cast<size_t>(nodes[0])];
    double sx = 0.0;
    double sy = 0.0;

    for (int64_t i = 0; i < count; ++i)
    {
        const size_t a = static_cast<size_t>(nodes[i]);
        const size_t b = static_cast<size_t>(nodes[(i + 1) % count]);
        const double cross = (x[a] - ox) * (y[b] - oy) - (x[b] - ox) * (y[a] - oy);
        sx += (x[a] + x[b] - 2.0 * ox) * cross;
        sy += (y[a] + y[b] - 2.0 * oy) * cross;
    }

    cx = ox + sx / (6.0 * area);
    cy = oy + sy / (6.0 * area);
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

    for (const auto *values : {&definition.nodeX, &definition.nodeY, &definition.nodeZ})
        for (double value : *values)
            if (!std::isfinite(value)) {
                message = "horizontal mesh has non-finite coordinates or elevations";
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
        std::set<int64_t> unique;
        for (int64_t i = 0; i < count; ++i) {
            if (!unique.insert(nodes[i]).second) {
                message = "face " + std::to_string(f) + " repeats a node";
                return false;
            }
            const size_t a = static_cast<size_t>(nodes[i]);
            const size_t b = static_cast<size_t>(nodes[(i + 1) % count]);
            if (m_nodeX[a] == m_nodeX[b] && m_nodeY[a] == m_nodeY[b]) {
                message = "face " + std::to_string(f) + " has a zero-length edge";
                return false;
            }
        }
        auto orient = [&](int64_t a, int64_t b, int64_t c) {
            return (m_nodeX[b] - m_nodeX[a]) * (m_nodeY[c] - m_nodeY[a]) -
                   (m_nodeY[b] - m_nodeY[a]) * (m_nodeX[c] - m_nodeX[a]);
        };
        for (int64_t i = 0; i < count; ++i)
            for (int64_t j = i + 2; j < count; ++j) {
                if ((j + 1) % count == i) continue;
                const auto a = nodes[i], b = nodes[(i + 1) % count],
                           c = nodes[j], d = nodes[(j + 1) % count];
                const double p = orient(a, b, c), q = orient(a, b, d),
                             r = orient(c, d, a), t = orient(c, d, b);
                auto between = [&](int64_t a, int64_t b, int64_t c) {
                    return m_nodeX[c] >= std::min(m_nodeX[a], m_nodeX[b]) &&
                           m_nodeX[c] <= std::max(m_nodeX[a], m_nodeX[b]) &&
                           m_nodeY[c] >= std::min(m_nodeY[a], m_nodeY[b]) &&
                           m_nodeY[c] <= std::max(m_nodeY[a], m_nodeY[b]);
                };
                if (((p > 0.0) != (q > 0.0) && (r > 0.0) != (t > 0.0) &&
                     p != 0.0 && q != 0.0 && r != 0.0 && t != 0.0) ||
                    (p == 0.0 && between(a, b, c)) || (q == 0.0 && between(a, b, d)) ||
                    (r == 0.0 && between(c, d, a)) || (t == 0.0 && between(c, d, b))) {
                    message = "face " + std::to_string(f) + " is not a simple polygon";
                    return false;
                }
            }
        if (count == 4) {
            const double first = orient(nodes[0], nodes[1], nodes[2]);
            for (int64_t i = 1; i < count; ++i) {
                const double turn = orient(nodes[i], nodes[(i + 1) % count], nodes[(i + 2) % count]);
                if (first == 0.0 || turn == 0.0 || ((turn > 0.0) != (first > 0.0))) {
                    message = "face " + std::to_string(f) + " is not a strictly convex quadrilateral";
                    return false;
                }
            }
        }
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

                if (m_edgeNodes[static_cast<size_t>(2 * edge)] != b ||
                    m_edgeNodes[static_cast<size_t>(2 * edge + 1)] != a) {
                    message = "adjacent faces overlap or traverse a shared edge in the same direction";
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
    m_edgeNormalDistance.assign(static_cast<size_t>(edges), 0.0);

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
        const double nx = m_edgeNormalX[e], ny = m_edgeNormalY[e];
        const double dl = nx * (m_edgeMidX[e] - m_faceCentroidX[left]) +
                          ny * (m_edgeMidY[e] - m_faceCentroidY[left]);
        const double dr = right == kNoFace ? dl :
            nx * (m_faceCentroidX[right] - m_edgeMidX[e]) +
            ny * (m_faceCentroidY[right] - m_edgeMidY[e]);
        m_edgeNormalDistance[e] = right == kNoFace ? dl : dl + dr;
        if (!std::isfinite(dl) || !std::isfinite(dr) || dl <= 0.0 || dr <= 0.0 ||
            !std::isfinite(m_edgeCentroidDistance[e])) {
            message = "edge " + std::to_string(e) + " has a nonpositive centroid-to-face normal distance";
            return false;
        }
    }

    if (!validateGlobalGeometry(*this, message)) return false;
    message.clear();
    return true;
}

GeometryReport HorizontalMesh::geometryReport() const
{
    GeometryReport report;
    if (faceCount() == 0 || edgeCount() == 0) return report;
    report.minimumArea = *std::min_element(m_faceArea.begin(), m_faceArea.end());
    report.maximumArea = *std::max_element(m_faceArea.begin(), m_faceArea.end());
    report.minimumNormalDistance = std::numeric_limits<double>::infinity();
    std::vector<std::vector<int64_t>> boundary(nodeCount()), adjacent(faceCount());
    std::vector<double> closureX(faceCount(), 0.0), closureY(faceCount(), 0.0), perimeter(faceCount(), 0.0);
    for (int64_t f = 0; f < faceCount(); ++f) {
        report.triangles += faceNodeCount(f) == 3; report.quads += faceNodeCount(f) == 4;
        report.polygons += faceNodeCount(f) > 4;
    }
    for (int64_t e = 0; e < edgeCount(); ++e) {
        const auto l = m_edgeLeftFace[e], r = m_edgeRightFace[e];
        const double nx = m_edgeNormalX[e], ny = m_edgeNormalY[e], length = m_edgeLength[e];
        closureX[l] += length * nx; closureY[l] += length * ny; perimeter[l] += length;
        double dx = m_edgeMidX[e] - m_faceCentroidX[l], dy = m_edgeMidY[e] - m_faceCentroidY[l];
        if (r == kNoFace) {
            ++report.boundaryEdges;
            const auto a = m_edgeNodes[2*e], b = m_edgeNodes[2*e+1];
            boundary[a].push_back(b); boundary[b].push_back(a);
        } else {
            closureX[r] -= length * nx; closureY[r] -= length * ny; perimeter[r] += length;
            adjacent[l].push_back(r); adjacent[r].push_back(l);
            dx = m_faceCentroidX[r] - m_faceCentroidX[l]; dy = m_faceCentroidY[r] - m_faceCentroidY[l];
        }
        const double distance = nx * dx + ny * dy;
        report.minimumNormalDistance = std::min(report.minimumNormalDistance, distance);
        report.maximumNonOrthogonalityDegrees = std::max(report.maximumNonOrthogonalityDegrees,
            std::acos(std::clamp(distance / std::hypot(dx,dy), -1.0, 1.0)) * 180.0 / std::acos(-1.0));
    }
    for (int64_t f = 0; f < faceCount(); ++f)
        report.maximumRelativeClosureError = std::max(report.maximumRelativeClosureError,
            std::hypot(closureX[f],closureY[f]) / perimeter[f]);
    auto countComponents = [](const auto &graph, bool omitEmpty) {
        int64_t count = 0; std::vector<bool> seen(graph.size(), false);
        for (size_t start = 0; start < graph.size(); ++start) {
            if (seen[start] || (omitEmpty && graph[start].empty())) continue;
            ++count; std::vector<size_t> stack{start}; seen[start] = true;
            while (!stack.empty()) {
                const auto n = stack.back(); stack.pop_back();
                for (auto next : graph[n]) if (!seen[next]) { seen[next] = true; stack.push_back(next); }
            }
        }
        return count;
    };
    report.boundaryLoops = countComponents(boundary, true);
    report.connectedDomains = countComponents(adjacent, false);
    return report;
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
