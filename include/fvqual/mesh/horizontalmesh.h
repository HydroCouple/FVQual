// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   horizontalmesh.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  HorizontalMesh — the tri-quad horizontal footprint of an FVQual
 *         domain, with finite-volume geometry precomputed (F1.1).
 *
 * FVQual's 3-D mesh is the product of a 2-D horizontal mesh and a vertical
 * layer count: cell (column, k) has index `column * layerCount + k`. This
 * class owns the horizontal half — the part that never changes during a run —
 * and precomputes everything the finite-volume kernels need per face and per
 * edge, once, in structure-of-arrays form.
 *
 * Faces are the mesh polygons (triangles or quadrilaterals); each face is one
 * water column. Edges are unique: an interior edge appears once and knows both
 * adjacent faces, so flux loops visit each edge a single time and scatter the
 * equal-and-opposite contribution to both sides.
 *
 * \threadsafety Immutable after build(); concurrent reads are safe.
 */

#ifndef FVQual_MESH_HORIZONTALMESH_H
#define FVQual_MESH_HORIZONTALMESH_H

#include "fvqual/fvqual_export.h"

#include "hydrocouplesdk/io/meshdefinition.h"

#include <cstdint>
#include <string>
#include <vector>

namespace FVQual::Mesh {

/*!
 * \brief Sentinel for the outward side of a boundary edge.
 */
inline constexpr int64_t kNoFace = -1;

/*!
 * \brief Tri-quad horizontal mesh with precomputed finite-volume geometry.
 *
 * All arrays are structure-of-arrays and allocation-stable after build(), so
 * they can be mirrored to a device buffer without repacking.
 */
class FVQual_EXPORT HorizontalMesh
{
public:
    HorizontalMesh() = default;

    /*!
     * \brief Builds the mesh from a UGRID-congruent mesh definition.
     *
     * Validates the definition, extracts unique edges, and precomputes face
     * areas/centroids, edge lengths/normals/midpoints, and the face-neighbour
     * table. Faces must be simple polygons of 3 or more nodes wound
     * consistently; winding is normalized to counter-clockwise (positive
     * signed area) so edge normals point out of the left face.
     *
     * \param[in]  definition Horizontal mesh (nodes + CSR face connectivity).
     * \param[out] message    Failure description when returning false.
     * \return true when the mesh was built.
     */
    [[nodiscard]] bool build(const HydroCouple::SDK::IO::MeshDefinition &definition,
                             std::string &message);

    // ── Counts ────────────────────────────────────────────────────────────

    [[nodiscard]] int64_t nodeCount() const noexcept
    {
        return static_cast<int64_t>(m_nodeX.size());
    }

    /*!
     * \brief Number of faces — equivalently, the number of water columns.
     */
    [[nodiscard]] int64_t faceCount() const noexcept
    {
        return m_faceNodeOffsets.empty()
                   ? 0
                   : static_cast<int64_t>(m_faceNodeOffsets.size()) - 1;
    }

    [[nodiscard]] int64_t edgeCount() const noexcept
    {
        return static_cast<int64_t>(m_edgeLength.size());
    }

    // ── Node coordinates ──────────────────────────────────────────────────

    [[nodiscard]] const std::vector<double> &nodeX() const noexcept { return m_nodeX; }
    [[nodiscard]] const std::vector<double> &nodeY() const noexcept { return m_nodeY; }

    /*!
     * \brief Bed elevation per node (positive up, metres). Populated from the
     *        mesh definition's nodeZ when present, else zero-filled.
     */
    [[nodiscard]] const std::vector<double> &nodeBed() const noexcept { return m_nodeBed; }
    [[nodiscard]] std::vector<double> &nodeBed() noexcept { return m_nodeBed; }

    // ── Face connectivity and geometry ────────────────────────────────────

    [[nodiscard]] const std::vector<int64_t> &faceNodeOffsets() const noexcept
    {
        return m_faceNodeOffsets;
    }
    [[nodiscard]] const std::vector<int64_t> &faceNodes() const noexcept
    {
        return m_faceNodes;
    }

    /*!
     * \brief Node count of \p face.
     */
    [[nodiscard]] int64_t faceNodeCount(int64_t face) const noexcept
    {
        return m_faceNodeOffsets[static_cast<size_t>(face) + 1] -
               m_faceNodeOffsets[static_cast<size_t>(face)];
    }

    /*!
     * \brief Plan-view area of each face (m²).
     */
    [[nodiscard]] const std::vector<double> &faceArea() const noexcept { return m_faceArea; }

    [[nodiscard]] const std::vector<double> &faceCentroidX() const noexcept
    {
        return m_faceCentroidX;
    }
    [[nodiscard]] const std::vector<double> &faceCentroidY() const noexcept
    {
        return m_faceCentroidY;
    }

    /*!
     * \brief Bed elevation at each face centroid (m, positive up) — the
     *        column's bed, averaged from its nodes.
     */
    [[nodiscard]] const std::vector<double> &faceBed() const noexcept { return m_faceBed; }
    [[nodiscard]] std::vector<double> &faceBed() noexcept { return m_faceBed; }

    // ── Unique-edge geometry ──────────────────────────────────────────────

    /*!
     * \brief Endpoint node indexes of each edge (2 per edge, flattened).
     */
    [[nodiscard]] const std::vector<int64_t> &edgeNodes() const noexcept { return m_edgeNodes; }

    /*!
     * \brief Face on the left of each edge — always a real face.
     */
    [[nodiscard]] const std::vector<int64_t> &edgeLeftFace() const noexcept
    {
        return m_edgeLeftFace;
    }

    /*!
     * \brief Face on the right of each edge, or kNoFace on the boundary.
     */
    [[nodiscard]] const std::vector<int64_t> &edgeRightFace() const noexcept
    {
        return m_edgeRightFace;
    }

    [[nodiscard]] const std::vector<double> &edgeLength() const noexcept { return m_edgeLength; }

    /*!
     * \brief Unit normal pointing from the left face toward the right face.
     */
    [[nodiscard]] const std::vector<double> &edgeNormalX() const noexcept
    {
        return m_edgeNormalX;
    }
    [[nodiscard]] const std::vector<double> &edgeNormalY() const noexcept
    {
        return m_edgeNormalY;
    }

    [[nodiscard]] const std::vector<double> &edgeMidX() const noexcept { return m_edgeMidX; }
    [[nodiscard]] const std::vector<double> &edgeMidY() const noexcept { return m_edgeMidY; }

    /*!
     * \brief Distance between the centroids either side of an edge; for a
     *        boundary edge, centroid-to-edge-midpoint. Used as the gradient
     *        denominator in two-point flux approximations.
     */
    [[nodiscard]] const std::vector<double> &edgeCentroidDistance() const noexcept
    {
        return m_edgeCentroidDistance;
    }

    /*!
     * \brief True when the edge has no right-hand face.
     */
    [[nodiscard]] bool isBoundaryEdge(int64_t edge) const noexcept
    {
        return m_edgeRightFace[static_cast<size_t>(edge)] == kNoFace;
    }

    /*!
     * \brief Indexes of the edges bounding \p face (CSR: offsets then values).
     */
    [[nodiscard]] const std::vector<int64_t> &faceEdgeOffsets() const noexcept
    {
        return m_faceEdgeOffsets;
    }
    [[nodiscard]] const std::vector<int64_t> &faceEdges() const noexcept { return m_faceEdges; }

    /*!
     * \brief Sum of all face areas (m²) — the domain's plan-view footprint.
     */
    [[nodiscard]] double totalArea() const noexcept;

    /*!
     * \brief Shortest distance between adjacent face centroids (m); the
     *        length scale that sets the horizontal CFL limit.
     */
    [[nodiscard]] double minimumCentroidDistance() const noexcept;

private:
    std::vector<double> m_nodeX;
    std::vector<double> m_nodeY;
    std::vector<double> m_nodeBed;

    std::vector<int64_t> m_faceNodeOffsets;
    std::vector<int64_t> m_faceNodes;
    std::vector<double>  m_faceArea;
    std::vector<double>  m_faceCentroidX;
    std::vector<double>  m_faceCentroidY;
    std::vector<double>  m_faceBed;

    std::vector<int64_t> m_edgeNodes;
    std::vector<int64_t> m_edgeLeftFace;
    std::vector<int64_t> m_edgeRightFace;
    std::vector<double>  m_edgeLength;
    std::vector<double>  m_edgeNormalX;
    std::vector<double>  m_edgeNormalY;
    std::vector<double>  m_edgeMidX;
    std::vector<double>  m_edgeMidY;
    std::vector<double>  m_edgeCentroidDistance;

    std::vector<int64_t> m_faceEdgeOffsets;
    std::vector<int64_t> m_faceEdges;
};

} // namespace FVQual::Mesh

#endif // FVQual_MESH_HORIZONTALMESH_H
