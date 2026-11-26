/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Caleb Buahin
 *
 * fvqual_mesh.h — horizontal mesh construction and geometry queries.
 *
 * A mesh is built once, queried freely, and then handed to one or more
 * simulations. Simulations take a private copy of it, so a mesh may be
 * destroyed or modified afterwards without disturbing a running model — the C
 * API deliberately avoids lifetime coupling between handles, because a
 * dangling reference is the one error class a binding cannot recover from.
 */

#ifndef FVQual_CAPI_MESH_H
#define FVQual_CAPI_MESH_H

#include "fvqual/capi/fvqual_engine.h"

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief Opaque handle to a horizontal tri-quad mesh.
 */
typedef struct FVQ_Mesh_ FVQ_Mesh;

/* ── Construction ──────────────────────────────────────────────────────── */

/*!
 * \brief Builds a mesh from node coordinates and CSR face connectivity.
 *
 * Faces may mix triangles and quadrilaterals. Winding is normalized
 * internally, so either orientation is accepted.
 *
 * \param[in]  nodeX        Node x coordinates (m), length \p nodeCount.
 * \param[in]  nodeY        Node y coordinates (m), length \p nodeCount.
 * \param[in]  nodeZ        Bed elevation at nodes (m, positive up); NULL for a
 *                          flat bed at zero, which fvq_mesh_set_bed can then
 *                          overwrite per face.
 * \param[in]  nodeCount    Number of nodes.
 * \param[in]  faceOffsets  CSR offsets, length \p faceCount + 1.
 * \param[in]  faceNodes    Concatenated node indexes.
 * \param[in]  faceCount    Number of faces (water columns).
 * \param[out] out          Receives the new mesh; destroy with fvq_mesh_destroy.
 */
FVQual_CAPI FVQ_ErrorCode fvq_mesh_create(const double *nodeX, const double *nodeY,
                                          const double *nodeZ, int64_t nodeCount,
                                          const int64_t *faceOffsets,
                                          const int64_t *faceNodes, int64_t faceCount,
                                          FVQ_Mesh **out);

/*!
 * \brief Builds a straight rectangular channel of quadrilateral cells.
 *
 * The workhorse of the verification suite and the fastest way to get a model
 * running from a script.
 *
 * \param[in]  columns Number of cells along the channel (>= 1).
 * \param[in]  length  Channel length (m).
 * \param[in]  width   Channel width (m).
 * \param[in]  bed     Uniform bed elevation (m, positive up; negative is below datum).
 * \param[out] out     Receives the new mesh.
 */
FVQual_CAPI FVQ_ErrorCode fvq_mesh_create_channel(int columns, double length, double width,
                                                  double bed, FVQ_Mesh **out);

/*!
 * \brief Reads a UGRID-1.0 two-dimensional mesh from a NetCDF-4 file.
 *
 * Returns FVQ_ERR_UNSUPPORTED when built without NetCDF.
 *
 * \param[in]  path     File to read.
 * \param[in]  meshName Mesh topology variable, or NULL for the first 2-D mesh.
 * \param[out] out      Receives the new mesh.
 */
FVQual_CAPI FVQ_ErrorCode fvq_mesh_read_ugrid(const char *path, const char *meshName,
                                              FVQ_Mesh **out);

/*! Reads a versioned SI mesh package with declared bed elevations.
 * Terrain elevations and unapplied non-wall boundaries are refused.
 * Initial depth and roughness are carried into simulations created from it.
 */
FVQual_CAPI FVQ_ErrorCode fvq_mesh_read_package(const char *path, FVQ_Mesh **out);

/*!
 * \brief Releases a mesh. Passing NULL is a no-op.
 */
FVQual_CAPI void fvq_mesh_destroy(FVQ_Mesh *mesh);

/* ── Geometry ──────────────────────────────────────────────────────────── */

FVQual_CAPI int64_t fvq_mesh_node_count(const FVQ_Mesh *mesh);
FVQual_CAPI int64_t fvq_mesh_face_count(const FVQ_Mesh *mesh);
FVQual_CAPI int64_t fvq_mesh_edge_count(const FVQ_Mesh *mesh);

/*!
 * \brief Total plan-view area of the mesh (m²).
 */
FVQual_CAPI double fvq_mesh_total_area(const FVQ_Mesh *mesh);

/*!
 * \brief Shortest distance between adjacent face centroids (m) — the length
 *        scale that sets the horizontal Courant limit.
 */
FVQual_CAPI double fvq_mesh_minimum_spacing(const FVQ_Mesh *mesh);

/*!
 * \brief Copies per-face geometry into caller-provided buffers.
 *
 * Any output pointer may be NULL to skip that quantity. \p count must equal
 * fvq_mesh_face_count().
 *
 * \param[out] centroidX Face centroid x (m).
 * \param[out] centroidY Face centroid y (m).
 * \param[out] area      Face plan area (m²).
 * \param[out] bed       Face bed elevation (m).
 */
FVQual_CAPI FVQ_ErrorCode fvq_mesh_face_geometry(const FVQ_Mesh *mesh, double *centroidX,
                                                 double *centroidY, double *area,
                                                 double *bed, int64_t count);

/*!
 * \brief Overwrites the bed elevation of every face (m, positive up).
 *
 * Bathymetry usually arrives from a raster sample or a survey rather than from
 * the mesh file, so it is settable independently of construction.
 *
 * \param[in] bed   One elevation per face.
 * \param[in] count Must equal fvq_mesh_face_count().
 */
FVQual_CAPI FVQ_ErrorCode fvq_mesh_set_bed(FVQ_Mesh *mesh, const double *bed,
                                           int64_t count);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* FVQual_CAPI_MESH_H */
