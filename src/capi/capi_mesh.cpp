// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "capi_internal.hpp"

#include "fvqual/capi/fvqual_mesh.h"
#include "fvqual/io/ugridreader.h"

#include <memory>
#include <new>

using FVQual::CApi::fail;
using FVQual::CApi::guard;
using FVQual::CApi::succeed;

namespace {

/*!
 * \brief Builds a rectangular channel mesh definition.
 *
 * Kept here rather than in the engine because it is a convenience for
 * scripting and testing, not a modelling capability — real domains come from
 * UGRID files or a mesher.
 */
HydroCouple::SDK::IO::MeshDefinition channelDefinition(int columns, double length,
                                                       double width, double bed)
{
    HydroCouple::SDK::IO::MeshDefinition mesh;
    mesh.meshName = "mesh2d";

    const double dx = length / static_cast<double>(columns);

    for (int i = 0; i <= columns; ++i)
    {
        const double x = static_cast<double>(i) * dx;

        mesh.nodeX.push_back(x);
        mesh.nodeY.push_back(0.0);
        mesh.nodeZ.push_back(bed);

        mesh.nodeX.push_back(x);
        mesh.nodeY.push_back(width);
        mesh.nodeZ.push_back(bed);
    }

    mesh.faceNodeOffsets.push_back(0);

    for (int i = 0; i < columns; ++i)
    {
        mesh.faceNodes.push_back(2 * i);
        mesh.faceNodes.push_back(2 * (i + 1));
        mesh.faceNodes.push_back(2 * (i + 1) + 1);
        mesh.faceNodes.push_back(2 * i + 1);
        mesh.faceNodeOffsets.push_back(static_cast<int64_t>(mesh.faceNodes.size()));
    }

    return mesh;
}

} // namespace

extern "C" {

FVQ_ErrorCode fvq_mesh_create(const double *nodeX, const double *nodeY,
                              const double *nodeZ, int64_t nodeCount,
                              const int64_t *faceOffsets, const int64_t *faceNodes,
                              int64_t faceCount, FVQ_Mesh **out)
{
    return guard([&]() -> FVQ_ErrorCode {
        if (nodeX == nullptr || nodeY == nullptr || faceOffsets == nullptr ||
            faceNodes == nullptr || out == nullptr)
        {
            return fail(FVQ_ERR_NULL_ARGUMENT, "mesh arrays and output must not be null");
        }

        if (nodeCount < 3 || faceCount < 1)
            return fail(FVQ_ERR_INVALID_ARGUMENT, "a mesh needs at least one face");

        HydroCouple::SDK::IO::MeshDefinition definition;
        definition.meshName = "mesh2d";
        definition.nodeX.assign(nodeX, nodeX + nodeCount);
        definition.nodeY.assign(nodeY, nodeY + nodeCount);

        if (nodeZ != nullptr)
            definition.nodeZ.assign(nodeZ, nodeZ + nodeCount);

        definition.faceNodeOffsets.assign(faceOffsets, faceOffsets + faceCount + 1);
        definition.faceNodes.assign(faceNodes, faceNodes + faceOffsets[faceCount]);

        auto        handle = std::make_unique<FVQ_Mesh_>();
        std::string message;

        if (!handle->mesh.build(definition, message))
            return fail(FVQ_ERR_MESH, message);

        *out = handle.release();
        return succeed();
    });
}

FVQ_ErrorCode fvq_mesh_create_channel(int columns, double length, double width,
                                      double bed, FVQ_Mesh **out)
{
    return guard([&]() -> FVQ_ErrorCode {
        if (out == nullptr)
            return fail(FVQ_ERR_NULL_ARGUMENT, "output handle is null");

        if (columns < 1 || !(length > 0.0) || !(width > 0.0))
            return fail(FVQ_ERR_INVALID_ARGUMENT,
                        "channel needs at least one column and positive extents");

        auto        handle = std::make_unique<FVQ_Mesh_>();
        std::string message;

        if (!handle->mesh.build(channelDefinition(columns, length, width, bed), message))
            return fail(FVQ_ERR_MESH, message);

        *out = handle.release();
        return succeed();
    });
}

FVQ_ErrorCode fvq_mesh_read_ugrid(const char *path, const char *meshName, FVQ_Mesh **out)
{
    return guard([&]() -> FVQ_ErrorCode {
        if (path == nullptr || out == nullptr)
            return fail(FVQ_ERR_NULL_ARGUMENT, "path and output must not be null");

        if (!FVQual::IO::UgridReader::available())
            return fail(FVQ_ERR_UNSUPPORTED,
                        "this build has NetCDF support disabled; "
                        "configure with -DFVQual_USE_NETCDF=ON");

        HydroCouple::SDK::IO::MeshDefinition definition;
        std::string                          message;

        if (!FVQual::IO::UgridReader::read(path, meshName == nullptr ? std::string{}
                                                                    : std::string(meshName),
                                           definition, message))
        {
            return fail(FVQ_ERR_FILE, message);
        }

        auto handle = std::make_unique<FVQ_Mesh_>();

        if (!handle->mesh.build(definition, message))
            return fail(FVQ_ERR_MESH, message);

        *out = handle.release();
        return succeed();
    });
}

void fvq_mesh_destroy(FVQ_Mesh *mesh)
{
    delete mesh;
}

int64_t fvq_mesh_node_count(const FVQ_Mesh *mesh)
{
    return mesh == nullptr ? 0 : mesh->mesh.nodeCount();
}

int64_t fvq_mesh_face_count(const FVQ_Mesh *mesh)
{
    return mesh == nullptr ? 0 : mesh->mesh.faceCount();
}

int64_t fvq_mesh_edge_count(const FVQ_Mesh *mesh)
{
    return mesh == nullptr ? 0 : mesh->mesh.edgeCount();
}

double fvq_mesh_total_area(const FVQ_Mesh *mesh)
{
    return mesh == nullptr ? 0.0 : mesh->mesh.totalArea();
}

double fvq_mesh_minimum_spacing(const FVQ_Mesh *mesh)
{
    return mesh == nullptr ? 0.0 : mesh->mesh.minimumCentroidDistance();
}

FVQ_ErrorCode fvq_mesh_face_geometry(const FVQ_Mesh *mesh, double *centroidX,
                                     double *centroidY, double *area, double *bed,
                                     int64_t count)
{
    if (mesh == nullptr)
        return fail(FVQ_ERR_NULL_ARGUMENT, "mesh is null");

    if (count != mesh->mesh.faceCount())
        return fail(FVQ_ERR_INVALID_ARGUMENT,
                    "buffer length must equal the face count");

    const int64_t faces = mesh->mesh.faceCount();

    for (int64_t f = 0; f < faces; ++f)
    {
        const size_t i = static_cast<size_t>(f);

        if (centroidX != nullptr)
            centroidX[i] = mesh->mesh.faceCentroidX()[i];
        if (centroidY != nullptr)
            centroidY[i] = mesh->mesh.faceCentroidY()[i];
        if (area != nullptr)
            area[i] = mesh->mesh.faceArea()[i];
        if (bed != nullptr)
            bed[i] = mesh->mesh.faceBed()[i];
    }

    return succeed();
}

FVQ_ErrorCode fvq_mesh_set_bed(FVQ_Mesh *mesh, const double *bed, int64_t count)
{
    if (mesh == nullptr || bed == nullptr)
        return fail(FVQ_ERR_NULL_ARGUMENT, "mesh and bed array must not be null");

    if (count != mesh->mesh.faceCount())
        return fail(FVQ_ERR_INVALID_ARGUMENT, "bed length must equal the face count");

    for (int64_t f = 0; f < count; ++f)
        mesh->mesh.faceBed()[static_cast<size_t>(f)] = bed[static_cast<size_t>(f)];

    return succeed();
}

} // extern "C"
