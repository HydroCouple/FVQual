// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/io/layeredugridwriter.h"

#include <algorithm>
#include <cstdint>
#include <vector>

#ifdef FVQual_HAVE_NETCDF
#include <netcdf.h>
#endif

namespace FVQual::IO {

namespace SDKIO = HydroCouple::SDK::IO;

namespace {

constexpr const char *kTopologyVar     = "mesh2d";
constexpr const char *kNodeXVar        = "mesh2d_node_x";
constexpr const char *kNodeYVar        = "mesh2d_node_y";
constexpr const char *kNodeZVar        = "mesh2d_node_z";
constexpr const char *kFaceNodesVar    = "mesh2d_face_nodes";
constexpr const char *kSigmaLayerVar   = "mesh2d_layer_sigma";
constexpr const char *kSigmaIfaceVar   = "mesh2d_interface_sigma";
constexpr const char *kDepthVar        = "mesh2d_depth";
constexpr const char *kTimeVar         = "time";

constexpr int64_t kConnectivityFill = -1;

/*!
 * \brief One registered output field.
 */
struct FieldSpec
{
    std::string itemId;
    std::string units;
    std::string standardName;
    bool        layered = false;
    int         varid   = -1;
};

} // namespace

// ── Implementation ────────────────────────────────────────────────────────

struct LayeredUgridWriter::Impl
{
    std::filesystem::path  path;
    SDKIO::MeshDefinition  mesh;
    LayeredUgridMetadata   vertical;
    std::vector<FieldSpec> fields;
    std::string            surfaceFieldId = "eta";
    std::string            error;
    int64_t                timeIndex = 0;

#ifdef FVQual_HAVE_NETCDF
    int ncid      = -1;
    int timeVarid = -1;
    std::vector<double> staging;

    /*!
     * \brief Records a NetCDF failure and returns its status code.
     */
    int fail(int status, const std::string &context)
    {
        error = context + ": " + nc_strerror(status);
        return status;
    }
#endif
};

LayeredUgridWriter::LayeredUgridWriter(std::filesystem::path path,
                                       SDKIO::MeshDefinition mesh,
                                       LayeredUgridMetadata vertical)
    : m_impl(std::make_unique<Impl>())
{
    m_impl->path     = std::move(path);
    m_impl->mesh     = std::move(mesh);
    m_impl->vertical = std::move(vertical);
}

LayeredUgridWriter::~LayeredUgridWriter()
{
#ifdef FVQual_HAVE_NETCDF
    if (m_impl && m_impl->ncid >= 0)
        nc_close(m_impl->ncid);
#endif
}

void LayeredUgridWriter::addLayeredField(const std::string &itemId,
                                         const std::string &units,
                                         const std::string &standardName)
{
    m_impl->fields.push_back(FieldSpec{itemId, units, standardName, true, -1});
}

void LayeredUgridWriter::addColumnField(const std::string &itemId,
                                        const std::string &units,
                                        const std::string &standardName)
{
    m_impl->fields.push_back(FieldSpec{itemId, units, standardName, false, -1});
}

void LayeredUgridWriter::setSurfaceFieldId(const std::string &itemId)
{
    m_impl->surfaceFieldId = itemId;
}

std::string LayeredUgridWriter::name() const
{
    return "layered_ugrid_netcdf";
}

std::string LayeredUgridWriter::errorMessage() const
{
    return m_impl->error;
}

const std::filesystem::path &LayeredUgridWriter::path() const noexcept
{
    return m_impl->path;
}

std::vector<SDKIO::ResultEntry> LayeredUgridWriter::catalog() const
{
    std::vector<SDKIO::ResultEntry> entries;
    entries.reserve(m_impl->fields.size());

    const int64_t faces  = m_impl->mesh.faceCount();
    const int64_t layers = m_impl->vertical.layerCount();

    for (const FieldSpec &field : m_impl->fields)
    {
        SDKIO::ResultEntry entry;
        entry.itemId   = field.itemId;
        entry.artifact = m_impl->path.filename().string();
        entry.format   = "netcdf_ugrid";
        entry.variable = field.itemId;
        entry.units    = field.units;
        entry.mesh     = m_impl->mesh.meshName;
        entry.location = SDKIO::MeshLocation::Face;

        if (field.layered)
        {
            entry.shape      = {m_impl->timeIndex, faces, layers};
            entry.dimensions = {"time", "mesh2d_nFaces", "mesh2d_nLayers"};
        }
        else
        {
            entry.shape      = {m_impl->timeIndex, faces};
            entry.dimensions = {"time", "mesh2d_nFaces"};
        }

        entry.time = SDKIO::ResultTimeAxis{};
        entries.push_back(std::move(entry));
    }

    return entries;
}

bool LayeredUgridWriter::available() noexcept
{
#ifdef FVQual_HAVE_NETCDF
    return true;
#else
    return false;
#endif
}

#ifdef FVQual_HAVE_NETCDF

int LayeredUgridWriter::initialize()
{
    Impl &impl = *m_impl;

    if (!impl.mesh.validate(impl.error))
        return NC_EINVAL;

    const int layers = impl.vertical.layerCount();

    if (layers < 1)
    {
        impl.error = "layered writer needs at least two sigma interfaces";
        return NC_EINVAL;
    }

    const int64_t faces = impl.mesh.faceCount();

    if (static_cast<int64_t>(impl.vertical.bedElevation.size()) != faces)
    {
        impl.error = "bed elevation array has " +
                     std::to_string(impl.vertical.bedElevation.size()) +
                     " entries; expected one per face (" + std::to_string(faces) + ")";
        return NC_EINVAL;
    }

    const bool hasSurfaceField =
        std::any_of(impl.fields.begin(), impl.fields.end(),
                    [&](const FieldSpec &f) {
                        return !f.layered && f.itemId == impl.surfaceFieldId;
                    });

    if (!hasSurfaceField)
    {
        impl.error = "the CF sigma formula needs a free-surface column field named '" +
                     impl.surfaceFieldId +
                     "'; register it with addColumnField() or call setSurfaceFieldId()";
        return NC_EINVAL;
    }

    std::error_code ec;
    if (impl.path.has_parent_path())
        std::filesystem::create_directories(impl.path.parent_path(), ec);

    int status = nc_create(impl.path.string().c_str(), NC_NETCDF4 | NC_CLOBBER, &impl.ncid);
    if (status != NC_NOERR)
        return impl.fail(status, "creating '" + impl.path.string() + "'");

    const int ncid = impl.ncid;

    // ── Dimensions ────────────────────────────────────────────────────────
    const int64_t nodes        = impl.mesh.nodeCount();
    const int64_t maxFaceNodes = impl.mesh.maxNodesPerFace();

    int dimNode = -1, dimFace = -1, dimFaceNodes = -1;
    int dimLayer = -1, dimInterface = -1, dimTime = -1;

    if ((status = nc_def_dim(ncid, "mesh2d_nNodes", static_cast<size_t>(nodes), &dimNode)) ||
        (status = nc_def_dim(ncid, "mesh2d_nFaces", static_cast<size_t>(faces), &dimFace)) ||
        (status = nc_def_dim(ncid, "mesh2d_nMaxFaceNodes",
                             static_cast<size_t>(maxFaceNodes), &dimFaceNodes)) ||
        (status = nc_def_dim(ncid, "mesh2d_nLayers", static_cast<size_t>(layers),
                             &dimLayer)) ||
        (status = nc_def_dim(ncid, "mesh2d_nInterfaces", static_cast<size_t>(layers + 1),
                             &dimInterface)) ||
        (status = nc_def_dim(ncid, "time", NC_UNLIMITED, &dimTime)))
    {
        return impl.fail(status, "defining dimensions");
    }

    // ── Global attributes ─────────────────────────────────────────────────
    const std::string conventions = "CF-1.11 UGRID-1.0";
    const std::string source      = "FVQual";

    if ((status = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", conventions.size(),
                                  conventions.c_str())) ||
        (status = nc_put_att_text(ncid, NC_GLOBAL, "source", source.size(),
                                  source.c_str())))
    {
        return impl.fail(status, "writing global attributes");
    }

    // ── Mesh topology dummy variable ──────────────────────────────────────
    int topologyVar = -1;
    if ((status = nc_def_var(ncid, kTopologyVar, NC_INT, 0, nullptr, &topologyVar)))
        return impl.fail(status, "defining mesh topology variable");

    // Carrying node elevations makes the file a complete model input, not
    // just a results archive: FVQual can re-read its own output and rebuild
    // the domain. UGRID permits a third node coordinate.
    const bool hasNodeZ =
        impl.mesh.nodeZ.size() == static_cast<size_t>(impl.mesh.nodeCount());

    const std::string role = "mesh_topology";
    const std::string nodeCoordinates =
        std::string(kNodeXVar) + " " + kNodeYVar + (hasNodeZ ? std::string(" ") + kNodeZVar : "");
    const std::string longName    = "Topology of the 2D unstructured mesh";
    const int         topologyDim = 2;

    if ((status = nc_put_att_text(ncid, topologyVar, "cf_role", role.size(),
                                  role.c_str())) ||
        (status = nc_put_att_int(ncid, topologyVar, "topology_dimension", NC_INT, 1,
                                 &topologyDim)) ||
        (status = nc_put_att_text(ncid, topologyVar, "node_coordinates",
                                  nodeCoordinates.size(), nodeCoordinates.c_str())) ||
        (status = nc_put_att_text(ncid, topologyVar, "face_node_connectivity",
                                  std::string(kFaceNodesVar).size(), kFaceNodesVar)) ||
        (status = nc_put_att_text(ncid, topologyVar, "long_name", longName.size(),
                                  longName.c_str())))
    {
        return impl.fail(status, "writing mesh topology attributes");
    }

    // ── Node coordinates ──────────────────────────────────────────────────
    int nodeXVar = -1, nodeYVar = -1, nodeZVar = -1;
    if ((status = nc_def_var(ncid, kNodeXVar, NC_DOUBLE, 1, &dimNode, &nodeXVar)) ||
        (status = nc_def_var(ncid, kNodeYVar, NC_DOUBLE, 1, &dimNode, &nodeYVar)))
    {
        return impl.fail(status, "defining node coordinate variables");
    }

    const std::string metres         = "m";
    const std::string positiveUpText = "up";

    for (int var : {nodeXVar, nodeYVar})
    {
        if ((status = nc_put_att_text(ncid, var, "units", metres.size(), metres.c_str())))
            return impl.fail(status, "writing node coordinate units");
    }

    if (hasNodeZ)
    {
        const std::string bedLongName = "bed elevation at mesh nodes";

        if ((status = nc_def_var(ncid, kNodeZVar, NC_DOUBLE, 1, &dimNode, &nodeZVar)) ||
            (status = nc_put_att_text(ncid, nodeZVar, "units", metres.size(),
                                      metres.c_str())) ||
            (status = nc_put_att_text(ncid, nodeZVar, "positive", positiveUpText.size(),
                                      positiveUpText.c_str())) ||
            (status = nc_put_att_text(ncid, nodeZVar, "long_name", bedLongName.size(),
                                      bedLongName.c_str())))
        {
            return impl.fail(status, "defining node elevations");
        }
    }

    // ── Face-node connectivity ────────────────────────────────────────────
    int       faceNodesVar = -1;
    const int connectivityDims[2] = {dimFace, dimFaceNodes};

    if ((status = nc_def_var(ncid, kFaceNodesVar, NC_INT64, 2, connectivityDims,
                             &faceNodesVar)))
    {
        return impl.fail(status, "defining face-node connectivity");
    }

    const std::string connectivityRole = "face_node_connectivity";
    const long long   startIndex       = 0;
    const long long   fillValue        = kConnectivityFill;

    if ((status = nc_put_att_text(ncid, faceNodesVar, "cf_role", connectivityRole.size(),
                                  connectivityRole.c_str())) ||
        (status = nc_put_att_longlong(ncid, faceNodesVar, "start_index", NC_INT64, 1,
                                      &startIndex)) ||
        (status = nc_put_att_longlong(ncid, faceNodesVar, "_FillValue", NC_INT64, 1,
                                      &fillValue)))
    {
        return impl.fail(status, "writing connectivity attributes");
    }

    // ── Vertical coordinate ───────────────────────────────────────────────
    //
    // CF ocean_sigma_coordinate runs 0 at the surface to -1 at the bed and
    // pairs with a positive-down depth; FVQual's internal sigma runs 0 to 1
    // downward against positive-up elevations. Negate on the way out.
    int sigmaLayerVar = -1, sigmaInterfaceVar = -1, depthVar = -1;

    if ((status = nc_def_var(ncid, kSigmaLayerVar, NC_DOUBLE, 1, &dimLayer,
                             &sigmaLayerVar)) ||
        (status = nc_def_var(ncid, kSigmaIfaceVar, NC_DOUBLE, 1, &dimInterface,
                             &sigmaInterfaceVar)) ||
        (status = nc_def_var(ncid, kDepthVar, NC_DOUBLE, 1, &dimFace, &depthVar)))
    {
        return impl.fail(status, "defining vertical coordinate variables");
    }

    const std::string sigmaStandard = "ocean_sigma_coordinate";
    const std::string formulaTerms  = "sigma: " + std::string(kSigmaLayerVar) +
                                     " eta: " + impl.surfaceFieldId +
                                     " depth: " + kDepthVar;

    if ((status = nc_put_att_text(ncid, sigmaLayerVar, "standard_name",
                                  sigmaStandard.size(), sigmaStandard.c_str())) ||
        (status = nc_put_att_text(ncid, sigmaLayerVar, "formula_terms",
                                  formulaTerms.size(), formulaTerms.c_str())) ||
        (status = nc_put_att_text(ncid, sigmaLayerVar, "positive", positiveUpText.size(),
                                  positiveUpText.c_str())))
    {
        return impl.fail(status, "writing sigma coordinate attributes");
    }

    const std::string interfaceLongName = "sigma at layer interfaces";
    if ((status = nc_put_att_text(ncid, sigmaInterfaceVar, "long_name",
                                  interfaceLongName.size(), interfaceLongName.c_str())))
    {
        return impl.fail(status, "writing interface sigma attributes");
    }

    const std::string depthStandard = "sea_floor_depth_below_geoid";
    const std::string positiveDown  = "down";

    if ((status = nc_put_att_text(ncid, depthVar, "standard_name", depthStandard.size(),
                                  depthStandard.c_str())) ||
        (status = nc_put_att_text(ncid, depthVar, "units", metres.size(),
                                  metres.c_str())) ||
        (status = nc_put_att_text(ncid, depthVar, "positive", positiveDown.size(),
                                  positiveDown.c_str())))
    {
        return impl.fail(status, "writing depth attributes");
    }

    // ── Time ──────────────────────────────────────────────────────────────
    if ((status = nc_def_var(ncid, kTimeVar, NC_DOUBLE, 1, &dimTime, &impl.timeVarid)))
        return impl.fail(status, "defining time variable");

    const std::string timeUnits    = "days since 1858-11-17 00:00:00";
    const std::string timeStandard = "time";

    if ((status = nc_put_att_text(ncid, impl.timeVarid, "units", timeUnits.size(),
                                  timeUnits.c_str())) ||
        (status = nc_put_att_text(ncid, impl.timeVarid, "standard_name",
                                  timeStandard.size(), timeStandard.c_str())))
    {
        return impl.fail(status, "writing time attributes");
    }

    // ── Field variables ───────────────────────────────────────────────────
    const std::string meshName     = impl.mesh.meshName;
    const std::string faceLocation = "face";

    for (FieldSpec &field : impl.fields)
    {
        const int layeredDims[3] = {dimTime, dimFace, dimLayer};
        const int columnDims[2]  = {dimTime, dimFace};

        const int  rank = field.layered ? 3 : 2;
        const int *dims = field.layered ? layeredDims : columnDims;

        if ((status = nc_def_var(ncid, field.itemId.c_str(), NC_DOUBLE, rank, dims,
                                 &field.varid)))
        {
            return impl.fail(status, "defining variable '" + field.itemId + "'");
        }

        const size_t chunkFaces = static_cast<size_t>(std::min<int64_t>(faces, 4096));
        const size_t layeredChunks[3] = {1, chunkFaces, static_cast<size_t>(layers)};
        const size_t columnChunks[2]  = {1, chunkFaces};

        if ((status = nc_def_var_chunking(ncid, field.varid, NC_CHUNKED,
                                          field.layered ? layeredChunks : columnChunks)) ||
            (status = nc_def_var_deflate(ncid, field.varid, 0, 1, 4)))
        {
            return impl.fail(status, "configuring storage for '" + field.itemId + "'");
        }

        if ((status = nc_put_att_text(ncid, field.varid, "mesh", meshName.size(),
                                      meshName.c_str())) ||
            (status = nc_put_att_text(ncid, field.varid, "location",
                                      faceLocation.size(), faceLocation.c_str())))
        {
            return impl.fail(status, "tagging '" + field.itemId + "' with mesh metadata");
        }

        if (!field.units.empty() &&
            (status = nc_put_att_text(ncid, field.varid, "units", field.units.size(),
                                      field.units.c_str())))
        {
            return impl.fail(status, "writing units of '" + field.itemId + "'");
        }

        if (!field.standardName.empty() &&
            (status = nc_put_att_text(ncid, field.varid, "standard_name",
                                      field.standardName.size(),
                                      field.standardName.c_str())))
        {
            return impl.fail(status, "writing standard_name of '" + field.itemId + "'");
        }

        if (field.layered)
        {
            const std::string coordinates = kSigmaLayerVar;
            if ((status = nc_put_att_text(ncid, field.varid, "coordinates",
                                          coordinates.size(), coordinates.c_str())))
            {
                return impl.fail(status, "linking '" + field.itemId +
                                             "' to the sigma coordinate");
            }
        }
    }

    if ((status = nc_enddef(ncid)))
        return impl.fail(status, "ending define mode");

    // ── Static payloads ───────────────────────────────────────────────────
    if ((status = nc_put_var_double(ncid, nodeXVar, impl.mesh.nodeX.data())) ||
        (status = nc_put_var_double(ncid, nodeYVar, impl.mesh.nodeY.data())))
    {
        return impl.fail(status, "writing node coordinates");
    }

    if (hasNodeZ &&
        (status = nc_put_var_double(ncid, nodeZVar, impl.mesh.nodeZ.data())))
    {
        return impl.fail(status, "writing node elevations");
    }

    std::vector<int64_t> padded(static_cast<size_t>(faces * maxFaceNodes),
                                kConnectivityFill);

    for (int64_t f = 0; f < faces; ++f)
    {
        const int64_t begin = impl.mesh.faceNodeOffsets[static_cast<size_t>(f)];
        const int64_t end   = impl.mesh.faceNodeOffsets[static_cast<size_t>(f) + 1];

        for (int64_t slot = 0; slot < end - begin; ++slot)
        {
            padded[static_cast<size_t>(f * maxFaceNodes + slot)] =
                impl.mesh.faceNodes[static_cast<size_t>(begin + slot)];
        }
    }

    if ((status = nc_put_var_longlong(ncid, faceNodesVar,
                                      reinterpret_cast<const long long *>(padded.data()))))
    {
        return impl.fail(status, "writing face-node connectivity");
    }

    std::vector<double> cfInterfaces(static_cast<size_t>(layers) + 1, 0.0);
    for (size_t k = 0; k < cfInterfaces.size(); ++k)
        cfInterfaces[k] = -impl.vertical.sigmaInterfaces[k];

    std::vector<double> cfLayers(static_cast<size_t>(layers), 0.0);
    for (size_t k = 0; k < cfLayers.size(); ++k)
        cfLayers[k] = 0.5 * (cfInterfaces[k] + cfInterfaces[k + 1]);

    std::vector<double> depth(static_cast<size_t>(faces), 0.0);
    for (size_t f = 0; f < depth.size(); ++f)
        depth[f] = -impl.vertical.bedElevation[f];

    if ((status = nc_put_var_double(ncid, sigmaLayerVar, cfLayers.data())) ||
        (status = nc_put_var_double(ncid, sigmaInterfaceVar, cfInterfaces.data())) ||
        (status = nc_put_var_double(ncid, depthVar, depth.data())))
    {
        return impl.fail(status, "writing vertical coordinate payloads");
    }

    impl.error.clear();
    return NC_NOERR;
}

int LayeredUgridWriter::write(const SDKIO::Snapshot &snapshot)
{
    Impl &impl = *m_impl;

    if (impl.ncid < 0)
    {
        impl.error = "writer was not initialized";
        return NC_EBADID;
    }

    const int64_t faces  = impl.mesh.faceCount();
    const int64_t layers = impl.vertical.layerCount();
    const size_t  index  = static_cast<size_t>(impl.timeIndex);

    const double time   = snapshot.time();
    const size_t start  = index;
    const size_t count  = 1;
    int          status = nc_put_vara_double(impl.ncid, impl.timeVarid, &start, &count,
                                             &time);

    if (status != NC_NOERR)
        return impl.fail(status, "appending a time level");

    for (const FieldSpec &field : impl.fields)
    {
        const auto found = std::find_if(snapshot.fields().begin(), snapshot.fields().end(),
                                        [&](const SDKIO::FieldSlice &slice) {
                                            return slice.itemId == field.itemId;
                                        });

        if (found == snapshot.fields().end())
        {
            impl.error = "snapshot is missing registered field '" + field.itemId + "'";
            return NC_EINVAL;
        }

        const int64_t expected = field.layered ? faces * layers : faces;

        if (found->elementCount() != expected)
        {
            impl.error = "field '" + field.itemId + "' has " +
                         std::to_string(found->elementCount()) + " values; expected " +
                         std::to_string(expected);
            return NC_EINVAL;
        }

        if (found->kind != HydroCouple::DataKind::Float64)
        {
            impl.error = "field '" + field.itemId +
                         "' must be Float64 (staging conversion is not implemented yet)";
            return NC_EINVAL;
        }

        const size_t layeredStart[3] = {index, 0, 0};
        const size_t layeredCount[3] = {1, static_cast<size_t>(faces),
                                        static_cast<size_t>(layers)};
        const size_t columnStart[2]  = {index, 0};
        const size_t columnCount[2]  = {1, static_cast<size_t>(faces)};

        status = nc_put_vara_double(impl.ncid, field.varid,
                                    field.layered ? layeredStart : columnStart,
                                    field.layered ? layeredCount : columnCount,
                                    found->values<double>());

        if (status != NC_NOERR)
            return impl.fail(status, "writing '" + field.itemId + "'");
    }

    ++impl.timeIndex;
    impl.error.clear();
    return NC_NOERR;
}

int LayeredUgridWriter::finalize()
{
    Impl &impl = *m_impl;

    if (impl.ncid < 0)
        return NC_NOERR;

    const int status = nc_close(impl.ncid);
    impl.ncid        = -1;

    if (status != NC_NOERR)
    {
        impl.error = std::string("closing the file: ") + nc_strerror(status);
        return status;
    }

    return NC_NOERR;
}

#else // FVQual_HAVE_NETCDF

namespace {
constexpr int kNoNetCDF = -1;

const char *kNoNetCDFMessage =
    "this FVQual build has NetCDF support disabled "
    "(configure with -DFVQual_USE_NETCDF=ON)";
} // namespace

int LayeredUgridWriter::initialize()
{
    m_impl->error = kNoNetCDFMessage;
    return kNoNetCDF;
}

int LayeredUgridWriter::write(const SDKIO::Snapshot &)
{
    m_impl->error = kNoNetCDFMessage;
    return kNoNetCDF;
}

int LayeredUgridWriter::finalize()
{
    return 0;
}

#endif // FVQual_HAVE_NETCDF

} // namespace FVQual::IO
