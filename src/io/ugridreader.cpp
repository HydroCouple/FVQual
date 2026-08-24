// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/io/ugridreader.h"

#ifdef FVQual_HAVE_NETCDF

#include <netcdf.h>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

namespace FVQual::IO {

namespace {

namespace SDKIO = HydroCouple::SDK::IO;

/*!
 * \brief RAII handle for an open NetCDF file.
 */
class NcFile
{
public:
    ~NcFile()
    {
        if (m_id >= 0)
            nc_close(m_id);
    }

    [[nodiscard]] int open(const std::string &path)
    {
        return nc_open(path.c_str(), NC_NOWRITE, &m_id);
    }

    [[nodiscard]] int id() const noexcept { return m_id; }

private:
    int m_id = -1;
};

/*!
 * \brief Reads a text attribute; returns false when it is absent.
 */
bool textAttribute(int ncid, int varid, const char *name, std::string &value)
{
    size_t length = 0;
    if (nc_inq_attlen(ncid, varid, name, &length) != NC_NOERR)
        return false;

    value.assign(length, '\0');
    if (nc_get_att_text(ncid, varid, name, value.data()) != NC_NOERR)
        return false;

    // NetCDF text attributes are not necessarily NUL-terminated, but may be.
    while (!value.empty() && value.back() == '\0')
        value.pop_back();

    return true;
}

std::vector<std::string> splitWords(const std::string &text)
{
    std::vector<std::string> words;
    std::istringstream       stream(text);
    std::string              word;

    while (stream >> word)
        words.push_back(word);

    return words;
}

/*!
 * \brief Locates the mesh-topology dummy variable.
 */
bool findTopology(int ncid, const std::string &requestedName, int &varid,
                  std::string &foundName, std::string &message)
{
    int variableCount = 0;
    if (nc_inq_nvars(ncid, &variableCount) != NC_NOERR)
    {
        message = "could not enumerate variables";
        return false;
    }

    for (int v = 0; v < variableCount; ++v)
    {
        char name[NC_MAX_NAME + 1] = {};
        if (nc_inq_varname(ncid, v, name) != NC_NOERR)
            continue;

        std::string role;
        if (!textAttribute(ncid, v, "cf_role", role) || role != "mesh_topology")
            continue;

        int dimension = 0;
        if (nc_get_att_int(ncid, v, "topology_dimension", &dimension) != NC_NOERR ||
            dimension != 2)
        {
            continue;
        }

        if (!requestedName.empty() && requestedName != name)
            continue;

        varid     = v;
        foundName = name;
        return true;
    }

    message = requestedName.empty()
                  ? "no 2-D UGRID mesh topology variable found "
                    "(needs cf_role=\"mesh_topology\" and topology_dimension=2)"
                  : "mesh topology variable '" + requestedName + "' not found";
    return false;
}

/*!
 * \brief Reads a 1-D double variable in full.
 */
bool readDoubles(int ncid, const std::string &variable, std::vector<double> &values,
                 std::string &message)
{
    int varid = -1;
    if (nc_inq_varid(ncid, variable.c_str(), &varid) != NC_NOERR)
    {
        message = "variable '" + variable + "' not found";
        return false;
    }

    int dimensionCount = 0;
    if (nc_inq_varndims(ncid, varid, &dimensionCount) != NC_NOERR || dimensionCount != 1)
    {
        message = "variable '" + variable + "' is not one-dimensional";
        return false;
    }

    int dimid = -1;
    if (nc_inq_vardimid(ncid, varid, &dimid) != NC_NOERR)
    {
        message = "could not inspect dimensions of '" + variable + "'";
        return false;
    }

    size_t length = 0;
    if (nc_inq_dimlen(ncid, dimid, &length) != NC_NOERR)
    {
        message = "could not read length of '" + variable + "'";
        return false;
    }

    values.assign(length, 0.0);
    if (nc_get_var_double(ncid, varid, values.data()) != NC_NOERR)
    {
        message = "could not read values of '" + variable + "'";
        return false;
    }

    return true;
}

} // namespace

bool UgridReader::available() noexcept
{
    return true;
}

bool UgridReader::read(const std::filesystem::path &path, const std::string &meshName,
                       SDKIO::MeshDefinition &mesh, std::string &message)
{
    NcFile file;

    if (file.open(path.string()) != NC_NOERR)
    {
        message = "could not open '" + path.string() + "' as a NetCDF file";
        return false;
    }

    const int ncid = file.id();

    int         topologyVar = -1;
    std::string topologyName;
    if (!findTopology(ncid, meshName, topologyVar, topologyName, message))
        return false;

    mesh = SDKIO::MeshDefinition{};
    mesh.meshName = topologyName;

    // ── Node coordinates ──────────────────────────────────────────────────
    std::string nodeCoordinates;
    if (!textAttribute(ncid, topologyVar, "node_coordinates", nodeCoordinates))
    {
        message = "mesh '" + topologyName + "' has no node_coordinates attribute";
        return false;
    }

    const std::vector<std::string> coordinateNames = splitWords(nodeCoordinates);
    if (coordinateNames.size() < 2)
    {
        message = "mesh '" + topologyName +
                  "' node_coordinates must name at least an x and a y variable";
        return false;
    }

    if (!readDoubles(ncid, coordinateNames[0], mesh.nodeX, message) ||
        !readDoubles(ncid, coordinateNames[1], mesh.nodeY, message))
    {
        return false;
    }

    if (mesh.nodeX.size() != mesh.nodeY.size())
    {
        message = "node x and y coordinate variables have different lengths";
        return false;
    }

    // ── Face-node connectivity ────────────────────────────────────────────
    std::string connectivityName;
    if (!textAttribute(ncid, topologyVar, "face_node_connectivity", connectivityName))
    {
        message = "mesh '" + topologyName + "' has no face_node_connectivity attribute";
        return false;
    }

    int connectivityVar = -1;
    if (nc_inq_varid(ncid, connectivityName.c_str(), &connectivityVar) != NC_NOERR)
    {
        message = "connectivity variable '" + connectivityName + "' not found";
        return false;
    }

    int connectivityDims = 0;
    if (nc_inq_varndims(ncid, connectivityVar, &connectivityDims) != NC_NOERR ||
        connectivityDims != 2)
    {
        message = "'" + connectivityName + "' must be two-dimensional { face, node }";
        return false;
    }

    int dimIds[2] = {-1, -1};
    if (nc_inq_vardimid(ncid, connectivityVar, dimIds) != NC_NOERR)
    {
        message = "could not inspect dimensions of '" + connectivityName + "'";
        return false;
    }

    size_t faceCount    = 0;
    size_t maxFaceNodes = 0;
    if (nc_inq_dimlen(ncid, dimIds[0], &faceCount) != NC_NOERR ||
        nc_inq_dimlen(ncid, dimIds[1], &maxFaceNodes) != NC_NOERR)
    {
        message = "could not read the shape of '" + connectivityName + "'";
        return false;
    }

    std::vector<long long> raw(faceCount * maxFaceNodes, 0);
    if (nc_get_var_longlong(ncid, connectivityVar, raw.data()) != NC_NOERR)
    {
        message = "could not read values of '" + connectivityName + "'";
        return false;
    }

    long long startIndex = 0;
    nc_get_att_longlong(ncid, connectivityVar, "start_index", &startIndex);

    long long fillValue    = 0;
    const bool hasFillValue = nc_get_att_longlong(ncid, connectivityVar, "_FillValue",
                                                  &fillValue) == NC_NOERR;

    const auto nodeCount = static_cast<long long>(mesh.nodeX.size());

    mesh.faceNodeOffsets.assign(faceCount + 1, 0);
    mesh.faceNodes.clear();
    mesh.faceNodes.reserve(faceCount * maxFaceNodes);

    for (size_t f = 0; f < faceCount; ++f)
    {
        for (size_t slot = 0; slot < maxFaceNodes; ++slot)
        {
            const long long value = raw[f * maxFaceNodes + slot];

            if (hasFillValue && value == fillValue)
                break;

            const long long index = value - startIndex;

            if (index < 0 || index >= nodeCount)
            {
                // Absent an explicit _FillValue, an out-of-range entry is the
                // padding of a face with fewer nodes than the widest one.
                if (!hasFillValue)
                    break;

                message = "face " + std::to_string(f) + " references node " +
                          std::to_string(index) + ", which is out of range";
                return false;
            }

            mesh.faceNodes.push_back(static_cast<int64_t>(index));
        }

        mesh.faceNodeOffsets[f + 1] = static_cast<int64_t>(mesh.faceNodes.size());
    }

    // ── Optional node elevations ──────────────────────────────────────────
    std::string elevationName;

    if (!textAttribute(ncid, topologyVar, "node_z", elevationName))
        textAttribute(ncid, topologyVar, "node_elevation", elevationName);

    if (!elevationName.empty())
    {
        std::string ignored;
        if (readDoubles(ncid, elevationName, mesh.nodeZ, ignored) &&
            mesh.nodeZ.size() != mesh.nodeX.size())
        {
            message = "node elevation variable '" + elevationName +
                      "' has a different length than the node coordinates";
            return false;
        }
    }
    else if (coordinateNames.size() >= 3)
    {
        std::string ignored;
        if (readDoubles(ncid, coordinateNames[2], mesh.nodeZ, ignored) &&
            mesh.nodeZ.size() != mesh.nodeX.size())
        {
            mesh.nodeZ.clear();
        }
    }

    if (!mesh.validate(message))
        return false;

    message.clear();
    return true;
}

} // namespace FVQual::IO

#else // FVQual_HAVE_NETCDF

namespace FVQual::IO {

bool UgridReader::available() noexcept
{
    return false;
}

bool UgridReader::read(const std::filesystem::path &path, const std::string &,
                       HydroCouple::SDK::IO::MeshDefinition &, std::string &message)
{
    message = "cannot read '" + path.string() +
              "': this FVQual build has NetCDF support disabled "
              "(configure with -DFVQual_USE_NETCDF=ON)";
    return false;
}

} // namespace FVQual::IO

#endif // FVQual_HAVE_NETCDF
