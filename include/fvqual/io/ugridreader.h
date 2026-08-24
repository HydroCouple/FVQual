// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   ugridreader.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  UgridReader — reads a UGRID-1.0 two-dimensional mesh from a
 *         NetCDF-4 file (F1.2).
 *
 * The SDK ships UGRID writers but no reader; FVQual needs one because a mesh
 * produced by any UGRID-aware tool (the SDK's own meshers, Delft3D-FM's grid
 * editor, QGIS plug-ins, xugrid) is the intended way to define a domain.
 *
 * Handles the mixed tri-quad case UGRID expresses as a padded connectivity
 * matrix: `face_node_connectivity` is shaped { face, maxNodesPerFace } with
 * unused slots set to the variable's `_FillValue`, so a triangle in a
 * quad-containing mesh carries one fill entry. Both `start_index` conventions
 * (0- and 1-based) are honoured.
 *
 * The NetCDF C API is fully encapsulated; this header is installable without
 * NetCDF present. When FVQual is built without NetCDF, read() fails with an
 * explanatory message rather than being absent.
 *
 * \threadsafety Stateless static entry point.
 */

#ifndef FVQual_IO_UGRIDREADER_H
#define FVQual_IO_UGRIDREADER_H

#include "fvqual/fvqual_export.h"

#include "hydrocouplesdk/io/meshdefinition.h"

#include <filesystem>
#include <string>

namespace FVQual::IO {

/*!
 * \brief UGRID-1.0 two-dimensional mesh reader.
 */
class FVQual_EXPORT UgridReader
{
public:
    /*!
     * \brief Reads the first (or a named) 2-D mesh topology from a file.
     *
     * Node z values are taken from a variable named in the topology's
     * optional `node_z` / `node_elevation` attribute, or from a variable
     * carrying `standard_name = "sea_floor_depth_below_geoid"` (negated, so
     * the result is an elevation, positive up). Absent any of these, node
     * elevations are left empty and the caller must supply bathymetry.
     *
     * \param[in]  path     NetCDF-4 file to read.
     * \param[in]  meshName Mesh topology variable to read; empty selects the
     *                      first 2-D mesh found.
     * \param[out] mesh     The mesh read.
     * \param[out] message  Failure description when returning false.
     */
    [[nodiscard]] static bool read(const std::filesystem::path &path,
                                   const std::string &meshName,
                                   HydroCouple::SDK::IO::MeshDefinition &mesh,
                                   std::string &message);

    /*!
     * \brief Convenience overload reading the first 2-D mesh in the file.
     */
    [[nodiscard]] static bool read(const std::filesystem::path &path,
                                   HydroCouple::SDK::IO::MeshDefinition &mesh,
                                   std::string &message)
    {
        return read(path, std::string{}, mesh, message);
    }

    /*!
     * \brief Whether this build can read NetCDF at all.
     */
    [[nodiscard]] static bool available() noexcept;
};

} // namespace FVQual::IO

#endif // FVQual_IO_UGRIDREADER_H
