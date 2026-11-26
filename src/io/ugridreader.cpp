// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
#include "fvqual/io/ugridreader.h"
#include "hydrocouplesdk/io/ugridreader.h"
namespace FVQual::IO {
bool UgridReader::available() noexcept {
#ifdef FVQual_HAVE_NETCDF
    return HydroCouple::SDK::IO::ugridReadSupported();
#else
    return false;
#endif
}
bool UgridReader::read(const std::filesystem::path &path, const std::string &meshName,
                      HydroCouple::SDK::IO::MeshDefinition &mesh, std::string &message) {
#ifdef FVQual_HAVE_NETCDF
    HydroCouple::SDK::IO::MeshDefinition candidate;
    if (!HydroCouple::SDK::IO::readUGRIDMesh(path, meshName, candidate, message)) return false;
    mesh = std::move(candidate);
    return true;
#else
    (void)meshName; (void)mesh;
    message = "cannot read '" + path.string() +
        "': this FVQual build has NetCDF support disabled (configure with -DFVQual_USE_NETCDF=ON)";
    return false;
#endif
}
} // namespace FVQual::IO
