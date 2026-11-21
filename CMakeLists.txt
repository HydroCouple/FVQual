# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Caleb Buahin
#
# CMakeLists.txt — FVQual top-level build configuration
#
# FVQual: a fully 3D, unstructured, finite-volume hydrodynamic and
# water-quality model — the HydroCouple-native successor to CE-QUAL-W2.
#
# Usage:
#   cmake --preset linux            (or macos / windows)
#   cmake --build build --parallel
#   ctest --test-dir build --output-on-failure

cmake_minimum_required(VERSION 3.21...4.2)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR
        "In-source builds are not allowed. Configure into a separate build directory.")
endif()

project(
    FVQual
    VERSION 0.1.0
    LANGUAGES C CXX
    HOMEPAGE_URL "https://hydrocouple.org/projects/fvqual"
    DESCRIPTION  "3D unstructured finite-volume hydrodynamic and water-quality model"
)

# Pre-release label appended wherever the full version string is used.
set(FVQual_PRERELEASE "alpha.1" CACHE STRING
    "Pre-release tag (e.g. alpha.1, beta.1, rc.1). Empty for a final release.")

if(FVQual_PRERELEASE)
    set(FVQual_FULL_VERSION "${PROJECT_VERSION}-${FVQual_PRERELEASE}")
else()
    set(FVQual_FULL_VERSION "${PROJECT_VERSION}")
endif()

# ── Language standard ─────────────────────────────────────────────────────
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(GNUInstallDirs)
include(GenerateExportHeader)

# ── Options ───────────────────────────────────────────────────────────────
option(FVQual_BUILD_TESTS  "Build the unit/verification test suite" ON)
option(FVQual_BUILD_APPS   "Build the fvqual command-line runner"   ON)
option(FVQual_BUILD_VERIFICATION
       "Build fvqual-verify, which regenerates the animated verification page" ON)
option(FVQual_USE_NETCDF   "Enable NetCDF-4 UGRID input/output"     OFF)
option(FVQual_USE_OPENMP   "Enable OpenMP threading"                OFF)
option(FVQual_BUILD_DOCS   "Build Doxygen documentation"             OFF)

# ── HydroCoupleSDK (required) ─────────────────────────────────────────────
#
# Preferred: an installed SDK found via CMAKE_PREFIX_PATH.
# Fallback:  point FVQual_SDK_SOURCE_DIR at a source checkout and it is
#            added to the build directly. The SDK is not published to a
#            package registry yet, so the fallback is the usual path today.
set(FVQual_SDK_SOURCE_DIR "" CACHE PATH
    "HydroCoupleSDK source checkout to build in-tree when no installed package is found")

find_package(HydroCoupleSDK 2.0 CONFIG QUIET)

if(NOT HydroCoupleSDK_FOUND)
    if(FVQual_SDK_SOURCE_DIR AND EXISTS "${FVQual_SDK_SOURCE_DIR}/CMakeLists.txt")
        message(STATUS "FVQual: building HydroCoupleSDK in-tree from ${FVQual_SDK_SOURCE_DIR}")
        add_subdirectory("${FVQual_SDK_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/_hydrocouplesdk")
    else()
        message(FATAL_ERROR
            "HydroCoupleSDK 2.0 was not found.\n"
            "  Either install it and pass -DCMAKE_PREFIX_PATH=<sdk-install-prefix>,\n"
            "  or pass -DFVQual_SDK_SOURCE_DIR=<path-to-HydroCoupleSDK> to build it in-tree.")
    endif()
endif()

# Both discovery paths must yield the same target name.
if(NOT TARGET HydroCoupleSDK::HydroCoupleSDK AND TARGET HydroCoupleSDK)
    add_library(HydroCoupleSDK::HydroCoupleSDK ALIAS HydroCoupleSDK)
endif()

if(NOT TARGET HydroCoupleSDK::HydroCoupleTools AND TARGET HydroCoupleTools)
    add_library(HydroCoupleSDK::HydroCoupleTools ALIAS HydroCoupleTools)
endif()

if(NOT TARGET HydroCoupleSDK::HydroCoupleTools)
    message(FATAL_ERROR
        "FVQual requires the HydroCoupleSDK tools library (sigma grids, meshing).\n"
        "  Rebuild the SDK with -DBUILD_TOOLS=ON.")
endif()

# ── Optional dependencies ─────────────────────────────────────────────────
if(FVQual_USE_NETCDF)
    find_package(netCDF CONFIG QUIET)
    if(NOT netCDF_FOUND)
        find_package(PkgConfig QUIET)
        if(PkgConfig_FOUND)
            pkg_check_modules(NETCDF IMPORTED_TARGET netcdf)
        endif()
        if(NOT NETCDF_FOUND)
            message(FATAL_ERROR "FVQual_USE_NETCDF=ON but NetCDF was not found.")
        endif()
    endif()
endif()

if(FVQual_USE_OPENMP)
    find_package(OpenMP REQUIRED)
endif()

# ── Version header ────────────────────────────────────────────────────────
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/generated/fvqual/version.h"
    @ONLY)

# ── Targets ───────────────────────────────────────────────────────────────
add_subdirectory(src)

if(FVQual_BUILD_APPS)
    add_subdirectory(apps)
endif()

if(FVQual_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(FVQual_BUILD_VERIFICATION)
    add_subdirectory(verification)
add_subdirectory(benchmarks)
endif()

if(FVQual_BUILD_DOCS)
    find_package(Doxygen REQUIRED)
    add_subdirectory(docs)
endif()

# ── Install / export ──────────────────────────────────────────────────────
include(CMakePackageConfigHelpers)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

install(EXPORT FVQualTargets
        FILE        FVQualTargets.cmake
        NAMESPACE   FVQual::
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/FVQual")

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FVQualConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/FVQualConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/FVQual")

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/FVQualConfigVersion.cmake"
    VERSION       "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/FVQualConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/FVQualConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/FVQual")

# ── Summary ───────────────────────────────────────────────────────────────
message(STATUS "")
message(STATUS "FVQual ${FVQual_FULL_VERSION}")
message(STATUS "  C++ standard ......... ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build type ........... ${CMAKE_BUILD_TYPE}")
message(STATUS "  Tests ................ ${FVQual_BUILD_TESTS}")
message(STATUS "  Apps ................. ${FVQual_BUILD_APPS}")
message(STATUS "  Verification tool .... ${FVQual_BUILD_VERIFICATION}")
message(STATUS "  NetCDF UGRID IO ...... ${FVQual_USE_NETCDF}")
message(STATUS "  OpenMP ............... ${FVQual_USE_OPENMP}")
message(STATUS "")
