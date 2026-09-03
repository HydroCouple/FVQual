# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Caleb Buahin
#
# libfvqual-core — mesh, state, and IO for the FVQual model.

add_library(fvqual_core
    mesh/horizontalmesh.cpp
    mesh/layeredmesh.cpp
    mesh/sigmaspec.cpp
    hydro/barotropicsolver.cpp
    hydro/equationofstate.cpp
    hydro/pressureclosure.cpp
    hydro/internalmode.cpp
    hydro/heatbudget.cpp
    hydro/turbulence.cpp
    transport/euleriantransport.cpp
    forcing/forcingfield.cpp
    forcing/sources.cpp
    forcing/structures.cpp
    forcing/outletgroups.cpp
    kinetics/expression.cpp
    kinetics/odesolver.cpp
    kinetics/reactionnetwork.cpp
    kinetics/lightfield.cpp
    kinetics/templates.cpp
    kinetics/reactionfile.cpp
    transport/settling.cpp
    forcing/transfers.cpp
    state/speciesregistry.cpp
    state/statestore.cpp
    io/w2bathymetry.cpp
    io/w2timeseries.cpp
    io/ugridreader.cpp
    io/layeredugridwriter.cpp
    component/sigmaverticalcoordinate.cpp
    component/fvqualcomponent.cpp
    component/fvqualcomponentinfo.cpp
    capi/capi_engine.cpp
    capi/capi_mesh.cpp
    capi/capi_sim.cpp
)

add_library(FVQual::fvqual_core ALIAS fvqual_core)

set_target_properties(fvqual_core PROPERTIES
    OUTPUT_NAME   fvqual-core
    VERSION       "${PROJECT_VERSION}"
    SOVERSION     "${PROJECT_VERSION_MAJOR}"
    EXPORT_NAME   fvqual_core
)

# The macro names are given explicitly because generate_export_header()
# upper-cases BASE_NAME internally, which would emit FVQUAL_EXPORT and no
# longer match the mixed-case spelling used throughout the headers.
generate_export_header(fvqual_core
    BASE_NAME             FVQual
    EXPORT_MACRO_NAME     FVQual_EXPORT
    NO_EXPORT_MACRO_NAME  FVQual_NO_EXPORT
    DEPRECATED_MACRO_NAME FVQual_DEPRECATED
    EXPORT_FILE_NAME      "${CMAKE_BINARY_DIR}/generated/fvqual/fvqual_export.h")

target_include_directories(fvqual_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/generated>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(fvqual_core
    PUBLIC
        HydroCoupleSDK::HydroCoupleSDK
        HydroCoupleSDK::HydroCoupleTools
)

target_compile_features(fvqual_core PUBLIC cxx_std_20)

# The C API is exported from this library; the build-side define flips the
# Windows dllexport/dllimport switch (fvqual_engine.h).
target_compile_definitions(fvqual_core PRIVATE FVQual_CAPI_BUILD)

if(FVQual_USE_NETCDF)
    target_compile_definitions(fvqual_core PRIVATE FVQual_HAVE_NETCDF)

    if(TARGET netCDF::netcdf)
        target_link_libraries(fvqual_core PRIVATE netCDF::netcdf)
    elseif(TARGET PkgConfig::NETCDF)
        target_link_libraries(fvqual_core PRIVATE PkgConfig::NETCDF)
    endif()
endif()

if(FVQual_USE_OPENMP)
    target_link_libraries(fvqual_core PUBLIC OpenMP::OpenMP_CXX)
endif()

if(MSVC)
    target_compile_options(fvqual_core PRIVATE /W4 /permissive-)
else()
    target_compile_options(fvqual_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

install(TARGETS fvqual_core
        EXPORT  FVQualTargets
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")

# ── The loadable component plugin ─────────────────────────────────────────
#
# FVQual published through the SDK's stamped ABI convention
# (hydrocouple_component_abi_v1 / hydrocouple_component_info_v1), so
# HydroCoupleComposer's registry, the headless runner, and the Python
# loader all pick it up identically. The core stays a plain library an
# embedder links; the entry points live here, where a loader looks.
add_library(fvqual_component MODULE component/fvqualplugin.cpp)
add_library(FVQual::fvqual_component ALIAS fvqual_component)
target_link_libraries(fvqual_component PRIVATE FVQual::fvqual_core)
set_target_properties(fvqual_component PROPERTIES
    OUTPUT_NAME fvqual-component
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/components")
install(TARGETS fvqual_component
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/hydrocouple/components)
