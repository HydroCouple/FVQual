# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Caleb Buahin

find_package(GTest CONFIG QUIET)

if(NOT GTest_FOUND)
    find_package(GTest QUIET)
endif()

if(NOT GTest_FOUND AND NOT TARGET GTest::gtest_main)
    message(WARNING
        "GoogleTest was not found — the FVQual test suite will not be built. "
        "Install it (vcpkg feature 'tests') or set -DFVQual_BUILD_TESTS=OFF.")
    return()
endif()

add_executable(fvqual_tests
    test_horizontalmesh.cpp
    test_sigmaspec.cpp
    test_layeredmesh.cpp
    test_barotropic.cpp
    test_equationofstate.cpp
    test_internalmode.cpp
    test_transport.cpp
    test_sources.cpp
    test_structures.cpp
    test_outletgroups.cpp
    test_transfers.cpp
    test_openboundaries.cpp
    test_macdonald.cpp
    test_kinetics_numerics.cpp
    test_reactionnetwork.cpp
    test_statestore.cpp
    test_ugrid_roundtrip.cpp
)

target_link_libraries(fvqual_tests
    PRIVATE
        FVQual::fvqual_core
        GTest::gtest
        GTest::gtest_main
)

target_compile_features(fvqual_tests PRIVATE cxx_std_20)

# Tests write reviewable artifacts here rather than to a temporary directory,
# so a failed run can be inspected afterwards (CONTRIBUTING: reviewable IO).
set(FVQual_TEST_OUTPUT_DIR "${CMAKE_SOURCE_DIR}/tests/output")
file(MAKE_DIRECTORY "${FVQual_TEST_OUTPUT_DIR}")

target_compile_definitions(fvqual_tests
    PRIVATE FVQual_TEST_OUTPUT_DIR="${FVQual_TEST_OUTPUT_DIR}")

if(FVQual_USE_NETCDF)
    target_compile_definitions(fvqual_tests PRIVATE FVQual_HAVE_NETCDF)
endif()

add_subdirectory(capi)

include(GoogleTest)
gtest_discover_tests(fvqual_tests
    DISCOVERY_MODE PRE_TEST
    PROPERTIES TIMEOUT 120)
