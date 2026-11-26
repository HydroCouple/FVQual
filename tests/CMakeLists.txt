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
    test_meshqualification.cpp
    test_completion.cpp
    test_sigmaspec.cpp
    test_layeredmesh.cpp
    test_hypsography.cpp
    test_w2bathymetry.cpp
    test_barotropic.cpp
    test_equationofstate.cpp
    test_internalmode.cpp
    test_momentumgeometry.cpp
    test_scalartrace.cpp
    test_reactivetrace.cpp
    test_thermalcolumn.cpp
    test_transport.cpp
    test_sources.cpp
    test_structures.cpp
    test_outletgroups.cpp
    test_transfers.cpp
    test_openboundaries.cpp
    test_boundarytransport.cpp
    test_heatbudget.cpp
    test_forcingfield.cpp
    test_turbulence.cpp
    test_verticalcoordinate.cpp
    test_fvqualcomponent.cpp
    test_fvqualdifferentiation.cpp
    test_reservoircontract.cpp
    test_liveinflows.cpp
    test_operations.cpp
    test_grossturnover.cpp
    test_fvqualplugin.cpp
    test_fvqualmeshargument.cpp
    test_fvqualforcingarguments.cpp
    test_fvqualruncli.cpp
    test_conveyance.cpp
    test_macdonald.cpp
    test_kinetics_numerics.cpp
    test_reactionnetwork.cpp
    test_lightfield.cpp
    test_settling.cpp
    test_templates.cpp
    test_reactionfile.cpp
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
    PRIVATE FVQual_TEST_OUTPUT_DIR="${FVQual_TEST_OUTPUT_DIR}"
    PRIVATE FVQUAL_PLUGIN_PATH="$<TARGET_FILE:fvqual_component>"
    PRIVATE FVQUAL_CLI_PATH="$<TARGET_FILE:fvqual>"
    PRIVATE FVQUAL_FIXTURES_DIR="${CMAKE_SOURCE_DIR}/tests/fixtures"
    PRIVATE FVQUAL_COMPOSITIONS_DIR="${CMAKE_SOURCE_DIR}/tests/compositions")

if(FVQual_USE_NETCDF)
    target_compile_definitions(fvqual_tests PRIVATE FVQual_HAVE_NETCDF)
endif()

fvqual_colocate_runtime_dlls(fvqual_tests)

add_subdirectory(capi)

# ctypes exercises the shipped shared ABI, including handle ownership and
# recorded/live transfer accounting, without a Python extension build.
if(BUILD_SHARED_LIBS)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        add_test(NAME fvqual_python_workflow
            COMMAND ${CMAKE_COMMAND} -E env
                "FVQUAL_LIBRARY=$<TARGET_FILE:fvqual_core>"
                "PYTHONPATH=${CMAKE_SOURCE_DIR}/python"
                ${Python3_EXECUTABLE} -m unittest discover
                -s "${CMAKE_SOURCE_DIR}/python/tests" -v)
        set_tests_properties(fvqual_python_workflow PROPERTIES TIMEOUT 120)
        add_test(NAME fvqual_degray_observation_operator
            COMMAND ${Python3_EXECUTABLE} -m unittest discover
                -s "${CMAKE_SOURCE_DIR}/verification/degray" -p "test_*.py" -v)
        set_tests_properties(fvqual_degray_observation_operator PROPERTIES TIMEOUT 30)
    endif()
endif()

include(GoogleTest)
add_dependencies(fvqual_tests fvqual_component fvqual)

gtest_discover_tests(fvqual_tests
    DISCOVERY_MODE PRE_TEST
    PROPERTIES TIMEOUT 120)
