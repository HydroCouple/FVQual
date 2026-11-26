# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Caleb Buahin
#
# The PY0 gate: the C API exercised from C, by a C compiler.

add_executable(fvqual_capi_tests test_capi.c c89_header_check.c)

target_link_libraries(fvqual_capi_tests PRIVATE FVQual::fvqual_core)
fvqual_colocate_runtime_dlls(fvqual_capi_tests)

if(NOT MSVC)
    target_link_libraries(fvqual_capi_tests PRIVATE m)
endif()

set_target_properties(fvqual_capi_tests PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON)

# The header check is the point of the exercise, so hold it to C89 strictly.
# MSVC has no C89 conformance mode worth the name; the GCC/Clang jobs carry it.
if(NOT MSVC)
    set_source_files_properties(c89_header_check.c
        PROPERTIES COMPILE_OPTIONS "-std=c89;-pedantic-errors;-Wall")
endif()

add_test(NAME capi.smoke COMMAND fvqual_capi_tests)
set_tests_properties(capi.smoke PROPERTIES TIMEOUT 120)
