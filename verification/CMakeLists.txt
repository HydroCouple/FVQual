# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Caleb Buahin
#
# fvqual-verify — runs the analytical benchmarks and emits animation frames
# for verification/index.html.

add_executable(fvqual-verify scenarios.cpp)

target_link_libraries(fvqual-verify PRIVATE FVQual::fvqual_core)
target_compile_features(fvqual-verify PRIVATE cxx_std_20)

# Not installed: this is a development and documentation tool, not part of the
# shipped model.
set_target_properties(fvqual-verify PROPERTIES EXCLUDE_FROM_ALL FALSE)
