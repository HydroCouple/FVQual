# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Caleb Buahin

add_executable(fvqual-bench fvqual_bench.cpp)

target_link_libraries(fvqual-bench PRIVATE FVQual::fvqual_core)
target_compile_features(fvqual-bench PRIVATE cxx_std_20)
