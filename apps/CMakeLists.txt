# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Caleb Buahin

add_executable(fvqual fvqual/main.cpp)

target_link_libraries(fvqual PRIVATE FVQual::fvqual_core)
target_compile_features(fvqual PRIVATE cxx_std_20)

if(FVQual_USE_NETCDF)
    target_compile_definitions(fvqual PRIVATE FVQual_HAVE_NETCDF)
endif()

install(TARGETS fvqual RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
