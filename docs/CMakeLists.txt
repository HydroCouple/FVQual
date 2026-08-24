# SPDX-License-Identifier: MIT
# docs/CMakeLists.txt — Doxygen documentation target
#
# Mirrors the HydroCoupleSDK pipeline so the two sets of documentation are
# built and themed the same way.

set(DOXYGEN_IN  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

add_custom_target(docs ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating Doxygen documentation for FVQual"
    VERBATIM
)
