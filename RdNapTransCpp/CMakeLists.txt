# CMakeList.txt : CMake project for RdNapTransCpp, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)
include(CMakeRC.cmake)
cmrc_add_resource_library(rdnaptransrc nlgeo04.b64 x2c.b64 y2c.b64)

# Add source to this project's executable.
add_executable ("RdNapTrans" "Cartesian.h" "Cartesian.cpp" "Geographic.h"  "Geographic.cpp""Constants.h" "Helpers.h" "Helpers.cpp" "Transformer.h" "Transformer.cpp" "GrdFile.h" "GrdFile.cpp" "Runner.cpp")
target_link_libraries(RdNapTrans PRIVATE rdnaptransrc)
# TODO: Add tests and install targets if needed.
