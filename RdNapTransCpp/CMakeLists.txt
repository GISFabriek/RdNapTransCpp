# CMakeList.txt : CMake project for RdNapTransCpp, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Add source to this project's executable.
add_executable ("RdNapTrans" "Cartesian.h" "Cartesian.cpp" "Geographic.h"  "Geographic.cpp""Constants.h" "Helpers.h" "Helpers.cpp" "Transformer.h" "Transformer.cpp" "GrdFile.h" "GrdFile.cpp" "Runner.cpp")

# TODO: Add tests and install targets if needed.
