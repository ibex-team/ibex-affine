cmake_minimum_required (VERSION 3.5.1)

project (IBEX_AFFINE VERSION 2.8.9 LANGUAGES CXX)
set (IBEX_AFFINE_DESCRIPTION "Affine arithmetic for the Ibex library")

find_package (IBEX 2.8.9 REQUIRED)
message (STATUS "Found Ibex version ${IBEX_VERSION}")
ibex_init_common () # Ibex should have installed this function

################################################################################
# Compile sources
################################################################################
add_subdirectory (src)

################################################################################
# Tests
################################################################################
# The following command assumes that the tests are in a tests/ subdirectory. So
# the tests/ subdirectory must exist and must contain a CMakeLists.txt file.
# If this is not the case, comment the next line.
add_make_target_for_ctest (check)
