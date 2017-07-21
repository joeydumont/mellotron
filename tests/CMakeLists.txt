# --------------------------------------------------------------------------- #
# Author:       Denis Gagnon                  <denis.gagnon@emt.inrs.ca>      #
# Date created: 2017-07-21                                                    #
# Description:  CMake compilation instructions for MELLOTRON tests            #
# ----------------------------------------------------------------------------#

# ----------------------------------------------------------------- #
# --               Configuration and Dependencies                -- #
# ----------------------------------------------------------------- #
# -- CMake version and installation directory.
# CMake version
cmake_minimum_required(VERSION 3.1)


# -- Dependency (Google Test)
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GTEST_LIBRARIES})

# -- Output tests in directory
add_executable("UnitTest-Particle.o" "UnitTest-Particle.cpp")
target_link_libraries("UnitTest-Particle.o" ${LIBS})
add_test(NAME "UnitTest-Particle" COMMAND "UnitTest-Particle.o")

add_executable("UnitTest-eDipoles.o" "UnitTest-eDipoles.cpp")
target_link_libraries("UnitTest-eDipoles.o" ${LIBS})
target_link_libraries("UnitTest-eDipoles.o" Cubature)
add_test(NAME "UnitTest-eDipoles" COMMAND "UnitTest-eDipoles.o")

add_executable("UnitTest-SalaminTightlyFocused.o" "UnitTest-SalaminTightlyFocused.cpp")
target_link_libraries("UnitTest-SalaminTightlyFocused.o" ${LIBS} ${CUBALIBS})
target_link_libraries("UnitTest-SalaminTightlyFocused.o" Cubature)
target_link_libraries("UnitTest-SalaminTightlyFocused.o" ${cuba_dir}/libcuba.a)
add_test(NAME "UnitTest-SalaminTightlyFocused" COMMAND "UnitTest-SalaminTightlyFocused.o")


