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
add_executable("UnitTest-Particle" "UnitTest-Particle.cpp")
target_link_libraries("UnitTest-Particle" ${LIBS})
add_test(NAME "UnitTest-Particle" COMMAND "UnitTest-Particle")

add_executable("UnitTest-eDipoles" "UnitTest-eDipoles.cpp")
target_link_libraries("UnitTest-eDipoles" ${LIBS})
target_link_libraries("UnitTest-eDipoles" Cubature)
add_test(NAME "UnitTest-eDipoles" COMMAND "UnitTest-eDipoles")

add_executable("UnitTest-ParticleObserverLienardWiechert" "UnitTest-ParticleObserverLienardWiechert.cpp")
target_link_libraries("UnitTest-ParticleObserverLienardWiechert" ${LIBS})
add_test(NAME "UnitTest-ParticleObserverLienardWiechert" COMMAND "UnitTest-ParticleObserverLienardWiechert")

add_executable("UnitTest-SalaminTightlyFocused" "UnitTest-SalaminTightlyFocused.cpp")
target_link_libraries("UnitTest-SalaminTightlyFocused" ${LIBS} ${CUBALIBS})
target_link_libraries("UnitTest-SalaminTightlyFocused" Cubature)
target_link_libraries("UnitTest-SalaminTightlyFocused" Cuba)
add_test(NAME "UnitTest-SalaminTightlyFocused" COMMAND "UnitTest-SalaminTightlyFocused" WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/tests/")

if (BUILD_STRATTO_DRIVERS)
  add_executable("UnitTest-StrattoCalculatorWrapper" "UnitTest-StrattoCalculatorWrapper.cpp")
  target_link_libraries("UnitTest-StrattoCalculatorWrapper" ${LIBS} ${STRATTODEPS})
  add_test(NAME "UnitTest-StrattoCalculatorWrapper" COMMAND "UnitTest-StrattoCalculatorWrapper")
endif()
