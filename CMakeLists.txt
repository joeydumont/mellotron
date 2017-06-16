# --------------------------------------------------------------------------- #
# Author:       Denis Gagnon                  <denis.gagnon@emt.inrs.ca>      #
# Date created: 2017-06-13                                                    #
# Description:  CMake compilation instructions for MELLOTRON                  #
# ----------------------------------------------------------------------------#

# ----------------------------------------------------------------- #
# --                 Name and version of library                 -- #
# ----------------------------------------------------------------- #
project(mellotron)
set (mellotron_VERSION_MAJOR 0)
set (mellotron_VERSION_MINOR 5)
set (mellotron_VERSION_RELEASE 0)

# ----------------------------------------------------------------- #
# --               Configuration and Dependencies                -- #
# ----------------------------------------------------------------- #
# -- CMake version and installation directory.
# CMake version
cmake_minimum_required(VERSION 3.1)

# C++14 Standard required
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX /usr)
endif()
LIST (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
MESSAGE( STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH})

# -- Required dependency: HDF5.
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
set (LIBS ${LIBS} ${HDF5_LIBRARIES})

# -- Required dependency: Boost.
find_package(Boost 1.48.0 COMPONENTS program_options REQUIRED )
include_directories(${Boost_INCLUDE_DIRS})
SET (LIBS ${LIBS} ${Boost_LIBRARIES} )

# -- Required dependency: Armadillo
find_package(armadillo REQUIRED)
include_directories(${armadillo_INCLUDE_DIRS})
set(LIBS ${LIBS} ${armadillo_LIBRARIES})

# -- Required dependency: Cubature (can be in submodules, i.e. external/Cubature)
add_subdirectory (external/Cubature/)
include_directories (external/Cubature/include)

#find_package(Cubature REQUIRED)
#include_directories(${Cubature_INCLUDE_DIRS})
#set (CUBALIBS ${CUBALIBS} ${Cubature_LIBRARIES})
#message(STATUS "Cubature include directory: " ${Cubature_INCLUDE_DIR})
#message(STATUS "Cubature library: " ${Cubature_LIBRARY})

# -- Required dependency: cuba (can be in submodules, i.e. external/Cuba)
find_package(cuba)

if (${cuba_FOUND})
    # Do nothing
else()
    # Build cuba and then attempt to find it again
    message(STATUS "Building cuba in submodules directory")
    execute_process(COMMAND ./configure
                    COMMAND make
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external/Cuba)
endif()

find_package(cuba REQUIRED)
include_directories(${cuba_INCLUDE_DIRS})
set (CUBALIBS ${CUBALIBS} ${cuba_LIBRARIES})
message(STATUS "cuba.h include directory: " ${cuba_INCLUDE_DIR})
message(STATUS "cuba library: " ${cuba_LIBRARY})

# -- Compiler flags
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -O3")
set (CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -fPIC -O3")


# ----------------------------------------------------------------- #
# --                  Compilation Instructions                   -- #
# ----------------------------------------------------------------- #
# -- Included files
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include)

# -- Install mellotron headers
install (DIRECTORY include/          DESTINATION include)

# -- Output binaries in directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/simulations)

# -- Rules for every target

add_executable("IntegrationSalamin.o" "simulations/IntegrationSalamin.cpp")
target_link_libraries("IntegrationSalamin.o" ${LIBS} ${CUBALIBS})
target_link_libraries("IntegrationSalamin.o" Cubature)
install(TARGETS "IntegrationSalamin.o" RUNTIME DESTINATION bin)

add_executable("ComputeNormalizationConstantSalaminLinear.o" "simulations/ComputeNormalizationConstantSalaminLinear.cpp")
target_link_libraries("ComputeNormalizationConstantSalaminLinear.o" ${LIBS} ${CUBALIBS})
target_link_libraries("ComputeNormalizationConstantSalaminLinear.o"  Cubature)
install(TARGETS "ComputeNormalizationConstantSalaminLinear.o" RUNTIME DESTINATION bin)

