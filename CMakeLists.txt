# --------------------------------------------------------------------------- #
# Author:       Denis Gagnon                  <denis.gagnon@emt.inrs.ca>      #
# Date created: 2017-06-13                                                    #
# Description:  CMake compilation instructions for MELLOTRON                  #
# ----------------------------------------------------------------------------#

# ----------------------------------------------------------------- #
# --                 Name and version of library                 -- #
# ----------------------------------------------------------------- #
project(mellotron)
set (mellotron_VERSION_MAJOR 1)
set (mellotron_VERSION_MINOR 0)
set (mellotron_VERSION_RELEASE 1)

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

# -- Compiler flags
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -O3")
set (CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -fPIC -O3")

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX /usr)
endif()
LIST (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
MESSAGE( STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH})

# -- Required dependency: HDF5.
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
set (LIBS ${LIBS} ${HDF5_LIBRARIES})

# -- Required dependency: MeshPI.
find_package(MeshPI REQUIRED)
include_directories(${meshpi_INCLUDE_DIR})
set (LIBS ${LIBS} ${meshpi_LIBRARY})
message(STATUS "MeshPI library: " ${meshpi_LIBRARIES})

# -- Required dependency: StrattoCalculator.
find_package(StrattoCalculator REQUIRED)
include_directories(${strattocalculator_INCLUDE_DIR})
set (STRATTOLIBS ${STRATTOLIBS} ${strattocalculator_LIBRARIES})

# -- Required dependency: Boost.
find_package(Boost 1.58.0 COMPONENTS program_options REQUIRED )
include_directories(${Boost_INCLUDE_DIRS})
SET (LIBS ${LIBS} ${Boost_LIBRARIES} )

# -- Required dependency: Armadillo
find_package(armadillo REQUIRED)
include_directories(${armadillo_INCLUDE_DIRS})
set(LIBS ${LIBS} ${armadillo_LIBRARIES})

# -- Required dependency: MPI
find_package(MPI REQUIRED)
set(LIBS ${LIBS} ${MPI_LIBRARIES})

# -- Required dependency: Cubature (can be in submodules, i.e. external/Cubature)
add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/external/Cubature/)
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/external/Cubature/include)

# -- Required dependency: cuba (MUST be in submodules, i.e. external/Cuba)
include(ExternalProject)
ExternalProject_Add(
  libcuba.a
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/Cuba
  CONFIGURE_COMMAND ${CMAKE_SOURCE_DIR}/external/Cuba/configure --prefix=${CMAKE_CURRENT_SOURCE_DIR}/external/Cuba
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/external/Cuba
  BUILD_COMMAND make
  BUILD_IN_SOURCE 1
)

ExternalProject_Get_Property(libcuba.a install_dir)
set (cuba_dir ${install_dir})
include_directories(${cuba_dir}/include)

# -- Mellotron tests
include(CTest)
add_subdirectory(tests)

# -- Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# -- Compiler config for RELEASE build.
# Configuration for the GCC compiler.
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -O3")
  set (CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE}    -O3")

elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -O3")
  set (CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE}    -O3")
endif()

# -- Compiler config for DEBUG build.
# Configuration for the GCC compiler.
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -pg -g -Wall")
  set (CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -O0 -pg -g -Wall")

# Configuration for the Intel compiler.
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -debug all")
  set (CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -O0 -g -debug all")
endif()

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
target_link_libraries("IntegrationSalamin.o" ${cuba_dir}/libcuba.a)
install(TARGETS "IntegrationSalamin.o" RUNTIME DESTINATION bin)

add_executable("ComputeNormalizationConstantSalaminLinear.o" "simulations/ComputeNormalizationConstantSalaminLinear.cpp")
target_link_libraries("ComputeNormalizationConstantSalaminLinear.o" ${LIBS} ${CUBALIBS})
target_link_libraries("ComputeNormalizationConstantSalaminLinear.o"  Cubature)
target_link_libraries("ComputeNormalizationConstantSalaminLinear.o"  ${cuba_dir}/libcuba.a)
install(TARGETS "ComputeNormalizationConstantSalaminLinear.o" RUNTIME DESTINATION bin)

add_executable("IntegrationStrattoLinear.o" "simulations/IntegrationStrattoLinear.cpp")
target_link_libraries("IntegrationStrattoLinear.o" ${LIBS} ${CUBALIBS} ${STRATTOLIBS})
target_link_libraries("IntegrationStrattoLinear.o"  Cubature)
target_link_libraries("IntegrationStrattoLinear.o"  ${cuba_dir}/libcuba.a)
install(TARGETS "IntegrationStrattoLinear.o" RUNTIME DESTINATION bin)
