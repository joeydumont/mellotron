# --------------------------------------------------------------------------- #
# Author:       Denis Gagnon                  <denis.gagnon@emt.inrs.ca>      #
#               Joey Dumont                   <joey.dumontn@gmail.com>        #
# Date created: 2017-06-13                                                    #
# Description:  CMake compilation instructions for MELLOTRON                  #
# ----------------------------------------------------------------------------#

# ----------------------------------------------------------------- #
# --                 Name and version of library                 -- #
# ----------------------------------------------------------------- #
project(mellotron)
set (mellotron_VERSION_MAJOR 1)
set (mellotron_VERSION_MINOR 1)
set (mellotron_VERSION_RELEASE 0)

# ----------------------------------------------------------------- #
# --               Configuration and Dependencies                -- #
# ----------------------------------------------------------------- #
# -- CMake version and installation directory.
# CMake version
cmake_minimum_required(VERSION 3.1)

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
message(STATUS "StrattoCalculator library: " ${STRATTOLIBS})
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

# -- Required dependency: GSL.
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

# -- Required dependency: Cubature (can be in submodules, i.e. external/Cubature)
add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/external/Cubature/)
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/external/Cubature/include)

# -- Required dependency: cuba (MUST be in submodules, i.e. external/Cuba)
include(ExternalProject)
ExternalProject_Add(
  libcuba.a
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/Cuba
  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/Cuba/configure --prefix=${CMAKE_CURRENT_SOURCE_DIR}/external/Cuba
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/external/Cuba
  BUILD_COMMAND make
  BUILD_IN_SOURCE 1
)

ExternalProject_Get_Property(libcuba.a install_dir)
set (cuba_dir ${install_dir})
include_directories(${cuba_dir}/include)

# ----------------------------------------------------------------- #
# --                    Compiler Configuration                   -- #
# ----------------------------------------------------------------- #
# -- Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# -- Macro definitions
if (CMAKE_BUILD_TYPE MATCHES RELEASE)
  # -- Remove bounds checking for performance.
  add_definitions(-DBOOST_DISABLE_ASSERTS)
  add_definitions(-DARMA_NO_DEBUG)

  # -- Tell GSL that the compiler supports inlining.
  add_definitions(-DHAVE_INLINE)
endif()

# C++14 Standard required
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Configuration for the GCC compiler.
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -O3")
  set (CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE}    -O3")

  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -pg -g -Wall -DNDEBUG")
  set (CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -O0 -pg -g -Wall -DNDEBUG")

  set (CMAKE_C_FLAGS_RELWITHDEBINFO   "${CMAKE_C_FLAGS_RELWITHDEBINFO}   -Wall -pg -O3 -DNDEBUG")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -Wall -pg -O3 -DNDEBUG")

# Configuration for the Intel compiler.
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -O3")
  set (CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE}    -O3")

  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -debug all")
  set (CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -O0 -g -debug all")

  set (CMAKE_C_FLAGS_RELWITHDEBINFO   "${CMAKE_C_FLAGS_RELWITHDEBINFO}   -O3 -g -debug all")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3 -g -debug all")
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

# -- Mellotron tests
include(CTest)
add_subdirectory(tests)

# -- Rules for every target
add_executable("IntegrationSalamin" "simulations/IntegrationSalamin.cpp")
target_link_libraries("IntegrationSalamin" ${LIBS} ${CUBALIBS})
target_link_libraries("IntegrationSalamin" Cubature)
target_link_libraries("IntegrationSalamin" ${cuba_dir}/libcuba.a)
install(TARGETS "IntegrationSalamin" RUNTIME DESTINATION bin)

add_executable("IntegrationSalaminIonized" "simulations/IntegrationSalaminIonized.cpp")
target_link_libraries("IntegrationSalaminIonized" ${LIBS} ${CUBALIBS})
target_link_libraries("IntegrationSalaminIonized" Cubature)
install(TARGETS "IntegrationSalaminIonized" RUNTIME DESTINATION bin)

add_executable("ComputeNormalizationConstantSalaminLinear" "simulations/ComputeNormalizationConstantSalaminLinear.cpp")
target_link_libraries("ComputeNormalizationConstantSalaminLinear" ${LIBS} ${CUBALIBS})
target_link_libraries("ComputeNormalizationConstantSalaminLinear"  Cubature)
target_link_libraries("ComputeNormalizationConstantSalaminLinear"  ${cuba_dir}/libcuba.a)
install(TARGETS "ComputeNormalizationConstantSalaminLinear" RUNTIME DESTINATION bin)

add_executable("IntegrationStrattoLinear" "simulations/IntegrationStrattoLinear.cpp")
target_link_libraries("IntegrationStrattoLinear" ${LIBS} ${CUBALIBS} ${STRATTOLIBS})
target_link_libraries("IntegrationStrattoLinear"  Cubature)
target_link_libraries("IntegrationStrattoLinear"  ${cuba_dir}/libcuba.a)
install(TARGETS "IntegrationStrattoLinear" RUNTIME DESTINATION bin)

add_executable("IntegrationStrattoMosaic" "simulations/IntegrationStrattoMosaic.cpp")
target_link_libraries("IntegrationStrattoMosaic" ${LIBS} ${CUBALIBS} ${STRATTOLIBS})
target_link_libraries("IntegrationStrattoMosaic"  Cubature)
target_link_libraries("IntegrationStrattoMosaic"  ${cuba_dir}/libcuba.a)
install(TARGETS "IntegrationStrattoMosaic" RUNTIME DESTINATION bin)

add_executable("IntegrationStrattoRadial" "simulations/IntegrationStrattoRadial.cpp")
target_link_libraries("IntegrationStrattoRadial" ${LIBS} ${CUBALIBS} ${STRATTOLIBS})
target_link_libraries("IntegrationStrattoRadial"  Cubature)
target_link_libraries("IntegrationStrattoRadial"  ${cuba_dir}/libcuba.a)
install(TARGETS "IntegrationStrattoRadial" RUNTIME DESTINATION bin)

add_executable("IntegrationStrattoLinearSG" "simulations/IntegrationStrattoLinearSG.cpp")
target_link_libraries("IntegrationStrattoLinearSG" ${LIBS} ${CUBALIBS} ${STRATTOLIBS})
target_link_libraries("IntegrationStrattoLinearSG"  Cubature)
target_link_libraries("IntegrationStrattoLinearSG"  ${cuba_dir}/libcuba.a)
install(TARGETS "IntegrationStrattoLinearSG" RUNTIME DESTINATION bin)
