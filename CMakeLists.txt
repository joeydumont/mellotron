# --------------------------------------------------------------------------- #
# Author:       Joey Dumont                   <joey.dumont@gmail.com>         #
# Date created: 2016-09-29                                                    #
# Description:  CMake compilation instructions for MELLOTRON.                 #
# ----------------------------------------------------------------------------#

# ----------------------------------------------------------------- #
# --                 Name and version of library                 -- #
# ----------------------------------------------------------------- #
project(mellotron)
set (mellotron_VERSION_MAJOR 0)
set (mellotron_VERSION_MINOR 0)
set (mellotron_VERSION_RELEASE 0)

# ----------------------------------------------------------------- #
# --               Configuration and Dependencies                -- #
# ----------------------------------------------------------------- #
# -- CMake version and installation directory.
cmake_minimum_required(VERSION 2.8)
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
find_package(Boost 1.48.0 REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})
SET (LIBS ${LIBS} ${BOOST_LIBRARIES})

# -- Required dependency: Cubature.
find_package(Cubature REQUIRED)
include_directories(${CUBATURE_INCLUDE_DIRS})
SET (LIBS ${LIBS} ${CUBATURE_LIBRARIES})

# -- Compiler config for RELEASE build.
# Configuration for the GCC compiler.
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  MESSAGE( STATUS "We use the GNU compiler settings.")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -O3 -Wall -march=native -std=c++14")
  set (CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -fPIC -O3 -Wall -march=native")
endif()

# Configuration for the Intel compiler.
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++14")
  set (CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -O3")
endif()

# -- Compiler config for DEBUG build.
# Configuration for the GCC compiler.
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  MESSAGE( STATUS "We use the GNU compiler settings.")
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -fPIC -O0 -pg -g -Wall -march=native -std=c++14")
  set (CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS}   -fPIC -O0 -pg -g -Wall -march=native")
endif()

# Configuration for the Intel compiler.
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g -debug all -std=c++14")
  set (CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS}   -O0 -g -debug all")
endif()

# ----------------------------------------------------------------- #
# --                  Compilation Instructions                   -- #
# ----------------------------------------------------------------- #
# -- Included files
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include INC_LIST)

install (DIRECTORY include/ DESTINATION include)

# ----------------------------------------------------------------- #
# --                   Generate Documentation                    -- #
# ----------------------------------------------------------------- #
# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND AND DOXYGEN_VERSION VERSION_GREATER 1.8)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(
  doc ALL ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND AND DOXYGEN_VERSION VERSION_GREATER 1.8)
