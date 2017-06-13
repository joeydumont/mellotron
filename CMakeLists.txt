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
# CMake version
cmake_minimum_required(VERSION 3.1)

# C++11 Standard required
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

# -- Required dependency: Cubature (should be in submodules)
find_package(Cubature REQUIRED)
include_directories(${Cubature_INCLUDE_DIRS})
set (CUBALIBS ${CUBALIBS} ${Cubature_LIBRARIES})
message(STATUS "Cubature include directory: " ${Cubature_INCLUDE_DIR})
message(STATUS "Cubature library: " ${Cubature_LIBRARY})


# -- Required dependency: cuba (should be in submodules)
find_package(cuba REQUIRED)
include_directories(${cuba_INCLUDE_DIRS})
set (CUBALIBS ${CUBALIBS} ${cuba_LIBRARIES})
message(STATUS "cuba.h include directory: " ${cuba_INCLUDE_DIR})
message(STATUS "cuba library: " ${cuba_LIBRARY})

# -- Required dependency: Armadillo
find_package(armadillo REQUIRED)
include_directories(${armadillo_INCLUDE_DIRS})
set(LIBS ${LIBS} ${armadillo_LIBRARIES})
message(STATUS "armadillo include directory: " ${armadillo_INCLUDE_DIR})

# -- Compiler flags
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -O3")
set (CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -fPIC -O3")


# ----------------------------------------------------------------- #
# --                  Compilation Instructions                   -- #
# ----------------------------------------------------------------- #
# -- Included files
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include)

# -- Install mellotron headers
install (DIRECTORY include/                  DESTINATION include)

# Output binaries in directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable("IntegrationSalamin" "simulations/IntegrationSalamin.cpp")
target_link_libraries("IntegrationSalamin" LINK_PUBLIC ${LIBS})
install(TARGETS IntegrationSalamin RUNTIME DESTINATION bin)

# ----------------------------------------------------------------- #
# --                   Generate Documentation                    -- #
# ----------------------------------------------------------------- #
# add a target to generate API documentation with Doxygen
#find_package(Doxygen)
#if(DOXYGEN_FOUND AND DOXYGEN_VERSION VERSION_GREATER 1.8)
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
#add_custom_target(
#  doc ALL ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
#      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#      COMMENT "Generating API documentation with Doxygen" VERBATIM
#)
#endif(DOXYGEN_FOUND AND DOXYGEN_VERSION VERSION_GREATER 1.8)
