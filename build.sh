# --------------------------------------------------------------------------- #
# Author:       Joey Dumont                   <joey.dumont@gmail.com>         #
# Date:         2016-10-18                                                    #
# Description:  Calls CMake from the build directory and compiles the library.#
# License:      CC0 - Public Domain                                           #
#                                                                             #
# This script simply calls CMake from the  build directory and compiles the   #
# library.                                                                    #
#                                                                             #
# Usage:                                                                      #
#     bash build.sh (release OR debug) [cluster]                              #
# where the optional argument cluster sets up the environment for compilation #
# in a given cluster.                                                         #
# --------------------------------------------------------------------------- #

# Check the number of arguments.
if [ $# -gt 2 ]; then
  printf "Usage is\n"
  printf "\tbash build.sh (release OR debug) [cluster]\n"
  exit 1
fi

# Change to file directory.
cd "$(dirname "$(readlink -f "$0")")";

# Check if build/ dir exists.
if [ ! -d build ]; then
    mkdir build
else
    rm -rf build
    mkdir  build
fi

# Change to build dir and compile the library.
cd build
CMAKE_FLAGS=
if [ $# -gt 0 ]; then
  if [ $1 == release ]; then
    CMAKE_FLAGS="${CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Release"
  elif [ $1 == debug ]; then
    CMAKE_FLAGS="${CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Debug"
  fi
fi

if [ $# -gt 1 ]; then
  if [ $2 == mammouth ]; then
    module purge

    if [ $? -ne 0 ]; then
      printf "module command could not be executed. Are you sure you are on a cluster?\n"
      #exit 1
    fi

    module use  /home/maclean_group/modulefiles/
    module load gcc/5.2.0
    module load boost64
    module load openmpi/1.10.0_gcc
    module load hdf5/1.8.15p1_openmpi_gcc5
    module load cmake/3.3.1
    module load muparser/2.2.3
    module load meshpi/1.1.0
    module load jsoncpp/1.6.2
    module load zernike/0.0.1
    module load strattocalculator/3.0.2
    module load armadillo

    export CMAKE_LIBRARY_PATH=$LIBRARY_PATH
    export CMAKE_INCLUDE_PATH=$INCLUDE_PATH
    export CC=gcc
    export CXX=g++
    CMAKE_FLAGS="${CMAKE_FLAGS} -DCMAKE_INSTALL_PREFIX=/home/maclean_group/software/mellotron/1.0.1"


  else
    printf "We do not know how to configure the build environment on this cluster.\n"
    printf "Available options are:\n"
    printf "\t-mammouth\n"
    exit 1
  fi
fi

cmake ${CMAKE_FLAGS} ..

make

exit 0
