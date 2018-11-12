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
#     bash build.sh (release OR debug)                                        #
#                                                                             #
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


PROJECT_DIR=~/projects/rrg-maclean-ab/maclean_group/
module use ${PROJECT_DIR}/modules
module load openmpi
module load hdf5-mpi
module load meshpi
module load boost-mpi
module load muparser
module load fftw-mpi
module load gtest
module load armadillo
module load strattocalculator
module load gsl

export CMAKE_LIBRARY_PATH=$LIBRARY_PATH
export CMAKE_INCLUDE_PATH=$INCLUDE_PATH

CMAKE_FLAGS="${CMAKE_FLAGS} -DCMAKE_INSTALL_PREFIX=${PROJECT_DIR}/software/mellotron/1.1.0"

cmake ${CMAKE_FLAGS} ..

make

exit 0
