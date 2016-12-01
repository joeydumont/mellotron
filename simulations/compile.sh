#!/bin/bash

# --------------------------------------------------------------------------- #
#  \file compile.sh                                                           #
#  \author Joey Dumont                  <joey.dumont@gmail.com>               #
#  \since 2015-10-19                                                          #
#                                                                             #
# This bash file compiles a given simulation program of the MELLOTRON.        #
#                                                                             #
# Usage: compile.sh <test-to-compile>                                         #
# --------------------------------------------------------------------------- #

# -- We first check that the script is given only a single argument.
if [[ $# -ne 1 ]]; then
  echo "Usage: compile.sh [-s] <program-to-compile>."
  exit 1
fi

# -- We check if the test exists.
if [ ! -f $1 ]; then
  echo "This test does not exist."
  exit 1
fi

# -- We prepare the compilation of the test.
filename=`basename $1`
filenameNoExt="${filename%.*}"

# -- We compile the test.
CXX_FLAGS="-Wall -std=c++14 -O3"
g++ ${CXX_FLAGS}                          \
      -I ../include                       \
      -o ${filenameNoExt}                 \
      ${filename}                         \
      -lhdf5 -lboost_program_options -lgtest -larmadillo
      
exit 0
