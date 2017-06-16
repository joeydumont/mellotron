#!/bin/bash

# --------------------------------------------------------------------------- #
#  \file compile_test.sh                                                      #
#  \author Joey Dumont                  <joey.dumont@gmail.com>               #
#  \since 2015-10-19                                                          #
#                                                                             #
# This bash file compiles a given unit test of the MELLOTRON.                 #
#                                                                             #
# Usage: compile_test.sh <test-to-compile>                                    #
# --------------------------------------------------------------------------- #

# -- We first check that the script is given only a single argument.
if [[ $# -ne 1 ]]; then
  echo "Usage: compile_test.sh [-s] <test-to-compile>."
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
      -o ${filenameNoExt}.o               \
      ${filename}                         \
      -lhdf5 -lgtest -lCubature -lcuba -larmadillo
exit 0
