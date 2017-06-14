#!/bin/bash

# --------------------------------------------------------------------------- #
#  \file compileAll.sh                                                        #
#  \author Justine Pepin                  <justine.pepin@emt.inrs.ca>         #
#  \since 2017-05-26                                                          #
#                                                                             #
# This bash file compile all required cpp files for the MELLOTRON.            #
#                                                                             #
# Usage: ./compileAll.sh                                                      #
# --------------------------------------------------------------------------- #

# -- Check if compile.sh exists
COMP="compile.sh"
if [ ! -f $COMP ]; then
echo "Usage: compileAll.sh needs compile.sh in same directory."
exit 1
fi
chmod +x compile.sh

# -- Compile IntegrationSalamin.cpp to have .o
INTSALCPP="IntegrationSalamin.cpp"
if [ ! -f $INTSALCPP ]; then
echo "Usage: compileAll.sh needs IntegrationSalamin.cpp in same directory."
exit 1
fi
./compile.sh $INTSALCPP
echo "Done: compile IntegrationSalamin.cpp."

# -- Compile ComputeNormalizationConstantsSalaminLinear.cpp to have .o
COMPUTECPP="ComputeNormalizationConstantSalaminLinear.cpp"
if [ ! -f $COMPUTECPP ]; then
echo "Usage: compileAll.sh needs ComputeNormalizationConstantSalaminLinear.cpp in same directory."
exit 1
fi
./compile.sh $COMPUTECPP
echo "Done: compile ComputeNormalizationConstantSalaminLinear.cpp."

# -- Check that the script can find IntegrationSalamin.o and ComputeNormalizationConstantSalaminLinear.o
INTSAL="IntegrationSalamin.o"
COMPUTE="ComputeNormalizationConstantSalaminLinear.o" 
if [ ! -f $INTSAL ]; then
echo "Usage: compileAll.sh needs IntegrationSalamin.o in same directory."
exit 1
fi
if [ ! -f $COMPUTE ]; then
echo "Usage: compileAll.sh needs ComputeNormalizationConstantSalaminLinear.o in same directory."
exit 1
fi



