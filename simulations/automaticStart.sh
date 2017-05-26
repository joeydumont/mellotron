#!/bin/bash

# --------------------------------------------------------------------------- #
#  \file automaticStart.sh                                                    #
#  \author Justine Pepin                  <justine.pepin@emt.inrs.ca>         #
#  \since 2017-05-19                                                          #
#                                                                             #
# This bash file start automatically the simulation of the MELLOTRON.         #
#                                                                             #
# Usage: automaticStart.sh                                                    #
# --------------------------------------------------------------------------- #

echo -e " \e[95m--- MELLOTRON SIMULATION AUTOMATIC START ---\e[39m"
read -r -p " Do you want to compile everything before starting simulation? [y/N] " response
    if [[ "$response" =~ ^([yY][eE][sS]|[yY])+$ ]]; then
        echo -e " \e[32m--- Starting compilation. ---\e[39m"

        # -- Check if compile.sh exists
        COMP="compile.sh"
        if [ ! -f $COMP ]; then
        echo "Usage: automaticStart.sh needs compile.sh in same directory."
        exit 1
        fi
        chmod +x compile.sh

        # -- Compile IntegrationSalamin.cpp to have .o
        INTSALCPP="IntegrationSalamin.cpp"
        if [ ! -f $INTSALCPP ]; then
        echo "Usage: automaticStart.sh needs IntegrationSalamin.cpp in same directory."
        exit 1
        fi
        ./compile.sh $INTSALCPP
        echo "Done: compile IntegrationSalamin.cpp."

        # -- Compile ComputeNormalizationConstantsSalaminLinear.cpp to have .o
        COMPUTECPP="ComputeNormalizationConstantSalaminLinear.cpp"
        if [ ! -f $COMPUTECPP ]; then
        echo "Usage: automaticStart.sh needs ComputeNormalizationConstantSalaminLinear.cpp in same directory."
        exit 1
        fi
        ./compile.sh $COMPUTECPP
        echo "Done: compile ComputeNormalizationConstantSalaminLinear.cpp."

        # -- Check that the script can find IntegrationSalamin.o, GenerateInitialConditions.py and ComputeNormalizationConstantSalaminLinear.o
        INTSAL="IntegrationSalamin.o"
        GENINIT="GenerateInitialConditions.py"
        COMPUTE="ComputeNormalizationConstantSalaminLinear.o" 
        if [ ! -f $INTSAL ]; then
        echo "Usage: automaticStart.sh needs IntegrationSalamin.o in same directory."
        exit 1
        fi
        if [ ! -f $GENINIT ]; then
        echo "Usage: automaticStart.sh needs GenerateInitialConditions.py in same directory."
        exit 1
        fi
        if [ ! -f $COMPUTE ]; then
        echo "Usage: automaticStart.sh needs ComputeNormalizationConstantSalaminLinear.o in same directory."
        exit 1
        fi

        # -- Compute normalization constant
        OUTNORMCONST="normalization_constant.txt"
        ./ComputeNormalizationConstantSalaminLinear.o --outfile $OUTNORMCONST
        echo "Done: compute normalization constant."

        echo -e " \e[32m--- Everything is compiled. Starting to generate initial conditions. ---\e[39m"

    else
        GENINIT="GenerateInitialConditions.py"
        OUTNORMCONST="normalization_constant.txt"
        echo -e " \e[32m--- Starting to generate initial conditions. ---\e[39m"

    fi


# -- Generate initial conditions
OUTINITCONDS="init_conds.txt"
read -r -p " Do you want to use default values? [y/N] " response
    if [[ "$response" =~ ^([yY][eE][sS]|[yY])+$ ]]; then
        python $GENINIT --outfile $OUTINITCONDS
        NUMBER="10"
    else
        # -- Number of initial conditions
        NUMBER=""
        while [[ ! $NUMBER =~ ^[0-9]+$ ]]; do
            echo "Please enter the number of initial conditions (integer)"
            read NUMBER
        done
        # -- Wavelength
        LENGTH=""
        while [[ ! $LENGTH =~ ^[0-9]+\.?[0-9]*$ ]]; do
            echo "Please enter the wavelength in m (float)"
            read LENGTH
        done
        # -- Momentum in z
        MOMENTUM=""
        while [[ ! $MOMENTUM =~ ^[-+]?[0-9]+\.?[0-9]*$ ]]; do
            echo "Please enter the initial value of momentum in z in eV/c (float)"
            read MOMENTUM
        done
        # -- Radius
        RADIUS=""
        while [[ ! $RADIUS =~ ^[0-9]+\.?[0-9]*$ ]]; do
            echo "Please enter the sphere radius in m (float)"
            read RADIUS
        done
        python $GENINIT --outfile $OUTINITCONDS --wavelength $LENGTH --pz $MOMENTUM --numpart $NUMBER --radius $RADIUS
fi
DIRNAME=$(date +%F--%T)
mkdir $DIRNAME
mv $OUTINITCONDS $DIRNAME
echo "Done: generate initial conditions."


# -- Calculate particles behavior
echo -e " \e[32m--- Starting to calculate particles behavior. ---\e[39m"
cd $DIRNAME
while IFS='' read -r LINE || [[ -n "$LINE" ]]; do
    ../IntegrationSalamin.o --norm_constant_file ../$OUTNORMCONST --init_conds $LINE
done < $OUTINITCONDS
cd ..
echo "Done: calculate particles behavior."

# -- Calculate particles behavior parallel
# -- echo -e " \e[32m--- Starting to calculate particles behavior. ---\e[39m"
# -- cd $DIRNAME
# -- while IFS='' read -r LINE || [[ -n "$LINE" ]]; do
# -- parallel -p 2 ../IntegrationSalamin.o --norm_constant_file {../$OUTNORMCONST} --init_conds $LINE {} ::: #1 ::: #2  
# -- done < $OUTINITCONDS
# -- cd ..
# -- echo "Done: calculate particles behavior."


# -- Manage outputs
echo -e " \e[32m--- Starting to manage the outputs. ---\e[39m"
python manageOutputs.py --nParticles $NUMBER --directory $DIRNAME
echo "Done: manage outputs."
exit 0