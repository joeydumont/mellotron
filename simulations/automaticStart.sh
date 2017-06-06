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
GENINIT="GenerateInitialConditions.py"
echo -e " \e[32m--- Starting to generate initial conditions. ---\e[39m"

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
        while [[ ! $LENGTH =~ ^[0-9]+\.?[0-9]+\e?[-+]?[0-9]*$ ]]; do
            echo "Please enter the wavelength in meters with the scientific notation (ex: 8.0e-7)"
            read LENGTH
        done
        # -- Momentum in z
        MOMENTUM=""
        while [[ ! $MOMENTUM =~ ^[-+]?[0-9]+\.?[0-9]*$ ]]; do
            echo "Please enter the initial value of the z momentum in eV/c (float)"
            read MOMENTUM
        done
        # -- Radius
        RADIUS=""
        while [[ ! $RADIUS =~ ^[0-9]+\.?[0-9]+\e?[-+]?[0-9]*$ ]]; do
            echo "Please enter the sphere radius in meters with the scientific notation (ex: 3.0e-6)"
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
read -r -p " Do you want to use default values? [y/N] " response
    if [[ "$response" =~ ^([yY][eE][sS]|[yY])+$ ]]; then
        cat $OUTINITCONDS | parallel -j 2 --colsep " " ../IntegrationSalamin.o --init_conds {1} {2} {3} {4} {5} {6}
    else
        # -- Energy
        ENERGY=""
        while [[ ! $ENERGY =~ ^[-+]?[0-9]+\.?[0-9]*$ ]]; do
            echo "Please enter the pulse energy in joules (float)"
            read ENERGY
        done
        # -- Wavelength
        LAMBDA=""
        while [[ ! $LAMBDA =~ ^[0-9]+\.?[0-9]+\e?[-+]?[0-9]*$ ]]; do
            echo "Please enter the wavelength in meters with the scientific notation (ex: 8.0e-7)"
            read LAMBDA
        done
        # -- w0
        WAIST=""
        while [[ ! $WAIST =~ ^[0-9]+\.?[0-9]+\e?[-+]?[0-9]*$ ]]; do
            echo "Please enter the beam waist in units of lambda with the scientific notation (ex: 0.2e+5)"
            read WAIST
        done
        # -- L
        AXIALWAIST=""
        while [[ ! $AXIALWAIST =~ ^[0-9]+\.?[0-9]+\e?[-+]?[0-9]*$ ]]; do
            echo "Please enter the axial length of beam in units of lambda with the scientific notation (ex: 7e-3)"
            read AXIALWAIST
        done
        # -- Mass
        MASS=""
        while [[ ! $MASS =~ ^[-+]?[0-9]+\.?[0-9]*$ ]]; do
            echo "Please enter the particle mass in units of electron mass (float)"
            read MASS
        done
        # -- Q
        Q=""
        while [[ ! $Q =~ ^[-+]?[0-9]+\.?[0-9]*$ ]]; do
            echo "Please enter the particle charge in units of electron charge (float)"
            read Q
        done
        # -- T_init
        TINIT=""
        while [[ ! $TINIT =~ ^[0-9]+\.?[0-9]+\e?[-+]?[0-9]*$ ]]; do
            echo "Please enter the initial time of simulation in seconds with the scientific notation (ex: 6.9e-1)"
            read TINIT
        done
        # -- Dt
        DT=""
        while [[ ! $DT =~ ^[0-9]+\.?[0-9]+\e?[-+]?[0-9]*$ ]]; do
            echo "Please enter the duration of a time step in seconds with the scientific notation (ex: 8.1e-15)"
            read DT
        done
        # -- Number of time steps
        NSTEPS=""
        while [[ ! $NSTEPS =~ ^[0-9]+$ ]]; do
            echo "Please enter the number of time steps (integer)"
            read NSTEPS
        done
        cat $OUTINITCONDS | parallel -j 2 --colsep " " ../IntegrationSalamin.o --init_conds {1} {2} {3} {4} {5} {6} --energy $ENERGY --lam $LAMBDA --w0 $WAIST --L $AXIALWAIST --mass $MASS --Q $Q --t_init $TINIT --dt $DT --nsteps $NSTEPS
fi
cd ..
echo "Done: calculate particles behavior."

# -- Manage outputs
echo -e " \e[32m--- Starting to manage the outputs. ---\e[39m"
python manageOutputs.py --nParticles $NUMBER --directory $DIRNAME
echo "Done: manage outputs."
exit 0