#!/bin/bash

# ------------------------------------------------------------------------------------- #
#  \file automaticStart.sh                                                              #
#  \author Justine Pepin                  <justine.pepin@emt.inrs.ca>                   #
#  \since 2017-05-19                                                                    #
#                                                                                       #
# This bash file start automatically the simulation of the MELLOTRON.                   #
#                                                                                       #
# Usage: ./automaticStart.sh -d <dirname> -c <name of config.xml> -s <shape> -j <njobs> #
#                             -t <BOOL:varInitTime>                                     #
#        Where dirname is the name of the directory containing a configSalamin.xml      #
# ------------------------------------------------------------------------------------- #

echo -e " \e[95m--- MELLOTRON SIMULATION AUTOMATIC START ---\e[39m"
DIR="./"
CONFIG="configSalamin.xml"
SHAPE="sphere"
NJOBS="8"
VARINITTIME=false
while getopts ":d:c:s:j:" opt; do
    case $opt in
        d)
            DIR=$OPTARG
            ;;
        c)
            CONFIG=$OPTARG
            ;;
        s)
            SHAPE=$OPTARG
            ;;
        j)
            NJOBS=$OPTARG
            ;;
        t)
            VARINITTIME=$OPTARG
            ;;
    esac
done

# Initialize the variables
GENINIT="GenerateInitialConditions.py"
INTEGSAL="IntegrationSalamin.o"
INTEGSTRATTOLIN="IntegrationStrattoLinear.o"
INTEGSTRATTOMOS="IntegrationStrattoMosaic.o"
INTEGSTRATTORAD="IntegrationStrattoRadial.o"
COMPNORMCONST="ComputeNormalizationConstantSalaminLinear.o"
MANAGEOUT="manageOutputs.py"
PRODUCEPLOTS="producePlots.py"
cp $GENINIT ./$DIR
cp $INTEGSAL ./$DIR
cp $INTEGSTRATTOLIN ./$DIR
cp $INTEGSTRATTOMOS ./$DIR
cp $INTEGSTRATTORAD ./$DIR
cp $COMPNORMCONST ./$DIR
cp $MANAGEOUT ./$DIR
cp $PRODUCEPLOTS ./$DIR
cd ./$DIR
OUTINITCONDS="init_conds.txt"
OUTNORMCONST="normalization_constant.txt"

# Check if config.xml is in current dir
if [ -f  ./$CONFIG ]; then
    echo -e " \e[32m--- Config file has been found. ---\e[39m"
else
    echo -e " \e[32m--- missing config.xml file. ---\e[39m"
    rm $GENINIT $INTEGSAL $INTEGSTRATTOLIN $INTEGSTRATTORAD $INTEGSTRATTOMOS $COMPNORMCONST $MANAGEOUT
    echo "Mellotron can not be run. Exiting. "
    exit 0
fi

# Generate initial conditions
if [ -f  ./$OUTINITCONDS ]; then
    echo -e " \e[32m--- Initial conditions has been found. ---\e[39m"
else
    echo -e " \e[32m--- Starting to generate initial conditions. ---\e[39m"
    python $GENINIT --shape $SHAPE --config $CONFIG
    echo "Done: generate initial conditions."
fi

CONFIGDEFAULT="configSalamin.xml"
if [ "$CONFIG" == "$CONFIGDEFAULT" ]; then
    # Compute normalization constant
    if [ -f  ./$OUTNORMCONST ]; then
        echo -e " \e[32m--- Normalization constant has been found. ---\e[39m"
    else
        echo -e " \e[32m--- Starting to generate normalization constant. ---\e[39m"
        ./$COMPNORMCONST
        echo "Done: compute normalization constant."
    fi
    INTEG=$INTEGSAL
elif [ "$CONFIG" == "configStrattoLinear.xml" ]; then
    INTEG=$INTEGSTRATTOLIN
elif [ "$CONFIG" == "configStrattoMosaic.xml" ]; then
    INTEG=$INTEGSTRATTOMOS
elif [ "$CONFIG" == "configStrattoRadial.xml" ]; then
    INTEG=$INTEGSTRATTORAD
fi

# -- Calculate particles behavior
echo -e " \e[32m--- Starting to calculate particles behavior. ---\e[39m"
if [[ "$VARINITTIME" == true ]]; then
    cat $OUTINITCONDS | parallel -j $NJOBS --colsep " " ./$INTEG --init_conds {1} {2} {3} {4} {5} {6}
else
    cat $OUTINITCONDS | parallel -j $NJOBS --colsep " " ./$INTEG --init_conds {1} {2} {3} {4} {5} {6} {7}
fi
echo "Done: calculate particles behavior."

# -- Manage outputs
echo -e " \e[32m--- Starting to manage the outputs. ---\e[39m"
NUMBER=$(ls -d *.hdf5 | wc -l)
python $MANAGEOUT --directory ./
echo "Done: manage outputs."

# -- Generate plots
echo -e " \e[32m--- Starting to produce the plots ---\e[39m"
python $PRODUCEPLOTS --directory ./
echo "Done: produce plots."

# -- Generate PDFs.
echo -e " \e[32m--- Convert EPS plots to PDF ---\e[39m"
for file in *.eps
do
    epstopdf $file
done
echo "Done: convert EPS files."

exit 0
