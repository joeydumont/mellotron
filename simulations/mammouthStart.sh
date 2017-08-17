
#!/bin/bash


#PBS -l walltime=120:00:00
#PBS -l nodes=16:ppn=1
#PBS -r n


#-------------------------------------------------------------#
# Module loading                                              #
#-------------------------------------------------------------#
module use /home/maclean_group/modulefiles/
module load gcc/5.2.0
module load openmpi/1.10.0_gcc
module load hdf5/1.8.15p1_openmpi_gcc5
module load python64/2.7.3
source /home/maclean_group/python_environment_openmpi/bin/activate
module load boost64
module load cmake/3.3.1
module load muparser/2.2.3
module load meshpi/1.1.0
module load jsoncpp/1.6.2
module load zernike/0.0.1
module load strattocalculator/3.0.2
module load armadillo
module load cubature
module load parallel
module load latex
module load mellotron
PYPATH=$( echo $PATH | grep -o '[^:]*mellotron[^:]*' )

#--------------------------------------------------------------#
# Arguments for the simulation                                 #
#--------------------------------------------------------------#
NNODES=16
CONFIG="configStrattoLinear.xml"
SHAPE="sphere"


# For mammouth : All processors of required nodes have their job.
NJOBS=$((24*$NNODES))
# To be sure the outputs will be in the scratch.
cd "${PBS_O_WORKDIR}"

# ------------------------------------------------------------ #
#  \file mammouthStart.sh                                      #
#  \author Justine Pepin         <justine.pepin@emt.inrs.ca>   #
#  \since 2017-08-11                                           #
#                                                              #
# This bash file start automatically the simulation of         #
# the MELLOTRON on the mammouth computer.                      #
#                                                              #
# Usage: qsub mammouthStart.sh                                 #
#                                                              #
# ------------------------------------------------------------ #

echo -e " \e[95m--- MELLOTRON SIMULATION AUTOMATIC START ---\e[39m"

# -- Initialize the variables
DIR="./"
GENINIT="GenerateInitialConditions.py"
INTEGSAL="IntegrationSalamin.o"
INTEGSTRATTOLIN="IntegrationStrattoLinear.o"
COMPNORMCONST="ComputeNormalizationConstantSalaminLinear.o"
MANAGEOUT="manageOutputs.py"
PRODUCEPLOTS="producePlots.py"
cp $PYPATH/$GENINIT $DIR
cp $PYPATH/$INTEGSAL $DIR
cp $PYPATH/$INTEGSTRATTOLIN $DIR
cp $PYPATH/$COMPNORMCONST $DIR
cp $PYPATH/$MANAGEOUT $DIR
cp $PYPATH/$PRODUCEPLOTS $DIR

OUTINITCONDS="init_conds.txt"
OUTNORMCONST="normalization_constant.txt"

# -- Check if config.xml is in current dir
if [ -f  ./$CONFIG ]; then
    echo -e " \e[32m--- Config file has been found. ---\e[39m"
else 
    echo -e " \e[32m--- missing config.xml file. ---\e[39m"
    rm $GENINIT $INTEGSAL $INTEGSTRATTOLIN $COMPNORMCONST $MANAGEOUT $PRODUCEPLOTS
    echo "Mellotron can not be run. Exiting. "
    exit 0
fi

# -- Generate initial conditions
if [ -f  ./$OUTINITCONDS ]; then
    echo -e " \e[32m--- Initial conditions has been found. ---\e[39m"
else
    echo -e " \e[32m--- Starting to generate initial conditions. ---\e[39m"
    python $GENINIT --shape $SHAPE --config $CONFIG
    echo "Done: generate initial conditions."
fi

# -- Generate the normalization constant if needed
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
else
    INTEG=$INTEGSTRATTOLIN
fi

# -- Calculate particles behavior
echo -e " \e[32m--- Starting to calculate particles behavior. ---\e[39m"
. /home/maclean_group/software/parallel/20170622/bin/env_parallel.bash
cat $OUTINITCONDS | env_parallel --sshloginfile $PBS_NODEFILE --sshdelay 0.1  --workdir=$PWD  --jobs $NJOBS --colsep " " ./$INTEG --init_conds {1} {2} {3} {4} {5} {6}
echo "Done: calculate particles behavior."

# -- Manage outputs
echo -e " \e[32m--- Starting to manage the outputs. ---\e[39m"
python $MANAGEOUT --directory $DIR
echo "Done: manage outputs."

# -- Generate plots
echo -e " \e[32m--- Starting to produce the plots ---\e[39m"
python $PRODUCEPLOTS --directory $DIR
echo "Done: produce plots."

# -- Clean dir
rm $GENINIT $INTEGSAL $INTEGSTRATTOLIN $COMPNORMCONST $MANAGEOUT $PRODUCEPLOTS
exit 0


