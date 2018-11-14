#!/bin/bash

# Print out the SLURM job information. Remove this if you don't need it.
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR

# -- Project directory.
PROJECT_DIR=~/projects/rrg-maclean-ab/maclean_group/modules/

# -- We import the proper modules.
module use ${PROJECT_DIR}
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

# -- We run the program.
srun <executable_in_path> <figure_out_what_to_put_here>
