#!/bin/bash

#SBATCH -J "HeatEqnOne_8"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=haswell
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=01:30:00
#SBATCH --exclusive

source ./tests/scaling/ZIH_TAURUS2_GCC_4_9_1/run-HeatEqnOne3D-timesteps_20-ADOLC_NO-MPI_ASYNC_YES.sh

#------------------------------------------------------------------------------#
################################################################################
export SCAFESRUN_N_PROCESSES_MPI_START=1
export SCAFESRUN_N_THREADS_OPENMP_START=8
export SCAFESRUN_N_NODES="512x512x512"
export SCAFESRUN_N_PARTITIONS="1x1x1"
date
scafesrun.sh
date

