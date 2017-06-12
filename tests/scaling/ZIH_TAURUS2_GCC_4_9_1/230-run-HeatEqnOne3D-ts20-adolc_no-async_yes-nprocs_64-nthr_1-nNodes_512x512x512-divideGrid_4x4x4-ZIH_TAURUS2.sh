#!/bin/bash

#SBATCH -J "HeatEqnOne_64"
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --partition=haswell
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=00:40:00
#SBATCH --exclusive

source ./tests/scaling/ZIH_TAURUS2_GCC_4_9_1/run-HeatEqnOne3D-timesteps_20-ADOLC_NO-MPI_ASYNC_YES.sh

#------------------------------------------------------------------------------#
################################################################################
export SCAFESRUN_N_PROCESSES_MPI_START=64
export SCAFESRUN_N_THREADS_OPENMP_START=1
export SCAFESRUN_N_NODES="512x512x512"
export SCAFESRUN_N_PARTITIONS="4x4x4"
date
scafesrun.sh
date
