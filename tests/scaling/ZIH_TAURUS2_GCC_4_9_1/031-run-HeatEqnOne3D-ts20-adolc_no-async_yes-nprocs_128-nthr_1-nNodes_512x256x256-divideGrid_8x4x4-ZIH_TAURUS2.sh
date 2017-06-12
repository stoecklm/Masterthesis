#!/bin/bash

#SBATCH -J "HeatEqnOne_128"
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --partition=haswell
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=00:02:00
#SBATCH --exclusive

source ./tests/scaling/ZIH_TAURUS2_GCC_4_9_1/run-HeatEqnOne3D-timesteps_20-ADOLC_NO-MPI_ASYNC_YES.sh

#------------------------------------------------------------------------------#
################################################################################
export SCAFESRUN_N_PROCESSES_MPI_START=128
export SCAFESRUN_N_THREADS_OPENMP_START=1
export SCAFESRUN_N_NODES="512x256x256"
export SCAFESRUN_N_PARTITIONS="8x4x4"
date
scafesrun.sh
date
