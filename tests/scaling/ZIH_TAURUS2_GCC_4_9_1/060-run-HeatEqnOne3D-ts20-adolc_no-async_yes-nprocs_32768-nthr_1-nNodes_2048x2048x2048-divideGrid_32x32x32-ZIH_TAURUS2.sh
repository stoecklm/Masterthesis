#!/bin/bash

#SBATCH -J "HeatEqnOne_32768"
#SBATCH --nodes=2048
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --partition=alle
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=00:02:00
#SBATCH --exclusive

source ./tests/scaling/ZIH_TAURUS2_GCC_4_9_1/run-HeatEqnOne3D-timesteps_20-ADOLC_NO-MPI_ASYNC_YES.sh

#------------------------------------------------------------------------------#
################################################################################
export SCAFESRUN_N_PROCESSES_MPI_START=32768
export SCAFESRUN_N_THREADS_OPENMP_START=1
export SCAFESRUN_N_NODES="2048x2048x2048"
export SCAFESRUN_N_PARTITIONS="32x32x32"
date
scafesrun.sh
date
