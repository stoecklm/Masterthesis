#!/bin/bash

#SBATCH -J "HeatEqnOne_48"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=ht
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=00:02:00
#SBATCH --exclusive

# ERROR: [KF, 06/03/0215]
# sbatch: error: Batch job submission failed: Requested node configuration is not available

source ./tests/scaling/ZIH_TAURUS2_GCC_4_9_1/run-HeatEqnOne3D-timesteps_20-ADOLC_NO-MPI_ASYNC_YES.sh

#------------------------------------------------------------------------------#
################################################################################
export SCAFESRUN_N_PROCESSES_MPI_START=1
export SCAFESRUN_N_THREADS_OPENMP_START=48
export SCAFESRUN_N_NODES="256x192x256"
export SCAFESRUN_N_PARTITIONS="1x1x1"
date
scafesrun.sh
date


