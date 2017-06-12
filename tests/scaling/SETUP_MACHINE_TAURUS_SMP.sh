#!/bin/bash
# Common machine settings.
# Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
# For details, see the files COPYING and LICENSE in the base directory
# of the package.

export SCAFESRUN_RUN_MODE="SLURM_WITH_SBATCH_CALL";
export SCAFESRUN_MACHINE_NAME="ZIH_TAURUS_GCC_4_9_1"
export SCAFESRUN_MACHINE_NAME_PROJECT="zihforschung";
export SCAFESRUN_MACHINE_NAME_PARTITION="smp"
export SCAFESRUN_MACHINE_MEMORY_PER_CORE=30000;
export SCAFESRUN_MACHINE_N_NODES=2;
export SCAFESRUN_MACHINE_N_CORES_PER_NODE=32;
export SCAFESRUN_MACHINE_TIME="00:05:00";
