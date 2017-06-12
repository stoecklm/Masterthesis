#!/bin/bash
# Main script.
# Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
# For details, see the files COPYING and LICENSE in the base directory
# of the package.

(
set -u
set -e

source ./tests/scaling/SETUP_KERNEL_HEATEQNFDM.sh

source ./tests/scaling/SETUP_MACHINE_TAURUS2_THROUGHPUT.sh

export SCAFESRUN_N_THREADS_OPENMP_START=1;
export SCAFESRUN_N_TESTS_OPENMP=1;
export SCAFESRUN_N_NODES="64x64x64";
export SCAFESRUN_N_PARTITIONS="1x1x1";
source ./tests/scaling/SETUP_SCALING_TESTS_WEAK_OPENMP.sh

)
