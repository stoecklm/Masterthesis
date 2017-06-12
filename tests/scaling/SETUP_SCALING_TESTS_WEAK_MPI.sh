#!/bin/bash
# Common settings for weak scaling tests wrt. MPI.
# Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
# For details, see the files COPYING and LICENSE in the base directory
# of the package.

export SCAFESRUN_TYPE_SCALINGTEST="WEAK_MPI";
export SCAFESRUN_N_THREADS_OPENMP_START=1;
export SCAFESRUN_N_TESTS_OPENMP=1;

export SCAFESRUN_ASYNCHRON_MODE="NO"
export SCAFESRUN_ENABLE_ADOLC="NO"
export SCAFESRUN_COMP_GRADIENTS="NO"
scafesscalingtest.sh

export SCAFESRUN_ASYNCHRON_MODE="YES"
export SCAFESRUN_ENABLE_ADOLC="NO"
export SCAFESRUN_COMP_GRADIENTS="NO"
scafesscalingtest.sh
