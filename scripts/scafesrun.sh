#!/bin/bash
#
# Call: ./scafesrun.sh
# Shell script for running a ScaFES application which was configured via
# environment variables prefixed with SCAFESRUN.
#
# * NECESSARY: Name of ScaFES executable. <string>
#   SCAFESRUN_NAME_EXECUTABLE="$1";
#
# * NECESSARY: Space dimension of problem. <int>
#   SCAFESRUN_SPACE_DIM="3";
#
# * Number of grid nodes in each direction. n1xn2x...xnd <string>
#   SCAFESRUN_N_NODES="8x8x8";
#
# * Direction into which the grid can be partitioned when using the RCB
#   algorithm for the domain decomposition. d1xd2x...xdd <string>
#   SCAFESRUN_PARTITION_GRID="0x0x1";
#
# * Coordinates of left lower point of the global computational domain.
#   s1xs2x...xsd with s_i \in R <string>
#   SCAFESRUN_COORD_NODE_FIRST="0x0x0";
#
# * Coordinates of right upper point of the global computational domain.
#   e1xe2x...xed with e_i \in R <string>
#   SCAFESRUN_COORD_NODE_LAST="1x1x1";
#
# * Number of timesteps. <int>
#   SCAFESRUN_N_TIMESTEPS="1";
#
# * Number of snapshots when writing data files. <int>
#   SCAFESRUN_N_SNAPSHOTS="1";
#
# * Start time t_S of time interval [t_S; t_E]. <double>
#   SCAFESRUN_START_TIME="0";
#
# * End time t_E of time interval [t_S; t_E]. <double>
#   SCAFESRUN_END_TIME="1";
#
# * Output level <int>
#   SCAFESRUN_OUTPUT_LEVEL="5";
#
# * Number of layers at the global boundary <int>
#   SCAFESRUN_NLAYERS_AT_BORDER="1";
#
# * Name of kind file if boundary description(=kind field) should be read in.
#   <string>
#   SCAFESRUN_NAME_KINDFILE="";
#
# * Should the kind field be written to a file? <NO/YES>
#   SCAFESRUN_WRITE_KINDFILE="NO";
#
#   Name of the configuration file (ini format) <string>
#   SCAFESRUN_NAME_CONFIGFILE="";
#
# * Enable ADOL-C for evaluating function (and eventually computing
#   derivatives)? <NO/YES>
#   SCAFESRUN_ENABLE_ADOLC="NO";
#
# * Write partition file resulting from the MPI domoain decomposition? <NO/YES>
#   SCAFESRUN_WRITE_PARTITIONFILE="NO";
#
# * Enable asynchronous MPI communication for synchronisation of fields?
#   <NO/YES>
#   SCAFESRUN_ASYNCHRON_MODE="YES";
#
# * Use Boost.MPI skeleton concept for the MPI communication? <NO/YES>
#   SCAFESRUN_USE_BOOST_MPI_SKELETON_CONCEPT="YES";
#
# * Compute gradients using ADOL-C? <NO/YES>
#   SCAFESRUN_COMP_GRADIENTS="NO"
#
# * Name of directory in which results should be written <string>.
#   SCAFESRUN_RESULTS_DIR="./";
#
# * Check results of ScaFES application? <NO/YES>
#   SCAFESRUN_CHECK_RESULTS="NO";
#
# * Run mode of scafesrun.sh script?
#   DEBUG, DEBUG_MPI, NORMAL, LSF, SLURM <string>
#   SCAFESRUN_RUN_MODE="DEBUG";
#
# * Number of runs with different number of MPI processes. <int>
#   SCAFESRUN_N_TESTS_MPI="1";
#
# * Number of runs with different number of OpenMP threads. <int>
#   SCAFESRUN_N_TESTS_OPENMP="1";
#
# * Number of MPI processes during the first run. <int>
#   SCAFESRUN_N_PROCESSES_MPI_START="1";
#
# * Number of OpenMP threads during the first run. <int>
#   SCAFESRUN_N_THREADS_OPENMP_START="1";
#
# * Type of scaling test <string> <STRONG_MPI/WEAK_MPI/STRONG_OPENMP/WEAK_OPENMP>
#   SCAFESRUN_TYPE_SCALINGTEST="STRONG_MPI";
#
# * Name of machine on which application will be run <string>
#   SCAFESRUN_MACHINE_NAME="LOCAL";
#
# * Number of nodes of the machine. <int>
#   SCAFESRUN_MACHINE_N_NODES="1";
#
# * Number of cores per node of the machine. <int>
#   SCAFESRUN_MACHINE_N_CORES_PER_NODE="1";
#
# * Name of the project of the machine. <string>
#   SCAFESRUN_MACHINE_NAME_PROJECT="";
#
# * Name of reservation on the machine. <string>
#   SCAFESRUN_MACHINE_NAME_RESERVATION="";
#
# * Name of the partition of the machine. <string>
#   SCAFESRUN_MACHINE_NAME_PARTITION="";
#
# * Desired memory per core at the machine in MB. <int>
#   SCAFESRUN_MACHINE_MEMORY_PER_CORE="900M";
#
# * Maximal run time in hh:mm:ss. <string>
#   SCAFESRUN_MACHINE_TIME="";
#
# * Name of the init file (netCDF file). <string>
#   SCAFESRUN_NAME_INITFILE="";
#
# * Threshold for convergence check. <double>
#   SCAFESRUN_THRESHOLD="";
#
# * Number of iteration when convergence should be checked
#   for the first time. <int>
#   SCAFESRUN_CHECK_CONV_FIRST_AT_ITER="";
#
# * Number of iterations between two convergence checks. <int>
#   SCAFESRUN_CHECK_CONV_AT_EVERY_N_ITER="";
#
# ScaFES
# Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
# For details, see the files COPYING and LICENSE in the base directory
# of the package.

(

#------------------------------------------------------------------------------#
# Variables which must be set.
#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_NAME_EXECUTABLE" == "x" ] ; then
    echo "* ERROR: SCAFESRUN_NAME_EXECUTABLE was not set."
    exit 22;
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_SPACE_DIM" == "x" ] ; then
    echo "* ERROR: SCAFESRUN_SPACE_DIM was not set."
    exit 22;
fi

#------------------------------------------------------------------------------#
declare    local nameProblem=${SCAFESRUN_NAME_EXECUTABLE};
declare -i local spaceDim=$SCAFESRUN_SPACE_DIM;

#------------------------------------------------------------------------------#
# Variables which can be set.
#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_N_TESTS_MPI" == "x" ] ; then
    SCAFESRUN_N_TESTS_MPI="1";
    echo "* WARNING: Use default value SCAFESRUN_N_TESTS_MPI=$SCAFESRUN_N_TESTS_MPI"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_N_TESTS_OPENMP" == "x" ] ; then
    SCAFESRUN_N_TESTS_OPENMP="1";
    echo "* WARNING: Use default value SCAFESRUN_N_TESTS_OPENMP=$SCAFESRUN_N_TESTS_OPENMP"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_TYPE_SCALINGTEST" == "x" ] ; then
    SCAFESRUN_TYPE_SCALINGTEST="STRONG_MPI";
    echo "* WARNING: Use default value SCAFESRUN_TYPE_SCALINGTEST=$SCAFESRUN_TYPE_SCALINGTEST"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_N_PROCESSES_MPI_START" == "x" ] ; then
    SCAFESRUN_N_PROCESSES_MPI_START="1";
    echo "* WARNING: Use default value SCAFESRUN_N_PROCESSES_MPI_START=$SCAFESRUN_N_PROCESSES_MPI_START"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_N_THREADS_OPENMP_START" == "x" ] ; then
    SCAFESRUN_N_THREADS_OPENMP_START="1";
    echo "* WARNING: Use default value SCAFESRUN_N_THREADS_OPENMP_START=$SCAFESRUN_N_THREADS_OPENMP_START"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_N_TIMESTEPS" == "x" ] ; then
    SCAFESRUN_N_TIMESTEPS="1";
    echo "* WARNING: Use default value SCAFESRUN_N_TESTS_MPI=$SCAFESRUN_N_TIMESTEPS"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_N_SNAPSHOTS" == "x" ] ; then
    SCAFESRUN_N_SNAPSHOTS="1";
    echo "* WARNING: Use default value SCAFESRUN_N_SNAPSHOTS=$SCAFESRUN_N_SNAPSHOTS"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_START_TIME" == "x" ] ; then
    SCAFESRUN_START_TIME="0";
    echo "* WARNING: Use default value SCAFESRUN_START_TIME=$SCAFESRUN_START_TIME"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_END_TIME" == "x" ] ; then
    SCAFESRUN_END_TIME="1";
    echo "* WARNING: Use default value SCAFESRUN_END_TIME=$SCAFESRUN_END_TIME"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_NAME_INITFILE" == "x" ] ; then
    SCAFESRUN_NAME_INITFILE=""
    echo "* WARNING: Use default value SCAFESRUN_NAME_INITFILE=$SCAFESRUN_NAME_INITFILE"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_NAME_CONFIGFILE" == "x" ] ; then
    SCAFESRUN_NAME_CONFIGFILE="";
    echo "* WARNING: Use default value SCAFESRUN_NAME_CONFIGFILE=$SCAFESRUN_NAME_CONFIGFILE"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_NAME_KINDFILE" == "x" ] ; then
    SCAFESRUN_NAME_KINDFILE="";
    echo "* WARNING: Use default value SCAFESRUN_NAME_KINDFILE=$SCAFESRUN_NAME_KINDFILE"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_WRITE_KINDFILE" == "x" ] ; then
    SCAFESRUN_WRITE_KINDFILE="NO";
    echo "* WARNING: Use default value SCAFESRUN_WRITE_KINDFILE=$SCAFESRUN_WRITE_KINDFILE"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_WRITE_PARTITIONFILE" == "x" ] ; then
    SCAFESRUN_WRITE_PARTITIONFILE="NO";
    echo "* WARNING: Use default value SCAFESRUN_WRITE_PARTITIONFILE=$SCAFESRUN_WRITE_PARTITIONFILE"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_ENABLE_ADOLC" == "x" ] ; then
    SCAFESRUN_ENABLE_ADOLC="NO";
    echo "* WARNING: Use default value SCAFESRUN_ENABLE_ADOLC=$SCAFESRUN_ENABLE_ADOLC"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_ASYNCHRON_MODE" == "x" ] ; then
    SCAFESRUN_ASYNCHRON_MODE="YES";
    echo "* WARNING: Use default value SCAFESRUN_ASYNCHRON_MODE=$SCAFESRUN_ASYNCHRON_MODE"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_USE_BOOST_MPI_SKELETON_CONCEPT" == "x" ] ; then
    SCAFESRUN_USE_BOOST_MPI_SKELETON_CONCEPT="YES";
    echo "* WARNING: Use default value SCAFESRUN_USE_BOOST_MPI_SKELETON_CONCEPT=$SCAFESRUN_USE_BOOST_MPI_SKELETON_CONCEPT"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_COMP_GRADIENTS" == "x" ] ; then
    SCAFESRUN_COMP_GRADIENTS="NO";
    echo "* WARNING: Use default value SCAFESRUN_COMP_GRADIENTS=$SCAFESRUN_COMP_GRADIENTS"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_RUN_MODE" == "x" ] ; then
    SCAFESRUN_RUN_MODE="NORMAL";
    echo "* WARNING: Use default value SCAFESRUN_RUN_MODE=$SCAFESRUN_RUN_MODE"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_OUTPUT_LEVEL" == "x" ] ; then
    SCAFESRUN_OUTPUT_LEVEL="3";
    echo "* WARNING: Use default value SCAFESRUN_OUTPUT_LEVEL=$SCAFESRUN_OUTPUT_LEVEL"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_MACHINE_N_NODES" == "x" ] ; then
    SCAFESRUN_MACHINE_N_NODES="1";
    echo "* WARNING: Use default value SCAFESRUN_MACHINE_N_NODES=$SCAFESRUN_MACHINE_N_NODES"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_MACHINE_NAME" == "x" ] ; then
    SCAFESRUN_MACHINE_NAME="LOCAL";
    echo "* WARNING: Use default value SCAFESRUN_MACHINE_NAME=$SCAFESRUN_MACHINE_NAME"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_MACHINE_N_CORES_PER_NODE" == "x" ] ; then
    SCAFESRUN_MACHINE_N_CORES_PER_NODE="1";
    echo "* WARNING: Use default value SCAFESRUN_MACHINE_N_CORES_PER_NODE=$SCAFESRUN_MACHINE_N_CORES_PER_NODE"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_NLAYERS_AT_BORDER" == "x" ] ; then
    SCAFESRUN_NLAYERS_AT_BORDER="1";
    echo "* WARNING: Use default value SCAFESRUN_NLAYERS_AT_BORDER=$SCAFESRUN_NLAYERS_AT_BORDER"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_RESULTS_DIR" == "x" ] ; then
    SCAFESRUN_RESULTS_DIR="$PWD";
    echo "* WARNING: Use default value SCAFESRUN_RESULTS_DIR=$SCAFESRUN_RESULTS_DIR"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_CHECK_RESULTS" == "x" ] ; then
    SCAFESRUN_CHECK_RESULTS="NO";
    echo "* WARNING: Use default value SCAFESRUN_CHECK_RESULTS=$SCAFESRUN_CHECK_RESULTS"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_N_PARTITIONS" == "x" ] ; then
    SCAFESRUN_N_PARTITIONS="$SCAFESRUN_N_PROCESSES_MPI_START";
    echo "* WARNING: Use default value SCAFESRUN_N_PARTITIONS=$SCAFESRUN_N_PARTITIONS"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_N_NODES" == "x" ] ; then
    SCAFESRUN_N_NODES="8";
    for idxDim in `seq 2 $spaceDim`; do
        tmpNnodes="${SCAFESRUN_N_NODES}x8";
        SCAFESRUN_N_NODES=$tmpNnodes;
    done
    echo "* WARNING: Use default value SCAFESRUN_N_NODES=$SCAFESRUN_N_NODES"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_COORD_NODE_FIRST" == "x" ] ; then
    SCAFESRUN_COORD_NODE_FIRST="0";
    for idxDim in `seq 2 $spaceDim`; do
        tmpCf="${SCAFESRUN_COORD_NODE_FIRST}x0";
        SCAFESRUN_COORD_NODE_FIRST=$tmpCf;
    done
    echo "* WARNING: Use default value SCAFESRUN_COORD_NODE_FIRST=$SCAFESRUN_COORD_NODE_FIRST"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_COORD_NODE_LAST" == "x" ] ; then
    SCAFESRUN_COORD_NODE_LAST="1";
    for idxDim in `seq 2 $spaceDim`; do
        tmpCl="${SCAFESRUN_COORD_NODE_LAST}x1";
        SCAFESRUN_COORD_NODE_LAST=$tmpCl;
    done
    echo "* WARNING: Use default value SCAFESRUN_COORD_NODE_LAST=$SCAFESRUN_COORD_NODE_LAST"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_TYPE_DOMAIN_DECOMPOSITION" == "x" ] ; then
    SCAFESRUN_TYPE_DOMAIN_DECOMPOSITION="RCB";
    echo "* WARNING: Use default value SCAFESRUN_TYPE_DOMAIN_DECOMPOSITION=$SCAFESRUN_TYPE_DOMAIN_DECOMPOSITION"
fi

#------------------------------------------------------------------------------#
# Only relevant if SCAFESRUN_TYPE_DOMAIN_DECOMPOSITION="RCB"; is set.
if [ "x$SCAFESRUN_PARTITION_GRID" == "x" ] ; then
    SCAFESRUN_PARTITION_GRID="1";
    if (("$spaceDim" > 1)); then
        for idxDim in `seq 2 $spaceDim`; do
            tmpPg="${SCAFESRUN_PARTITION_GRID}x1";
            SCAFESRUN_PARTITION_GRID=$tmpPg;
        done
    fi
    echo "* WARNING: Use default value SCAFESRUN_PARTITION_GRID=$SCAFESRUN_PARTITION_GRID"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_MACHINE_MEMORY_PER_CORE" == "x" ] ; then
    if [ "x$SCAFESRUN_RUN_MODE" == "xLSF" -o "$SCAFESRUN_RUN_MODE" == "xSLURM" ] ; then
        echo "* ERROR: SCAFESRUN_MACHINE_MEMORY_PER_CORE was not set."
        exit 22;
    else
        SCAFESRUN_MACHINE_MEMORY_PER_CORE="900";
        echo "* WARNING: Use default value SCAFESRUN_MACHINE_MEMORY_PER_CORE=$SCAFESRUN_MACHINE_MEMORY_PER_CORE";
    fi
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_MACHINE_NAME_PROJECT" == "x" ] ; then
    if [ "x$SCAFESRUN_RUN_MODE" == "xLSF" -o "$SCAFESRUN_RUN_MODE" == "xSLURM" ] ; then
        echo "* ERROR: SCAFESRUN_MACHINE_NAME_PROJECT was not set."
        exit 22;
    else
        SCAFESRUN_MACHINE_NAME_PROJECT="";
        echo "* WARNING: Use default value SCAFESRUN_MACHINE_NAME_PROJECT=$SCAFESRUN_MACHINE_NAME_PROJECT"
    fi
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_MACHINE_NAME_RESERVATION" == "x" ] ; then
    if [ "x$SCAFESRUN_RUN_MODE" == "xLSF" -o "$SCAFESRUN_RUN_MODE" == "xSLURM" ] ; then
        echo "* ERROR: SCAFESRUN_MACHINE_NAME_RESERVATION was not set."
        exit 22;
    else
        SCAFESRUN_MACHINE_NAME_RESERVATION="";
        echo "* WARNING: Use default value SCAFESRUN_MACHINE_NAME_RESERVATION=$SCAFESRUN_MACHINE_NAME_RESERVATION"
    fi
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_MACHINE_NAME_PARTITION" == "x" ] ; then
    if [ "x$SCAFESRUN_RUN_MODE" == "xLSF" -o "$SCAFESRUN_RUN_MODE" == "xSLURM" ] ; then
        echo "* ERROR: SCAFESRUN_MACHINE_NAME_PARTITION was not set."
        exit 22;
    fi
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_MACHINE_TIME" == "x" ] ; then
    if [ "x$SCAFESRUN_RUN_MODE" == "xLSF" -o "$SCAFESRUN_RUN_MODE" == "xSLURM" ] ; then
        echo "* ERROR: SCAFESRUN_MACHINE_TIME was not set."
        exit 22;
    fi
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_THRESHOLD" == "x" ] ; then
    SCAFESRUN_THRESHOLD="0.00001";
    echo "* WARNING: Use default value SCAFESRUN_THRESHOLD=$SCAFESRUN_THRESHOLD"
fi
#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_CHECK_CONV_FIRST_AT_ITER" == "x" ] ; then
    SCAFESRUN_CHECK_CONV_FIRST_AT_ITER="1";
    echo "* WARNING: Use default value SCAFESRUN_CHECK_CONV_FIRST_AT_ITER=$SCAFESRUN_CHECK_CONV_FIRST_AT_ITER"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_CHECK_CONV_AT_EVERY_N_ITER" == "x" ] ; then
    SCAFESRUN_CHECK_CONV_AT_EVERY_N_ITER="1";
    echo "* WARNING: Use default value SCAFESRUN_CHECK_CONV_AT_EVERY_N_ITER=$SCAFESRUN_CHECK_CONV_AT_EVERY_N_ITER"
fi

#------------------------------------------------------------------------------#
if [ "x$SCAFESRUN_SOFTWARE_VERSION" == "x" ] ; then
    SCAFESRUN_SOFTWARE_VERSION="trunk";
    echo "* WARNING: Use default value SCAFESRUN_SOFTWARE_VERSION=$SCAFESRUN_SOFTWARE_VERSION"
fi

#------------------------------------------------------------------------------#
# From here on, all variables should be declared.
set -u
set -e

#------------------------------------------------------------------------------#
# All SCAFESRUN* environment variables were exported as strings.
# Thus, those SCAFESRUN* variables which should be used as integer variables
# have to be declared explicitly as local integers variables.
declare -i local nTestsMpi=$SCAFESRUN_N_TESTS_MPI;
declare -i -r local nTestsOpenMp=$SCAFESRUN_N_TESTS_OPENMP;
declare    -r local testForWeakScaling=$SCAFESRUN_TYPE_SCALINGTEST;
declare -i local nProcessesMpiStart=$SCAFESRUN_N_PROCESSES_MPI_START;
declare -i -r local nThreadsOpenMpStart=$SCAFESRUN_N_THREADS_OPENMP_START;
declare -i -r local nTimesteps=$SCAFESRUN_N_TIMESTEPS;
declare -i -r local nSnapshots=$SCAFESRUN_N_SNAPSHOTS;
declare -i -r local starttime=$SCAFESRUN_START_TIME;
declare -i local endtime=$SCAFESRUN_END_TIME; #TODO: Must be set by hand.
declare    -r local nameConfigfile=$SCAFESRUN_NAME_CONFIGFILE;
declare    -r local nameKindfile=$SCAFESRUN_NAME_KINDFILE;
declare    -r local writeKindfile=$SCAFESRUN_WRITE_KINDFILE;
declare    -r local enableAdolc=$SCAFESRUN_ENABLE_ADOLC;
declare    -r local asynchronMode=$SCAFESRUN_ASYNCHRON_MODE;
declare    -r local useSkeletonConcept=$SCAFESRUN_USE_BOOST_MPI_SKELETON_CONCEPT;
declare    -r local computeGradients=$SCAFESRUN_COMP_GRADIENTS;
declare    -r local runMode=$SCAFESRUN_RUN_MODE;
declare    -r local outputlevel=$SCAFESRUN_OUTPUT_LEVEL;
declare -i -r local oneAsInteger=1;
declare -i local nNodesMachineMAX=$SCAFESRUN_MACHINE_N_NODES;
declare -i local nCoresPerNodeMAX=$SCAFESRUN_MACHINE_N_CORES_PER_NODE;
declare    local nameMachine=${SCAFESRUN_MACHINE_NAME};
declare    -r local softwareVersion=${SCAFESRUN_SOFTWARE_VERSION};
declare    -r local nLayersAtBorder="${SCAFESRUN_NLAYERS_AT_BORDER}";
declare    -r local resDir=${SCAFESRUN_RESULTS_DIR};
declare    -r local checkResults=${SCAFESRUN_CHECK_RESULTS};
declare    local nPartitions=${SCAFESRUN_N_PARTITIONS};
declare    -r local typeDomainDecomposition=${SCAFESRUN_TYPE_DOMAIN_DECOMPOSITION};
declare    -r local nameProject=${SCAFESRUN_MACHINE_NAME_PROJECT};
declare    -r local nameReservation=${SCAFESRUN_MACHINE_NAME_RESERVATION};
declare    -r local memoryPerCore=${SCAFESRUN_MACHINE_MEMORY_PER_CORE};
declare    local jobexestring="";
declare    -r local nameInitfile=${SCAFESRUN_NAME_INITFILE};
declare    -r local valueThreshold=${SCAFESRUN_THRESHOLD};
declare -i -r local valueCheckConvFirstAtIter=${SCAFESRUN_CHECK_CONV_FIRST_AT_ITER};
declare -i -r local valueCheckConvAtEveryNIter=${SCAFESRUN_CHECK_CONV_AT_EVERY_N_ITER};
#------------------------------------------------------------------------------#
# Compute value of dependent variable.
declare -i local nCoresTotalMAX=${nNodesMachineMAX};
nCoresTotalMAX=$((nCoresTotalMAX*nCoresPerNodeMAX));

#------------------------------------------------------------------------------#
if [ "x$runMode" == "xDEBUG" ] ; then
    nTestsMpi=1;
    nProcessesMpiStart=1;
    #nThreadsOpenMpStart=1;
fi
#------------------------------------------------------------------------------#
if [ "x$nameKindfile" == "x" ] ; then
    optionKindfile="";
else
    optionKindfile="--kindfile=$nameKindfile";
fi
#------------------------------------------------------------------------------#
if [ "x$nameConfigfile" == "x" ] ; then
    optionConfigfile="";
else
    optionConfigfile="--configfile=$nameConfigfile";
fi
#------------------------------------------------------------------------------#
if [ "x$nameInitfile" == "x" ] ; then
    optionInitfile="";
else
    optionInitfile="--initfile=$nameInitfile";
fi
#------------------------------------------------------------------------------#
boolEnabledAdolc=0;
if [ "x$enableAdolc" == "xYES" ] ; then
    boolEnabledAdolc=1;
fi
boolAsynchronMode=0;
if [ "x$asynchronMode" == "xYES" ] ; then
    boolAsynchronMode=1;
fi
boolUseSkeletonConcept=0;
if [ "x$useSkeletonConcept" == "xYES" ] ; then
    boolUseSkeletonConcept=1;
fi
boolCompGradients=0;
if [ "x$computeGradients" == "xYES" ] ; then
    boolCompGradients=1;
fi
boolWriteKindfile=0;
if [ "x$writeKindfile" == "xYES" ] ; then
    boolWriteKindfile=1;
fi

declare -i local endTestMpi;
declare -i local endTestOpenMp;
declare    local currNameProblem;
declare -i local currIdxProblem;
declare -i local currTest;
declare -i local currIdTest;
declare -i local currIdTestMpi;
declare -i local currIdTestOpenMp;
declare -i local currNprocessesMpi;
declare -i local currNthreadsOpenMp;
declare -i local currEndTest;
declare    local currNameExecutable; # Relative path to executable
declare -i tmpNodesDir;
declare    local errCode=0;
endTestMpi=${nTestsMpi}-${oneAsInteger};
endTestOpenMp=${nTestsOpenMp}-${oneAsInteger};
currEndTest=${nTestsMpi}*${nTestsOpenMp};

timestamp=$(date +%Y%m%d%H%M);
declare local nNodesMachine=1;

currNameProblem=${nameProblem};
currNameExecutable="./${currNameProblem}";
currIdTest=$oneAsInteger;
(
#------------------------------------------------------------------------------#
# Perform a loop over all MPI tests.
currNprocessesMpi=${nProcessesMpiStart};
for idxTestMpi in `seq 0 $endTestMpi`; do
        currIdTestMpi=${idxTestMpi}+${oneAsInteger};

        #----------------------------------------------------------------------#
        # Perform a loop over all OpenMP tests.
        currNthreadsOpenMp=${nThreadsOpenMpStart};
        for idxTestOpenMp in `seq 0 $endTestOpenMp`; do
            currIdTestOpenMp=${idxTestOpenMp}+${oneAsInteger};
            #currNnodes=(${SCAFESRUN_N_NODES[*]})
            currStarttime=${starttime};
            currEndtime=${endtime};
            currNtimesteps=${nTimesteps};
            currNsnapshots=${nSnapshots};

            #------------------------------------------------------------------#
            divideGrid=${SCAFESRUN_PARTITION_GRID};
            currNnodes="${SCAFESRUN_N_NODES}";
            currCoordNodeFirst="${SCAFESRUN_COORD_NODE_FIRST}";
            currCoordNodeLast="${SCAFESRUN_COORD_NODE_LAST}";
            namePartfile="${currNameProblem}_${currNprocessesMpi}";
            namePartfile="${namePartfile}_${currNthreadsOpenMp}_${currNnodes}";
            nameDatafile=${namePartfile};

            #------------------------------------------------------------------#
            str_currNprocessesMpi=$(printf "%.5d" $currNprocessesMpi);
            str_currNthreadsOpenMp=$(printf "%.5d" $currNthreadsOpenMp);
            str_nTimesteps=$(printf "%.6d" $currNtimesteps);
            basefile="${resDir}/${currNameProblem}"
            basefile="$basefile-${nameMachine}"
            basefile="$basefile-${softwareVersion}"
            basefile="$basefile-mpi${str_currNprocessesMpi}"
            basefile="$basefile-omp${str_currNthreadsOpenMp}"
            basefile="$basefile-nNodes${currNnodes}"
#            basefile="$basefile-decomp${typeDomainDecomposition}"
            basefile="$basefile-nPart${nPartitions}"
#            basefile="$basefile-divide${divideGrid}"
             basefile="$basefile-adolc${enableAdolc}"
             basefile="$basefile-async${asynchronMode}"
             basefile="$basefile-skelet${useSkeletonConcept}"
#             basefile="$basefile-grad${computeGradients}"
#             basefile="$basefile-nTs${str_nTimesteps}"
#             basefile="$basefile-tEnd${endtime}"
#             basefile="$basefile-nLayers${nLayersAtBorder}"
            basefile="$basefile-${timestamp}"
            outputfile="$basefile.out"
            errorfile="$basefile.err"

            #------------------------------------------------------------------#
            if ( expr ${currNprocessesMpi} \< ${nNodesMachineMAX} > /dev/null ) ; then
                nNodesMachine=1;
            else
                nNodesMachine=$((currNprocessesMpi/nNodesMachineMAX));
            fi

            declare -i local nCoresTotal=$((currNprocessesMpi*currNthreadsOpenMp));
            if [ ${nCoresTotal} -gt ${nCoresTotalMAX} ]; then
                echo "* WARNING: nCoresTotal > nCoresTotalMAX. Continue loop."
                continue;
            fi

            #------------------------------------------------------------------#
            # Set default value if nPartitions was not given at the command
            # line: nPartitions=nProcsMpix1x...x1
            # default value: nProcsMpix1x...x1
            declare local isGivenInPartitionGridFormat="no";
            case "$nPartitions" in
               *x*) isGivenInPartitionGridFormat="yes" ;;
            esac

            if [ "$isGivenInPartitionGridFormat" == "no" ] ; then
                nPartitions="$currNprocessesMpi";
                if (("$spaceDim" > 1)); then
                    for idxDim in `seq 2 $spaceDim`; do
                        tmpPart="${nPartitions}x1"
                        nPartitions=$tmpPart;
                    done
                fi
            fi

            echo -e "\n\n  +-------------------------------+";
            echo -e "  | TEST: ${currIdTest} / ${currEndTest}";
            printf "  |     * Name of Problem     = %s\n" ${currNameProblem};
            printf "  |     * #{MPI processes}    = %s\n" ${currNprocessesMpi};
            printf "  |     * #{OpenMP threads}   = %s\n" ${currNthreadsOpenMp};
            printf "  |     * #{Nodes}            = %s\n" ${currNnodes};
            printf "  |     * #{Partitions}       = %s\n" ${nPartitions};
            printf "  |     * Type domain decomp. = %s\n" ${typeDomainDecomposition};
            printf "  |     * divideGrid          = %s\n" ${divideGrid};
            printf "  |     * enableAdolc         = %s\n" ${enableAdolc};
            printf "  |     * asynchronMode       = %s\n" ${asynchronMode};
            printf "  |     * useSkeletonConcept  = %s\n" ${useSkeletonConcept};
            printf "  |     * computeGradients    = %s\n" ${computeGradients};
            printf "  |     * #{Start time}       = %s\n" ${currStarttime};
            printf "  |     * #{End time}         = %s\n" ${currEndtime};
            printf "  |     * #{Time steps}       = %s\n" ${currNtimesteps};
            printf "  |     * #{Layers At Border} = %s\n" ${nLayersAtBorder};
            printf "  |     * #{Snap shots}       = %s\n  |\n" ${currNsnapshots};
            printf "  |     * #{NodesMachineUsed} = %s\n" ${nNodesMachine};
            printf "  |     * #{CoresTotalUsed}   = %s\n  |\n" \
                        $((currNprocessesMpi*currNthreadsOpenMp));

            printf "  |     * Name of Machine     = %s\n" ${nameMachine};
            printf "  |     * #{NodesMachineMAX}  = %s\n" ${nNodesMachineMAX};
            printf "  |     * #{CoresPerNodeMAX}  = %s\n" ${nCoresPerNodeMAX};
            printf "  |     * #{CoresTotalMAX}    = %s\n" ${nCoresTotalMAX};
            echo -e "  +-------------------------------+\n\n";

            #------------------------------------------------------------------#
            echo "--- Run ScaFES application. ---"
            optsString=" \
                    --dim=$spaceDim \
                    --typeDomainDecomposition=$typeDomainDecomposition \
                    --divideGrid=$divideGrid \
                    --nPartitions=$nPartitions \
                    --nNodes=$currNnodes \
                    --coordNodeFirst=$currCoordNodeFirst \
                    --coordNodeLast=$currCoordNodeLast \
                    --starttime=$currStarttime \
                    --endtime=$currEndtime \
                    --nTimesteps=$currNtimesteps \
                    --threshold=$valueThreshold \
                    --checkConvFirstAtIter=$valueCheckConvFirstAtIter \
                    --checkConvAtEveryNIter=$valueCheckConvAtEveryNIter \
                    --nSnapshots=$currNsnapshots \
                    ${optionKindfile} \
                    --writeKindFile=$boolWriteKindfile \
                    ${optionConfigfile} \
                    ${optionInitfile} \
                    --enabledAdolc=$boolEnabledAdolc \
                    --asynchronMode=$boolAsynchronMode \
                    --useBoostMpiSkeletonConcept=$boolUseSkeletonConcept \
                    --computeGradients=$boolCompGradients \
                    --outputfile=$nameDatafile \
                    --outputlevel=$outputlevel \
                    --nLayersAtBorder=$nLayersAtBorder"

            if [ "x$SCAFESRUN_WRITE_PARTITIONFILE" == "xYES" ] ; then
                optsString="$optsString --writePartitionfile=$namePartfile";
            fi
            exeString="$currNameExecutable $optsString"

            case $runMode in
            #------------------------------------------------------------------#
            DEBUG )
                export OMP_NUM_THREADS=$currNthreadsOpenMp;
                jobexestring="libtool --mode=execute gdb --args $exeString"
                ;;
            #------------------------------------------------------------------#
            DEBUG_MPI )
                export OMP_NUM_THREADS=$currNthreadsOpenMp;
                jobexestring="libtool --mode=execute mpirun --mca orte_base_helpggregate 0 \
                        -np $currNprocessesMpi \
                        xterm -e gdb --args $exeString"
                ;;
            #------------------------------------------------------------------#
            SLURM_WITH_SBATCH_CALL )
                export OMP_NUM_THREADS=$currNthreadsOpenMp;
                jobexestring="srun -v -v $exeString"
                #jobexestring="srun -v -v \
                #     --nodes=${nNodesMachine} \
                #     --ntasks=$currNprocessesMpi \
                #     --time=${SCAFESRUN_MACHINE_TIME} \
                #     --cpus-per-task=$currNthreadsOpenMp \
                #     --cpu-freq=High \
                #     --partition=${SCAFESRUN_MACHINE_NAME_PARTITION} \
                #    $exeString"
                ;;
            #------------------------------------------------------------------#
            # TODO: Remote debugging
            # http://www.redbooks.ibm.com/redbooks/pdfs/sg247948.pdf
            # [KF, 04/22/2015]
            # runjob  --block R00-M0-N01 --cwd `pwd` \
            # --start-tool /sbin/gdbtool --exp_env OMP_NUM_THREADS \
            # --ranks-per-node=2 \
            # --tool-args "--rank=2 --listen_port=10001"
            #
            # --ranks-per-node=$currNprocessesMpi
            # TODO: Currently, this holds for Juqueen, only.
            # [KF, 04/27/2015]
            LOADLEVELER_WITH_LLSUBMIT_CALL )

               runjob -p 4 \
                 --verbose 4 : /bgsys/local/samples/personality/personality.elf > partinfos_NMPI_${currNprocessesMpi}_NNODES_${currNnodes}_${timestamp}.txt

                # Choose 16 cores per compute node for pure MPI jobs at Juqueen!
                # [KF, 06/01/2015]
                # Juqueen requires: currNthreadsOpenMp * currNprocessesMpi = 64
                # Each compute node has 4 hardware hyperthreads.
                ranksperNode=$SCAFESRUN_MACHINE_N_CORES_PER_NODE;
                #if [ $currNthreadsOpenMp -gt 1 ] ; then
                #    #ranksperNode=$(($SCAFESRUN_MACHINE_N_CORES_PER_NODE/currNthreadsOpenMp));
                #fi
                export OMP_NUM_THREADS=$currNthreadsOpenMp;
                export BG_COREDUMPDISABLED=1;
                jobexestring="runjob \
                     --exp_env OMP_NUM_THREADS \
                     --exp_env BG_COREDUMPDISABLED \
                     --ranks-per-node=$ranksperNode \
                     --np=$currNprocessesMpi \
                     : ./$exeString"
                ;;
            #------------------------------------------------------------------#
            OPENMP_ONLY )
                export OMP_NUM_THREADS=$currNthreadsOpenMp;
                jobexestring="./$exeString"
                ;;
            #------------------------------------------------------------------#
            SERIAL )
                export OMP_NUM_THREADS=1; # Just to be sure.
                jobexestring="./$exeString"
                ;;
            #------------------------------------------------------------------#
            NORMAL )
                export OMP_NUM_THREADS=$currNthreadsOpenMp;
                jobexestring="mpirun -np $currNprocessesMpi $exeString"
                ;;
            #------------------------------------------------------------------#
            NORMAL_PERF_REPORT )
                 export OMP_NUM_THREADS=$currNthreadsOpenMp;
                 declare local isInstalledPerfreport=$(which perf-report);
                 if [ "$isInstalledPerfreport" == "" ] ; then
                     echo "perf-report not available / not in PATH".
                     exit 1;
                 fi
                 jobexestring="perf-report \
                    mpirun -np $currNprocessesMpi \
                    $exeString"
                ;;
            #------------------------------------------------------------------#
            SLURM_PERF_REPORT )
                 export DDT_NO_TIMEOUT=1;
                 export OMP_NUM_THREADS=$currNthreadsOpenMp;
                 # Attention: perf-report:
                 # For 128 MPI processes:
                 # Your current licence (xxx) does not permit as many
                 # processes as you requested.
                 declare local isInstalledPerfreport=$(which perf-report);
                 if [ "$isInstalledPerfreport" == "" ] ; then
                     echo "perf-report not available / not in PATH".
                     exit 1;
                 fi
                 jobexestring="perf-report \
                    -n $currNprocessesMpi --mpi="slurm" \
                    $exeString"
                ;;
            #------------------------------------------------------------------#
            LSF | LSF_ONE_NODE )
                # -npernode 64
                export OMP_NUM_THREADS=$currNthreadsOpenMp;
                if [ $currNthreadsOpenMp -gt 1 ] ; then
                    bsubOpts="-R span[hosts=1]";
                else
                    bsubOpts="";
                fi
                #  -U ${nameReservation}
                jobexestring="bsub $bsubOpts \
                  -x \
                  -P ${nameProject} \
                  -n ${nCoresTotal} \
                  -o ${outputfile}.lsf.out \
                  -e ${errorfile}.lsf.err \
                  -J ${currNameExecutable} \
                  -M ${memoryPerCore} \
                  mpirun -np $currNprocessesMpi \
                  $exeString"
                ;;
            #------------------------------------------------------------------#
            SLURM )
                ## --reservation=${nameReservation}
                jobexestring="srun -v -v \
                    --nodes=${nNodesMachine} \
                    --ntasks=${currNprocessesMpi} \
                    --reservation=${nameReservation} \
                    --time=${SCAFESRUN_MACHINE_TIME} \
                    -o ${outputfile}.slurm.out \
                    -e ${errorfile}.slurm.err \
                    -J ${currNameExecutable} \
                    -A ${nameProject} \
                    --mem-per-cpu=${memoryPerCore} \
                    --partition=${SCAFESRUN_MACHINE_NAME_PARTITION} \
                    $exeString"
                ;;
            #------------------------------------------------------------------#
            * )
                echo "ERROR: Mode not known." ;;
            esac
            echo "+++ OMP_NUM_THREADS=$OMP_NUM_THREADS +++"
            case $runMode in
            DEBUG* )
                echo "+++ $jobexestring +++"
                $jobexestring
                ;;
            * )
                echo "+++ $jobexestring >$outputfile 2> $errorfile +++"
                $jobexestring >$outputfile 2> $errorfile
                ;;
            esac
            errCode=$?
            echo -e "";
            if [ "${errCode}" != "0" ]; then
                echo -e "RUN TIME ERROR CODE = $errCode. Exit!";
                echo -e "RUN TIME ERROR = YES";
                break 2;
            else
                echo -e "* RUN SUCCESSFULLY?    YES";
            fi
            #------------------------------------------------------------------#
            if [ "x$checkResults" == "xYES" ]; then
                grep -h "Local Linf" ${outputfile} > tmp.dat
                correctErrCode=$?
                if [ "$correctErrCode" != "0" ]; then
                    echo "CORRECTNESS CHECK: FILE NOT FOUND. Exit!";
                    errCode=1;
                    break 2;
                fi

                #--------------------------------------------------------------#
                # Disable exit-at-error mode temporarily:
                # As the search in the grep is inverted
                # the  grep command will return with an error code != 0
                # if the results are correct.
                set +e
                grep -h -v "0.0000" tmp.dat
                correctErrCode=$?
                set -e
                # If lines with error!=0 were found.
                if [ "x$correctErrCode" == "x0" ]; then
                    echo "* ARE RESULTS CORRECT? NO";
                    errCode=1;
                    break 2;
                else
                    echo "* ARE RESULTS CORRECT? YES";
                    errCode=0;
                fi
            else
                echo "* ARE RESULTS CORRECT? UNKNOWN (DID NOT CHECK)";
            fi

            currNthreadsOpenMp=$currNthreadsOpenMp*2;
            currIdTest=$currIdTest+$oneAsInteger;
        done

        #------------------------------------------------------------------#
        # Test for weak scalability:
        # Problem size increases linearly with number of processes.
        #------------------------------------------------------------------#
        # Test for strong scalability:
        # Problem size remains constant.
        # The number of processes increases linearly.
        case $testForWeakScaling in
        WEAK* )
            echo "TODO: nNodesTotal has to be doubled!"; ;;
        * )
            currNprocessesMpi=$currNprocessesMpi*2;
            nPartitions="$currNprocessesMpi";
            for idxDim in `seq 2 $spaceDim`; do
                tmpPart="${nPartitions}x1"
                nPartitions=$tmpPart;
            done
            ;;
        esac
done
exit $errCode;
)
)
