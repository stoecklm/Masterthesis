#!/bin/bash
# Call: ./scafesscalingtest.sh
#
# Weak scaling test wrt. MPI or OpenMP for machines based on a job
# submission system.
# The run times all of tests should be approximately the same.
#
# If nProcessesMpiStart > 1, then
# nGridNodes and nPartitions must be set via SCAFESRUN* environment variables.
#
# ScaFES
# Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
# For details, see the files COPYING and LICENSE in the base directory
# of the package.

(
set -u
set -e

################################################################################
#------------------------------------------------------------------------------#
declare local JOBSYSTEM_JOBID=
declare local JOBSYSTEM_NODES=
declare local JOBSYSTEM_TASKS_PER_NODE=
declare local JOBSYSTEM_CPUS_PER_TASK=
declare local JOBSYSTEM_NAME=
declare local JOBSYSTEM_JOBFILE=

#------------------------------------------------------------------------------#
declare -r local dateCurr="$(date +%Y_%m_%d_%H%M)"
declare -i -r local nThreadsMax=$SCAFESRUN_MACHINE_N_CORES_PER_NODE
declare -i -r local nProcessesMax=$((SCAFESRUN_MACHINE_N_NODES*SCAFESRUN_MACHINE_N_CORES_PER_NODE));
declare -r local nNodesVectFull=${SCAFESRUN_N_NODES//x/ };
declare -a local nNodesVect=( ${nNodesVectFull} );
declare -r local nPartitionsVectFull=${SCAFESRUN_N_PARTITIONS//x/ };
declare -a local nPartitionsVect=( ${nPartitionsVectFull} );

declare local nIter=0;
declare local dir=$((nIter%SCAFESRUN_SPACE_DIM));
declare local tmpNodes=0;
declare local tmpStringNodes=;
declare local tmpPartitions=0;
declare local tmpStringPartitions=;
declare -r local scalingtestFactor=2;
declare -r local nScalingTestsMpi=$SCAFESRUN_N_TESTS_MPI;
declare -r local nScalingTestsOpenmp=$SCAFESRUN_N_TESTS_OPENMP;
declare -i local nProcessesCurr=$SCAFESRUN_N_PROCESSES_MPI_START;
declare -i local nThreadsCurr=$SCAFESRUN_N_THREADS_OPENMP_START

# Dirty hack: Reset both variables to 1 because there do exist two loops
# over all MPI and over all OpenMP tests within this file.
export SCAFESRUN_N_TESTS_MPI=1;
export SCAFESRUN_N_TESTS_OPENMP=1;

################################################################################
#------------------------------------------------------------------------------#
# Perform a loop over all MPI tests.
nProcessesCurr=$SCAFESRUN_N_PROCESSES_MPI_START;
for idxTestMpi in `seq 1 $nScalingTestsMpi`; do

    #--------------------------------------------------------------------------#
    # Perform a loop over all OpenMP tests.
    nThreadsCurr=$SCAFESRUN_N_THREADS_OPENMP_START;
    for idxTestOpenMp in `seq 1 $nScalingTestsOpenmp`; do

        if [ $nProcessesCurr -gt ${nProcessesMax} ] ; then
            echo "* WARNING: nProcessesCurr > nProcessesMax. Continue loop."
            continue;
        fi
        if [ $nThreadsCurr -gt ${nThreadsMax} ] ; then
            echo "* WARNING: nThreadsCurr > nThreadsMax. Continue loop."
            continue;
        fi

        export SCAFESRUN_N_THREADS_OPENMP_START=$nThreadsCurr;
        export SCAFESRUN_N_PROCESSES_MPI_START=$nProcessesCurr;
        #----------------------------------------------------------------------#
        if ( expr ${nProcessesCurr} \< ${SCAFESRUN_MACHINE_N_CORES_PER_NODE} > /dev/null ) ;
        then
            JOBSYSTEM_NODES=1;
            JOBSYSTEM_TASKS_PER_NODE=$nProcessesCurr;
        else
            JOBSYSTEM_NODES=$((nProcessesCurr/SCAFESRUN_MACHINE_N_CORES_PER_NODE));
            JOBSYSTEM_TASKS_PER_NODE=$SCAFESRUN_MACHINE_N_CORES_PER_NODE;
        fi
        JOBSYSTEM_CPUS_PER_TASK="$nThreadsCurr"
        JOBSYSTEM_NAME="${SCAFESRUN_NAME_EXECUTABLE}"
        JOBSYSTEM_NAME="${JOBSYSTEM_NODES}_${JOBSYSTEM_NAME}"
        JOBSYSTEM_NAME="${JOBSYSTEM_TASKS_PER_NODE}_${JOBSYSTEM_NAME}"
        JOBSYSTEM_NAME="${JOBSYSTEM_CPUS_PER_TASK}_${JOBSYSTEM_NAME}"
        JOBSYSTEM_JOBFILE="$JOBSYSTEM_NAME.$SCAFESRUN_TYPE_SCALINGTEST.$dateCurr.sh"

        #----------------------------------------------------------------------#
        # Find all software specific environment variables.
        declare -a local allvarsString=$(echo ${!SCAFESRUN*});
        declare -a local allvals=;
        declare local ctr=0;
        for v in ${allvarsString[@]}; do
            varV=$(eval echo "$v")
            valV=$(eval echo "$""$v")
            allvals[$ctr]=$valV;
            # Export all variables, just to be sure.
            export $varV=$valV;
            ctr=$((ctr+1));
        done
        declare -a local allvarsArray=( $(echo ${!SCAFESRUN*}) );
        declare local stringExportVars="";
        # allvarsString, allvarsArray, allvals do all have the same size.
        for ((i=0; i<${#allvarsArray[@]}; ++i)) ;
        do
            stringExportVars="$stringExportVars\nexport ${allvarsArray[i]}=\"${allvals[i]}\";"
        done
        declare local exportedVars=$(echo -e $stringExportVars);

        ########################################################################
        case $SCAFESRUN_RUN_MODE in
            ####################################################################
            LOADLEVELER* )
                cat > $JOBSYSTEM_JOBFILE <<EOF
# @ job_name = $JOBSYSTEM_NAME
# @ error = $JOBSYSTEM_NAME.$dateCurr.out
# @ input = /dev/null
# @ output = $JOBSYSTEM_NAME.$dateCurr.out
# @ shell = /bin/bash
# @ environment = COPY_ALL;
# @ wall_clock_limit = $SCAFESRUN_MACHINE_TIME
# @ notification = never
# @ job_type = bluegene
# @ bg_size = $JOBSYSTEM_NODES
# <Number of nodes>
# @ bg_rotate = TRUE
# <True, False>
# @ bg_connectivity = EITHER
# <TORUS|MESH|EITHER>
# @ queue
$exportedVars

scafesrun.sh
EOF
            set -x
            JOBSYSTEM_JOBID=$(llsubmit $JOBSYSTEM_JOBFILE)
            set +x
            ;;
            ####################################################################
            SLURM* )
                # --dependency=afterany:${JOBSYSTEM_JOBID##* } \
                jobopts="";
                if  [ $nIter > 1 ] ; then
                    jobopts="--dependency=afterany:${JOBSYSTEM_JOBID##* }"
                fi
                set -x
                JOBSYSTEM_JOBID=$(sbatch -A $SCAFESRUN_MACHINE_NAME_PROJECT \
                --job-name=$JOBSYSTEM_NAME \
                $jobopts \
                --mem-per-cpu=$SCAFESRUN_MACHINE_MEMORY_PER_CORE \
                --time=$SCAFESRUN_MACHINE_TIME \
                --nodes=$JOBSYSTEM_NODES \
                --tasks-per-node=$JOBSYSTEM_TASKS_PER_NODE \
                --cpus-per-task=$JOBSYSTEM_CPUS_PER_TASK \
                --partition=$SCAFESRUN_MACHINE_NAME_PARTITION \
                --exclusive \
                scafesrun.sh)
                set +x
                ;;
            ####################################################################
            LSF* )
            # Dirty hack: Change LSF -> NORMAL because job will be already
            # submitted using "bsub" within this section.
                export SCAFESRUN_RUN_MODE="NORMAL";
                jobopts="";
                case $SCAFESRUN_TYPE_SCALINGTEST in
                    *OPENMP ) jobopts="-R span[hosts=1]"; ;;
                esac
                set -x
                JOBSYSTEM_JOBID=$(bsub -P $SCAFESRUN_MACHINE_NAME_PROJECT \
                -J $JOBSYSTEM_NAME \
                -M $SCAFESRUN_MACHINE_MEMORY_PER_CORE \
                -W $SCAFESRUN_MACHINE_TIME \
                -n $JOBSYSTEM_CPUS_PER_TASK \
                $jobopts \
                -o ${JOBSYSTEM_NAME}.lsf.out \
                -e ${JOBSYSTEM_NAME}.lsf.err \
                -x \
                -m $SCAFESRUN_MACHINE_NAME_PARTITION \
                scafesrun.sh)
                set +x
                ;;
            ####################################################################
            * )
                echo " RUN_MODE=$SCAFESRUN_RUN_MODE not known."
                ;;
        esac

        #----------------------------------------------------------------------#
        # Increase number of OpenMP threads or number of MPI processes.
        case $SCAFESRUN_TYPE_SCALINGTEST in
            *OPENMP )
                nThreadsCurr=$((nThreadsCurr*scalingtestFactor)); ;;
            *MPI )
                nProcessesCurr=$((nProcessesCurr*scalingtestFactor)); ;;
        esac

        dir=$((nIter%SCAFESRUN_SPACE_DIM));

        #----------------------------------------------------------------------#
        # Increase number of grid nodes (only for weak scaling tests).
        case $SCAFESRUN_TYPE_SCALINGTEST in
        WEAK* )
            tmpNodes=${nNodesVect[${dir}]}
            nNodesVect[${dir}]=$((scalingtestFactor*tmpNodes));
            export SCAFESRUN_N_NODES="${nNodesVect[0]}";
            for i in `seq 2 $SCAFESRUN_SPACE_DIM`; do
                tmpStringNodes="${SCAFESRUN_N_NODES}x${nNodesVect[$((i-1))]}"
                SCAFESRUN_N_NODES=$tmpStringNodes;
            done
        esac
        # Increase number of MPI partitions (weak and strong scaling tests wrt. MPI)
        case $SCAFESRUN_TYPE_SCALINGTEST in
        *MPI )
            tmpPartitions=${nPartitionsVect[${dir}]}
            nPartitionsVect[${dir}]=$((scalingtestFactor*tmpPartitions));
            export SCAFESRUN_N_PARTITIONS="${nPartitionsVect[0]}";
            for i in `seq 2 $SCAFESRUN_SPACE_DIM`; do
                tmpStringPartitions="${SCAFESRUN_N_PARTITIONS}x${nPartitionsVect[$((i-1))]}"
                SCAFESRUN_N_PARTITIONS=$tmpStringPartitions;
            done
        esac

        nIter=$((nIter+1));
    done
done
)
