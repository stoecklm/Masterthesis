import configparser
import glob
import numpy as np
import os
import subprocess
import sys

from plotSurface import plot_results
from startSimulation import parse_config_file
from startSimulation import check_variables
from startSimulation import calc_variables
from startSimulation import check_stability
from startSimulation import create_init_file
from startSimulation import set_environment_variables
from startSimulation import start_simulation


def parse_config_file_scaling(params):
    print('Parsing {0} for scaling parameters.'.format(params['NAME_CONFIGFILE']))

    # Create configparser and open file.
    config = configparser.ConfigParser()
    config.read(params['NAME_CONFIGFILE'])
    # Get values from section 'Scaling'.
    params['TEST'] = config['Scaling'].get('TEST', fallback='Unknown')
    params['THREADS_OPENMP_START'] = config['Scaling'].getint('THREADS_OPENMP_START', fallback=1)
    params['THREADS_OPENMP'] = params['THREADS_OPENMP_START']
    params['PROCESSES_MPI_START'] = config['Scaling'].getint('PROCESSES_MPI_START', fallback=1)
    params['PROCESSES_MPI'] = params['PROCESSES_MPI_START']
    params['N_TASKS_PER_NODE'] = config['Scaling'].getint('N_TASKS_PER_NODE', fallback=24)
    SCALING = config['Scaling'].get('SCALING', fallback="1")
    params['SCALING'] = list(map(int, SCALING.split(' ')))
    # Get values from section 'Geometry'.
    # Number of nodes.
    N_NODES = config['Geometry'].get('N_NODES')
    params['N_NODES_ENV'] = N_NODES
    params['N_NODES_ENV_START'] = N_NODES
    N_NODES = list(map(int, N_NODES.split('x')))
    params['N_NODES_START'] = N_NODES

    print('Done.')

def check_variables_scaling(params):
    print('Checking scaling variables.')

    if params['TEST'] in ('MPI', 'OpenMP', 'Hybrid'):
        pass
    else:
        print('* ERROR:', params['TEST'], 'is not a valid type of test.')
        print('  Valid types: MPI, OpenMP, or Hybrid.')
        print('Aborting.')
        exit()

    print('Done.')

def calc_variables_scaling(factor, params):
    print('Calculating scaling variables.')

    start_factor = params['SCALING'][0]
    # Calculate new number of nodes.
    N_NODES = params['N_NODES_START']
    params['N_NODES_ENV'] = str(N_NODES[0]*factor) + 'x' + str(N_NODES[1]) + 'x' + str(N_NODES[2])
    params['N_NODES'] = list(map(int, params['N_NODES_ENV'].split('x')))
    # Calculate new number of threads or processes.
    if params['TEST'] == 'MPI' or params['TEST'] == 'Hybrid':
        params['PROCESSES_MPI'] = factor * params['PROCESSES_MPI_START']
    else: # TEST == OpenMP
        params['THREADS_OPENMP'] = factor * params['THREADS_OPENMP_START']
    # Set result dir.
    if params['TEST'] == 'MPI' or params['TEST'] == 'OpenMP':
        RESULT_DIR = './results/weak-scaling/' \
                     + params['TEST'] + '/' \
                     + params['NAME_CONFIGFILE'].split('_')[0] + '/' \
                     + params['N_NODES_ENV_START'] + 'x' \
                     + str(params['N_TIMESTEPS']) + '/' \
                     + str(factor)
    else: # TEST == Hybrid
        RESULT_DIR = './results/weak-scaling/' \
                     + params['TEST'] + '/' \
                     + params['NAME_CONFIGFILE'].split('_')[0] + '/' \
                     + params['N_NODES_ENV_START'] + 'x' \
                     + str(params['N_TIMESTEPS'])  + '/' \
                     + str(params['N_TASKS_PER_NODE']) + 'x' \
                     + str(params['THREADS_OPENMP_START']) + '/' \
                     + str(factor)
    params['RESULT_DIR'] = RESULT_DIR

    print('Done.')

def set_environment_variables_scaling(params):
    print('Setting environment scaling variables.')

    os.putenv('SCAFESRUN_RESULTS_DIR', str(params['RESULT_DIR']))
    print(params['RESULT_DIR'])
    os.putenv('SCAFESRUN_N_THREADS_OPENMP_START', str(params['THREADS_OPENMP']))
    os.putenv('SCAFESRUN_N_PROCESSES_MPI_START', str(params['PROCESSES_MPI']))
    if params['TEST'] == 'MPI' or params['TEST'] == 'Hybrid':
        os.putenv('SCAFESRUN_MACHINE_N_TASKS_PER_NODE', str(params['N_TASKS_PER_NODE']))
    else: # TEST == OpenMP
        os.putenv('SCAFESRUN_MACHINE_N_TASKS_PER_NODE', '1')

    print('Done.')

def main():
    params = {'NAME_CONFIGFILE' : ''}
    params['NAME_RESULTFILE'] = ''
    # Check if path to configfile is provided and if file exists.
    if len(sys.argv) > 1:
        if os.path.isfile(sys.argv[1]) == True:
            params['NAME_CONFIGFILE'] = sys.argv[1]
            run_script = 'RUN_WEAK_SCALINGTEST_HELPER.sh'
        else:
            print('* ERROR:', sys.argv[1], 'does not exist.')
        if len(sys.argv) > 2:
            if os.path.isfile(sys.argv[2]) == True:
                run_script = sys.argv[2]
            else:
                print('* ERROR: Optional run script', sys.argv[2], 'does not exist.')
                print('Aborting.')
                exit()
    else:
        print('* ERROR: No command line argument for configfile provided.')

    if params['NAME_CONFIGFILE'] == '':
        print('Usage: python3', sys.argv[0], '<PATH/TO/CONFIGFILE> [<PATH/TO/RUN/SCRIPT]')
        print('Aborting.')
        exit()

    print('Start weak scaling test.')

    parse_config_file(params)
    parse_config_file_scaling(params)
    check_variables(params)
    check_variables_scaling(params)

    print('Type of test: {0}.'.format(params['TEST']))

    for factor in params['SCALING']:
        calc_variables(params)
        calc_variables_scaling(factor, params)
        check_stability(params)
        if params['CREATE_INITFILE'] == True:
            create_init_file(params)
        set_environment_variables(params)
        set_environment_variables_scaling(params)
        start_simulation(params, run_script)
        if params['NAME_RESULTFILE'] != '' and params['SPACE_DIM'] == 3:
            plot_results(params['NAME_RESULTFILE'])

    print()
    print('Test finished.')

if __name__ == '__main__':
    main()
