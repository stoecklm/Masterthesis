import configparser
import glob
import os
import subprocess
import sys

import netCDF4 as nc
import numpy as np
import matplotlib.path as mpath
from scipy.interpolate import griddata

from plotSurface import plot_surface

from readMRIData import read_intra_op_points
from readMRIData import read_tumor_point
from readMRIData import rotate_points
from readMRIData import move_points
from readMRIData import interpolation
from readMRIData import get_interpolated_path
from readMRIData import get_path

from readMRIVolume import switch_space

from postProcessing import open_surface_temperatures
from postProcessing import tumor_temperatures
from postProcessing import tumor_near_surface_temperatures
from postProcessing import brain_temperatures
from postProcessing import domain_temperatures
from postProcessing import csv_result_temperatures
from postProcessing import vessels_temperatures
from postProcessing import non_vessels_temperatures
from postProcessing import calc_l2_norm

def parse_config_file(params):
    print('Parsing {0}.'.format(params['NAME_CONFIGFILE']))

    # Create configparser and open file.
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(params['NAME_CONFIGFILE'])
    # Get values from section 'Dimension'.
    try:
        params['SPACE_DIM'] = config['Dimension'].getint('SPACE_DIM', fallback=3)
    except KeyError:
        print('* ERROR:', params['NAME_CONFIGFILE'], 'does not contain section \'Dimension\'.')
        print(' ', params['NAME_CONFIGFILE'], 'may not be a config file.')
        print('Aborting.')
        exit()
    # Get values from section 'Geometry'.
    # Coordinates of first node.
    COORD_NODE_FIRST = config['Geometry'].get('COORD_NODE_FIRST')
    params['COORD_NODE_FIRST_ENV'] = COORD_NODE_FIRST
    COORD_NODE_FIRST = list(map(float, COORD_NODE_FIRST.split('x')))
    params['COORD_NODE_FIRST'] = COORD_NODE_FIRST
    # Coordinates of last node.
    COORD_NODE_LAST = config['Geometry'].get('COORD_NODE_LAST')
    params['COORD_NODE_LAST_ENV'] = COORD_NODE_LAST
    COORD_NODE_LAST = list(map(float, COORD_NODE_LAST.split('x')))
    params['COORD_NODE_LAST'] = COORD_NODE_LAST
    # Number of nodes.
    N_NODES = config['Geometry'].get('N_NODES')
    params['N_NODES_ENV'] = N_NODES
    N_NODES = list(map(int, N_NODES.split('x')))
    params['N_NODES'] = N_NODES
    # Get values from section 'Time'.
    params['START_TIME'] = config['Time'].getint('START_TIME', fallback=0)
    params['END_TIME']= config['Time'].getint('END_TIME', fallback=1)
    params['N_TIMESTEPS'] = config['Time'].getint('N_TIMESTEPS', fallback=0)
    # Get values from section 'Output'.
    params['N_SNAPSHOTS'] = config['Output'].getint('N_SNAPSHOTS')
    # Get values from section 'Input'.
    params['USE_MRI_FILE'] = config['Input'].getboolean('USE_MRI_FILE',
                                                        fallback=False)
    params['NAME_REGION_FILE'] = config['Input'].get('NAME_REGION_FILE',
                                                    fallback='region')
    params['NAME_INITFILE'] = config['Input'].get('NAME_INITFILE',
                                                  fallback='init')
    params['USE_INITFILE'] = config['Input'].getboolean('USE_INITFILE',
                                                        fallback=False)
    params['CREATE_INITFILE'] = config['Input'].getboolean('CREATE_INITFILE',
                                                           fallback=False)
    params['NAME_VESSELS_FILE'] = config['Input'].get('NAME_VESSELS_FILE',
                                                      fallback='vessels')
    params['CREATE_VESSELS_FILE'] = config['Input'].getboolean('CREATE_VESSELS_FILE',
                                                               fallback=True)
    params['THRESHOLD'] = config['Input'].getfloat('THRESHOLD',
                                                   fallback=0.00001)
    params['CHECK_CONV_FIRST_AT_ITER'] = config['Input'].getfloat('CHECK_CONV_FIRST_AT_ITER',
                                                                  fallback=1)
    params['CHECK_CONV_AT_EVERY_N_ITER'] = config['Input'].getfloat('CHECK_CONV_AT_EVERY_N_ITER',
                                                                    fallback=1)
    # Get values from section 'MRI'.
    mri_case = config['MRI'].get('CASE', fallback='')
    params['MRI_DATA_CASE'] = mri_case.split('_')[0]
    if params['MRI_DATA_CASE'] != '':
        mri_folder = glob.glob(params['MRI_DATA_CASE'] + '*/')
        if len(mri_folder) == 0:
            print('* ERROR: Folder for case', params['MRI_DATA_CASE'], 'does not exist.')
            print('Aborting.')
            exit()
        params['MRI_DATA_FOLDER'] = mri_folder[0]
    else:
        params['MRI_DATA_FOLDER'] = ''
    params['USE_VESSELS_SEGMENTATION'] = config['MRI'].getboolean('USE_VESSELS_SEGMENTATION',
                                                                  fallback=False)
    VARIABLES_VESSELS = config['MRI'].get('VARIABLES_VESSELS', fallback=list())
    if len(VARIABLES_VESSELS) > 0:
        params['VARIABLES_VESSELS'] = list(VARIABLES_VESSELS.split(' '))
    else:
        params['VARIABLES_VESSELS'] = VARIABLES_VESSELS
    VALUES_VESSELS = config['MRI'].get('VALUES_VESSELS', fallback=list())
    if len(VALUES_VESSELS) > 0:
        params['VALUES_VESSELS'] = list(map(float, VALUES_VESSELS.split(' ')))
    else:
        params['VALUES_VESSELS'] = VALUES_VESSELS
    VALUES_NON_VESSELS = config['MRI'].get('VALUES_NON_VESSELS', fallback=list())
    if len(VALUES_VESSELS) > 0:
        params['VALUES_NON_VESSELS'] = list(map(float, VALUES_NON_VESSELS.split(' ')))
    else:
        params['VALUES_NON_VESSELS'] = VALUES_NON_VESSELS
    params['VESSELS_DEPTH'] = config['MRI'].getint('DEPTH', fallback=1)
    # Get values from section 'Brain'.
    brain = dict(config.items('Brain'))
    for key in brain:
        brain[key] = float(brain[key])
    params['BRAIN'] = brain
    if params['USE_VESSELS_SEGMENTATION'] == True:
        vessels = dict(config.items('Brain'))
        for key in vessels:
            vessels[key] = float(vessels[key])
        params['VESSELS'] = vessels
        non_vessels = dict(config.items('Brain'))
        for key in non_vessels:
            non_vessels[key] = float(non_vessels[key])
        params['NON_VESSELS'] = non_vessels
    # Get values from section 'Tumor'.
    tumor = dict(config.items('Tumor'))
    for key in tumor:
        tumor[key] = float(tumor[key])
    params['TUMOR'] = tumor
    # Get values from section 'Parameters'.
    parameters = dict(config.items('Parameters'))
    for key in parameters:
        parameters[key] = float(parameters[key])
    try:
        parameters['DIAMETER'] = 2.0 * parameters['RADIUS']
    except KeyError:
        pass
    params['PARAMETERS'] = parameters
    # PyMC section.
    try:
        params['ITERATIONS'] = config['PyMC'].getint('ITERATIONS', fallback=5)
        params['BURNS'] = config['PyMC'].getint('BURNS', fallback=1)
        params['T_NORMAL'] = config['PyMC'].getfloat('T_NORMAL', fallback=-1.0)
        params['T_TUMOR'] = config['PyMC'].getfloat('T_TUMOR', fallback=-1.0)
        params['T_VESSEL'] = config['PyMC'].getfloat('T_VESSEL', fallback=-1.0)
    except KeyError:
        params['T_NORMAL'] = -1.0
        params['T_TUMOR'] = -1.0
        params['T_VESSEL'] = -1.0

    print('Done.')

def check_variables(params):
    print('Checking variables.')

    # Check if dimension makes sense and
    # some functions and variables only work for dimension 1, 2 or 3.
    SPACE_DIM = params['SPACE_DIM']
    if SPACE_DIM != 3:
        print('* ERROR: SPACE_DIM is {0}.'.format(SPACE_DIM))
        print('  SPACE_DIM must be 3.')
        print('Aborting.')
        exit()
    # Check if there are enough coordinates for first node.
    DIM_COORD_NODE_FIRST = len(params['COORD_NODE_FIRST'])
    if DIM_COORD_NODE_FIRST != SPACE_DIM:
        print('* ERROR: Dimension of COORD_NODE_FIRST has to be {0}.'.format(SPACE_DIM))
        print('  Dimension of COORD_NODE_FIRST is {0}.'.format(DIM_COORD_NODE_FIRST))
        print('Aborting.')
        exit()
    # Check if there are enough coordinates for last node.
    DIM_COORD_NODE_LAST = len(params['COORD_NODE_LAST'])
    if DIM_COORD_NODE_LAST != SPACE_DIM:
        print('* ERROR: Dimension of COORD_NODE_LAST has to be {0}.'.format(SPACE_DIM))
        print('  Dimension of COORD_NODE_LAST is {0}.'.format(DIM_COORD_NODE_LAST))
        print('Aborting.')
        exit()
    # Check if there are enough number of nodes.
    DIM_N_NODES = len(params['N_NODES'])
    if DIM_N_NODES != SPACE_DIM:
        print('* ERROR: Dimension of N_NODES has to be {0}.'.format(SPACE_DIM))
        print('  Dimension of N_NODES is {0}.'.format(DIM_N_NODES))
        print('Aborting.')
        exit()
    # Check if END_TIME is after START_TIME.
    START_TIME = params['START_TIME']
    END_TIME = params['END_TIME']
    if END_TIME < START_TIME:
        print('* ERROR: END_TIME is smaller than START_TIME.')
        print('  END_TIME must be greater than START_TIME.')
        print('Aborting.')
        exit()
    # Check if threshold is positive.
    if params['THRESHOLD'] < 0.0:
        print('* WARNING: THRESHOLD < 0.0.')
        params['THRESHOLD'] = abs(params['THRESHOLD'])
        print('  THRESHOLD was set to abs(THRESHOLD).')
    # Check if combinations of USE_INITFILE and CREATE_INITFILE makes sense.
    if params['USE_INITFILE'] == True and params['CREATE_INITFILE'] == False:
        if os.path.isfile(params['NAME_INITFILE'] + '.nc') == False:
            print('* ERROR: USE_INITFILE = True and CREATE_INITFILE = False,',
                  'but', params['NAME_INITFILE'] + '.nc', 'does not exist.')
            print('Aborting.')
            exit()
    if params['USE_INITFILE'] == False and params['CREATE_INITFILE'] == True:
        print('* WARNING: CREATE_INITFILE = True, but USE_INITFILE = False.')
    # Check CHECK_CONV parameters.
    if params['CHECK_CONV_FIRST_AT_ITER'] < 0:
        print('* WARNING: CHECK_CONV_FIRST_AT_ITER < 0.')
        params['CHECK_CONV_FIRST_AT_ITER'] = abs(params['CHECK_CONV_FIRST_AT_ITER'])
        print('  CHECK_CONV_FIRST_AT_ITER set to',
              'abs(CHECK_CONV_FIRST_AT_ITER).')
    if params['CHECK_CONV_AT_EVERY_N_ITER'] < 0:
        print('* WARNING: CHECK_CONV_AT_EVERY_N_ITER < 0.')
        params['CHECK_CONV_AT_EVERY_N_ITER'] = abs(params['CHECK_CONV_AT_EVERY_N_ITER'])
        print('  CHECK_CONV_AT_EVERY_N_ITER set to',
              'abs(CHECK_CONV_AT_EVERY_N_ITER).')
    if params['CHECK_CONV_FIRST_AT_ITER'] < 1:
        print('* WARNING: CHECK_CONV_FIRST_AT_ITER < 1.')
        print('  CHECK_CONV_FIRST_AT_ITER is assumend to be a ratio.')
    if params['CHECK_CONV_AT_EVERY_N_ITER'] < 1:
        print('* WARNING: CHECK_CONV_AT_EVERY_N_ITER < 1.')
        print('  CHECK_CONV_AT_EVERY_N_ITER is assumend to be a ratio.')
    # Check if executable exists.
    NAME_EXECUTABLE = os.path.basename(os.getcwd()) \
                      + str(params['SPACE_DIM']) + 'D'
    if os.path.isfile(NAME_EXECUTABLE) == False:
        print(NAME_EXECUTABLE, 'does not exist.')
        print('Aborting.')
        exit()
    params['NAME_EXECUTABLE'] = NAME_EXECUTABLE
    # Check if MRI data exist.
    # Check if path to folder (i.e. results) is provided,
    # and if folder does contain fiducials.csv.
    folder = params['MRI_DATA_FOLDER']
    if folder != '':
        if os.path.isdir(folder) == True:
            tmp1 = os.path.join(folder, 'fiducials.csv')
            tmp2 = os.path.join(folder, 'OpenIGTLink.fcsv')
            if os.path.isfile(tmp1) != True and os.path.isfile(tmp2) != True:
                print('* ERROR:', folder, 'does not contain fiducials.csv',
                      'or OpenIGTLink.fcsv.')
                print('Aborting.')
                exit()
        else:
            print('* ERROR:', folder, 'does not exist.')
            print('Aborting.')
            exit()
    if params['USE_VESSELS_SEGMENTATION'] == True:
        vessels_seg_path = os.path.join(folder, 'vessels_segmentation.csv')
        if os.path.isfile(vessels_seg_path) != True:
            print('* ERROR:', vessels_seg_path, 'does not exist.')
            print('Aborting.')
            exit()
    # Check if file for vessels exist if none shall be created.
    if params['USE_VESSELS_SEGMENTATION'] == True and params['CREATE_VESSELS_FILE'] == False:
        if os.path.isfile(params['NAME_VESSELS_FILE'] + '.nc') == False:
            print('* ERROR: File for vessels does not exist.')
            print('Aborting.')
            exit()
    # Check if names specified in VARIABLES for vessels are
    # variables known in ScaFES.
    names = ['rho', 'c', 'lambda', 'rho_blood', 'c_blood', 'omega', 'T_blood', \
             'q', 'T']
    for var in params['VARIABLES_VESSELS']:
        if var not in names:
            print('* ERROR:', var, 'in VARIABLES_VESSELS not known.')
            print('Aborting.')
            exit()
    if params['VESSELS_DEPTH'] > params['N_NODES'][2]:
        print('* WARNING: Depth for vessel segmentation is bigger than nNodes_2.')
        print('  VESSELS_DEPTH was set to {0}.'.format(params['N_NODES'][2]))
        params['VESSELS_DEPTH'] = params['N_NODES'][2]
    if len(params['VARIABLES_VESSELS']) != len(params['VALUES_VESSELS']):
        print('* ERROR: length of VARIABLES_VESSELS does not match length of',
              'VALUES_VESSELS.')
        print('Aborting.')
        exit()
    if len(params['VARIABLES_VESSELS']) != len(params['VALUES_NON_VESSELS']):
        print('* ERROR: length of VARIABLES_VESSELS does not match length of',
              'VALUES_NON_VESSELS.')
        print('Aborting.')
        exit()

    print('Done.')

def calc_delta_time_helper(params, material, parameters):
    RHO = material['RHO']
    C = material['C']
    LAMBDA = material['LAMBDA']
    RHO_BLOOD = material['RHO_BLOOD']
    C_BLOOD = material['C_BLOOD']
    OMEGA = material['OMEGA']
    T_I = material['T']
    Q = material['Q']

    H = parameters['H']
    EPSILON = parameters['EPSILON']

    GRIDSIZE = params['GRIDSIZE']
    SPACE_DIM = params['SPACE_DIM']

    SIGMA = 5.670367e-8
    T_MAX = T_I + Q/(RHO_BLOOD*C_BLOOD*OMEGA)

    # Pennes Bioheat Equation.
    tmp = 0
    for dim in range(0, SPACE_DIM):
        tmp += (2.0/(GRIDSIZE[dim]*GRIDSIZE[dim])) * (LAMBDA/(RHO*C))
    # Inner nodes.
    tmp += ((RHO_BLOOD*C_BLOOD)/(RHO*C)) * OMEGA
    if tmp != 0:
        DELTA_TIME = 1.0/tmp
    else:
        # If time is infinity,
        # it will later not be considered for min(delta_time).
        DELTA_TIME = float('Inf')

    # Border with convection and thermal radiation:
    # Convection.
    tmp += 2.0*(1.0/GRIDSIZE[SPACE_DIM-1]) * (H/(RHO*C))
    # Thermal radiation.
    tmp += 2.0 * (1.0/GRIDSIZE[SPACE_DIM-1]) \
                     * ((EPSILON*SIGMA)/(RHO*C)) \
                     * ((T_MAX + 273.15)**3)
    if tmp != 0:
        DELTA_TIME_BC = 1.0/tmp
    else:
        # If time is infinity,
        # it will later not be considered for min(delta_time).
        DELTA_TIME_BC = float('Inf')

    return DELTA_TIME, DELTA_TIME_BC

def calc_delta_time_inner_nodes(params, material, parameters):
    tmp,_ = calc_delta_time_helper(params, material, parameters)

    return tmp

def calc_delta_time_boundary_condition(params, material, parameters):
    _,tmp = calc_delta_time_helper(params, material, parameters)

    return tmp

def calc_variables(params):
    print('Calculating variables.')

    # Calculate gridsize in each dimension.
    GRIDSIZE = []
    for dim in range(0, params['SPACE_DIM']):
        GRIDSIZE.append((params['COORD_NODE_LAST'][dim] \
                         - params['COORD_NODE_FIRST'][dim])
                        / (params['N_NODES'][dim]-1))
    params['GRIDSIZE'] = GRIDSIZE
    # Create parameter collection for vessels.
    if params['USE_VESSELS_SEGMENTATION'] == True:
        VARIABLES_VESSELS = params['VARIABLES_VESSELS']
        for NAME_VARIABLE in params['VARIABLES_VESSELS']:
            if NAME_VARIABLE.upper() in params['VESSELS'].keys():
                params['VESSELS'][NAME_VARIABLE.upper()] = params['VALUES_VESSELS'][VARIABLES_VESSELS.index(NAME_VARIABLE)]
            if NAME_VARIABLE.upper() in params['NON_VESSELS'].keys():
                params['NON_VESSELS'][NAME_VARIABLE.upper()] = params['VALUES_NON_VESSELS'][VARIABLES_VESSELS.index(NAME_VARIABLE)]
    # Calculate delta time.
    if params['N_TIMESTEPS'] < 1:
        print('* WARNING: N_TIMESTEPS not specified.')
        print('  Calculate N_TIMESTEPS from stability criterion.')
        BRAIN = calc_delta_time_inner_nodes(params, params['BRAIN'],
                                            params['PARAMETERS'])
        BRAIN_BC = calc_delta_time_boundary_condition(params, params['BRAIN'],
                                                      params['PARAMETERS'])
        TUMOR = calc_delta_time_inner_nodes(params, params['TUMOR'],
                                            params['PARAMETERS'])
        TUMOR_BC = calc_delta_time_boundary_condition(params, params['TUMOR'],
                                                      params['PARAMETERS'])
        if params['USE_VESSELS_SEGMENTATION'] == True:
            VESSELS = calc_delta_time_inner_nodes(params, params['VESSELS'],
                                                  params['PARAMETERS'])
            VESSELS_BC = calc_delta_time_boundary_condition(params,
                                                            params['VESSELS'],
                                                            params['PARAMETERS'])
            NON_VESSELS = calc_delta_time_inner_nodes(params, params['NON_VESSELS'],
                                                      params['PARAMETERS'])
            NON_VESSELS_BC = calc_delta_time_boundary_condition(params,
                                                                params['NON_VESSELS'],
                                                                params['PARAMETERS'])
        else:
            VESSELS = float('Inf')
            VESSELS_BC = float('Inf')
            NON_VESSELS = float('Inf')
            NON_VESSELS_BC = float('Inf')

        # Get minimum for calculation of timesteps.
        DELTA_TIME_MIN = min((BRAIN, BRAIN_BC, TUMOR, TUMOR_BC,
                              VESSELS, VESSELS_BC, NON_VESSELS, NON_VESSELS_BC))
        # Add five percent for safety reasons.
        params['N_TIMESTEPS'] = int(((params['END_TIME'] \
                                      - params['START_TIME']) \
                                     / DELTA_TIME_MIN) * 1.05) + 1

    # Final calculation for delta time.
    params['DELTA_TIME'] = (params['END_TIME'] - params['START_TIME']) \
                           / params['N_TIMESTEPS']

    # Calculate location of tumor center.
    TUMOR_CENTER = []
    TUMOR_CENTER.append((params['COORD_NODE_LAST'][0] \
                         + params['COORD_NODE_FIRST'][0]) / 2.0)
    TUMOR_CENTER.append((params['COORD_NODE_LAST'][1] \
                         + params['COORD_NODE_FIRST'][1]) / 2.0)
    TUMOR_CENTER.append(params['COORD_NODE_LAST'][2]
                        - params['PARAMETERS']['DEPTH'])
    params['TUMOR_CENTER'] = TUMOR_CENTER

    # Calc CHECK_CONV parameters if they are a ratio.
    if params['CHECK_CONV_FIRST_AT_ITER'] < 1:
        params['CHECK_CONV_FIRST_AT_ITER'] = params['CHECK_CONV_FIRST_AT_ITER'] \
                                             * params['N_TIMESTEPS']
    params['CHECK_CONV_FIRST_AT_ITER'] = int(params['CHECK_CONV_FIRST_AT_ITER'])
    if params['CHECK_CONV_AT_EVERY_N_ITER'] < 1:
        params['CHECK_CONV_AT_EVERY_N_ITER'] = params['CHECK_CONV_AT_EVERY_N_ITER'] \
                                               * params['N_TIMESTEPS']
    params['CHECK_CONV_AT_EVERY_N_ITER'] = int(params['CHECK_CONV_AT_EVERY_N_ITER'])

    # Check if number of snapshots is possible.
    if params['N_SNAPSHOTS'] > params['N_TIMESTEPS']:
        print('* WARNING: N_SNAPSHOTS was bigger than N_TIMESTEPS.')
        params['N_SNAPSHOTS'] = params['N_TIMESTEPS']
        print('  N_SNAPSHOTS was set to N_TIMESTEPS.')

    print('Done.')

def check_stability(params):
    print('Checking stability.')

    BRAIN = calc_delta_time_inner_nodes(params, params['BRAIN'],
                                        params['PARAMETERS'])
    BRAIN_BC = calc_delta_time_boundary_condition(params, params['BRAIN'],
                                                  params['PARAMETERS'])
    TUMOR = calc_delta_time_inner_nodes(params, params['TUMOR'],
                                                  params['PARAMETERS'])
    TUMOR_BC = calc_delta_time_boundary_condition(params, params['TUMOR'],
                                                  params['PARAMETERS'])
    if params['USE_VESSELS_SEGMENTATION'] == True:
        VESSELS = calc_delta_time_inner_nodes(params, params['VESSELS'],
                                              params['PARAMETERS'])
        VESSELS_BC = calc_delta_time_boundary_condition(params,
                                                        params['VESSELS'],
                                                        params['PARAMETERS'])
        NON_VESSELS = calc_delta_time_inner_nodes(params, params['NON_VESSELS'],
                                                  params['PARAMETERS'])
        NON_VESSELS_BC = calc_delta_time_boundary_condition(params,
                                                            params['NON_VESSELS'],
                                                            params['PARAMETERS'])
    else:
        VESSELS = float('Inf')
        VESSELS_BC = float('Inf')
        NON_VESSELS = float('Inf')
        NON_VESSELS_BC = float('Inf')

    # Get minimum for calculation of timesteps.
    DELTA_TIME_MIN = min((BRAIN, BRAIN_BC, TUMOR, TUMOR_BC,
                          VESSELS, VESSELS_BC, NON_VESSELS, NON_VESSELS_BC))

    DELTA_TIME = params['DELTA_TIME']
    # Abort simulation if stability is not fulfilled.
    if DELTA_TIME > BRAIN:
        print('* ERROR: Stability not fulfilled in healthy brain region.')
        print('  DELTA_TIME = {0}, but has to be DELTA_TIME < {1}.'.format(DELTA_TIME,
                                                                           BRAIN))
        print('Aborting.')
        exit()
    if DELTA_TIME > TUMOR:
        print('* ERROR: Stability not fulfilled in tumor region.')
        print('  DELTA_TIME = {0}, but has to be DELTA_TIME < {1}.'.format(DELTA_TIME,
                                                                           TUMOR))
        print('Aborting.')
        exit()
    if DELTA_TIME > BRAIN_BC:
        print('* ERROR: Stability not fulfilled in healty brain region at \
              border with convection and thermal radiation.')
        print('  DELTA_TIME = {0}, but has to be DELTA_TIME < {1}.'.format(DELTA_TIME,
                                                                           BRAIN_BC))
        print('Aborting.')
        exit()
    if DELTA_TIME > TUMOR_BC:
        print('* ERROR: Stability not fulfilled in tumor region at border \
              with convection and thermal radiation.')
        print('  DELTA_TIME = {0}, but has to be DELTA_TIME < {1}.'.format(DELTA_TIME,
                                                                           TUMOR_BC))
        print('Aborting.')
        exit()
    if params['USE_VESSELS_SEGMENTATION'] == True:
        if DELTA_TIME > VESSELS:
            print('* ERROR: Stability not fulfilled in vessels region.')
            print('  DELTA_TIME = {0}, but has to be DELTA_TIME < {1}.'.format(DELTA_TIME,
                                                                               VESSELS))
            print('Aborting.')
            exit()
        if DELTA_TIME > NON_VESSELS:
            print('* ERROR: Stability not fulfilled in non-vessels region.')
            print('  DELTA_TIME = {0}, but has to be DELTA_TIME < {1}.'.format(DELTA_TIME,
                                                                               NON_VESSELS))
            print('Aborting.')
            exit()
        if DELTA_TIME > VESSELS_BC:
            print('* ERROR: Stability not fulfilled in vessels region at \
                  border with convection and thermal radiation.')
            print('  DELTA_TIME = {0}, but has to be DELTA_TIME < {1}.'.format(DELTA_TIME,
                                                                               VESSELS_BC))
            print('Aborting.')
            exit()
        if DELTA_TIME > NON_VESSELS_BC:
            print('* ERROR: Stability not fulfilled in non-vessels region at \
                  border with convection and thermal radiation.')
            print('  DELTA_TIME = {0}, but has to be DELTA_TIME < {1}.'.format(DELTA_TIME,
                                                                               NON_VESSELS_BC))
            print('Aborting.')
            exit()

    print('Done.')

def create_region_array(params, nc_file, BRAIN_VALUE, TUMOR_VALUE,
                        NAME_VARIABLE):
    RADIUS = params['PARAMETERS']['DIAMETER']/2
    # Get file/grid dimensions.
    dim0, dim1, dim2 = params['N_NODES']
    COORD_NODE_FIRST = params['COORD_NODE_FIRST']
    # Get tumor center location.
    TUMOR_CENTER = params['TUMOR_CENTER']
    num_elem = dim0 * dim1 * dim2
    values_array = BRAIN_VALUE \
                   * np.ones(num_elem, dtype=int).reshape(dim2, dim1, dim0)
    # Iterate through array.
    for elem_z in range(0, values_array.shape[0]):
        for elem_y in range(0, values_array.shape[1]):
            for elem_x in range(0, values_array.shape[2]):
                # Calculate location of current node.
                x = (elem_x * params['GRIDSIZE'][0]) + COORD_NODE_FIRST[0]
                y = (elem_y * params['GRIDSIZE'][1]) + COORD_NODE_FIRST[1]
                z = (elem_z * params['GRIDSIZE'][2]) + COORD_NODE_FIRST[2]
                # Calculate distance (squared) to tumor center.
                distance = (x - TUMOR_CENTER[0]) * (x - TUMOR_CENTER[0])
                distance += (y - TUMOR_CENTER[1]) * (y - TUMOR_CENTER[1])
                distance += (z - TUMOR_CENTER[2]) * (z - TUMOR_CENTER[2])
                # Check if current point is inside tumor.
                # If yes, set value to tumor specific value
                if distance <= RADIUS*RADIUS:
                    values_array[elem_z, elem_y, elem_x] = TUMOR_VALUE
    # Create netCDF variable.
    nNodes = []
    nNodes.append('time')
    for dim in range(len(values_array.shape), 0, -1):
        nNodes.append('nNodes_' + str(dim-1))
    init_values = nc_file.createVariable(NAME_VARIABLE, 'i', nNodes)
    # Write NumPy Array to file.
    init_values[0,] = values_array

def create_region_file(params):
    filepath = params['NAME_REGION_FILE'] + '.nc'
    SPACE_DIM = params['SPACE_DIM']
    print('Creating {0}.'.format(filepath))

    # Delete old region file.
    if os.path.isfile(filepath) == True:
        os.remove(filepath)

    nc_file = nc.Dataset(filepath, 'w', format='NETCDF3_CLASSIC')
    time = nc_file.createDimension('time')
    for dim in range(0, SPACE_DIM):
        nNodes = nc_file.createDimension('nNodes_' + str(dim),
                                         params['N_NODES'][dim])

    # 0 means brain, 1 means tumor.
    create_region_array(params, nc_file, 0, 1, 'region')

    nc_file.close()

    print('Done.')

def write_values_to_file(nc_file, values_array, NAME_VARIABLE):
    # Create netCDF variable.
    nNodes = []
    nNodes.append('time')
    for dim in range(len(values_array.shape), 0, -1):
        nNodes.append('nNodes_' + str(dim-1))
    init_values = nc_file.createVariable(NAME_VARIABLE, 'f8', nNodes)
    # Write NumPy Array to file.
    init_values[0,] = values_array

def create_vessels_array(params, surface):
    print('Creating {0}.nc.'.format(params['NAME_VESSELS_FILE']))
    vessels_small = read_vessels_segmentation(params)

    dim0, dim1, dim2 = params['N_NODES']
    num_elem = dim0 * dim1 * dim2
    # Special Case: No trepanation domain is set,
    # but vessel segmentation is read.
    if np.count_nonzero(surface) == 0:
        print('* WARNING: No trepanation area is set, but vessel segmentation is read.')
        print('  Vessels can only be created in trepanation area.')
        print('  File will contain no vessels.')
        surface[-1,:,:] = 0
    # Normal case: trepanation domain is set.
    # - 1 = grid node outside of trepanation domain
    # 0 = grid node inside trepanation domain, no vessel
    # 1 = grid node is vessel inside trepanation domain
    vessels_big = np.ones(dim1*dim0).reshape(dim1, dim0)
    vessels_big *= -1.0
    x_min = params['surface_cmin']
    x_max = params['surface_cmax']
    y_min = params['surface_rmin']
    y_max = params['surface_rmax']
    depth = params['VESSELS_DEPTH']
    surface = surface[-1,:,:]
    vessels_tmp = np.zeros(dim1*dim0).reshape(dim1, dim0)
    vessels_tmp[y_min:y_max+1,x_min:x_max+1] = vessels_small[:,:]
    vessels_big = np.where(surface == 1, vessels_tmp, vessels_big)
    vessels_big = np.repeat(vessels_big[np.newaxis,:,:], depth, axis=0)
    vessels = np.ones(dim2*dim1*dim0).reshape(dim2, dim1, dim0)
    vessels *= -1.0
    vessels[-depth:,:,:] = vessels_big

    # Create vessels file.
    filepath = params['NAME_VESSELS_FILE'] + '.nc'
    nc_file = nc.Dataset(filepath, 'w', format='NETCDF3_CLASSIC')
    time = nc_file.createDimension('time')
    for dim in range(0, params['SPACE_DIM']):
        nNodes = nc_file.createDimension('nNodes_' + str(dim),
                                         params['N_NODES'][dim])
    write_values_to_file(nc_file, vessels, 'vessels')

    nc_file.close()

    print('Done')

    return vessels

def create_init_array(params, nc_file, region, BRAIN_VALUE, TUMOR_VALUE,
                      NAME_VARIABLE, vessels, surface):
    dim0, dim1, dim2 = params['N_NODES']
    num_elem = dim0 * dim1 * dim2
    values_array = BRAIN_VALUE * np.ones(num_elem).reshape(dim2, dim1, dim0)

    if params['USE_VESSELS_SEGMENTATION'] == True:
        VARIABLES_VESSELS = params['VARIABLES_VESSELS']
        if NAME_VARIABLE in VARIABLES_VESSELS:
            VALUE_VESSEL = params['VALUES_VESSELS'][VARIABLES_VESSELS.index(NAME_VARIABLE)]
            VALUE_NON_VESSEL = params['VALUES_NON_VESSELS'][VARIABLES_VESSELS.index(NAME_VARIABLE)]
            values_array = np.where(vessels == 1, VALUE_VESSEL, values_array)
            values_array = np.where(vessels == 0, VALUE_NON_VESSEL, values_array)

    values_array = np.where(region == 1, TUMOR_VALUE, values_array)
    write_values_to_file(nc_file, values_array, NAME_VARIABLE)

def create_surface_array(params, nc_file, BRAIN_VALUE, TUMOR_VALUE,
                         NAME_VARIABLE):
    RADIUS = (params['PARAMETERS']['DIAMETER'] \
              * params['PARAMETERS']['HOLE_FACTOR'])/2
    # Get file/grid dimensions.
    dim0, dim1, dim2 = params['N_NODES']
    COORD_NODE_FIRST = params['COORD_NODE_FIRST']
    # Get tumor center location.
    TUMOR_CENTER = params['TUMOR_CENTER']
    # Resize array.
    num_elem = dim0 * dim1 * dim2
    values_array = BRAIN_VALUE \
                   * np.ones(num_elem, dtype=int).reshape(dim2, dim1, dim0)
    # Iterate through array.
    for elem_y in range(0, values_array.shape[1]):
        for elem_x in range(0, values_array.shape[2]):
            # Calculate location of current node.
            x = (elem_x * params['GRIDSIZE'][0]) + COORD_NODE_FIRST[0]
            y = (elem_y * params['GRIDSIZE'][1]) + COORD_NODE_FIRST[1]
            # Calculate distance (squared) to tumor center.
            distance = (x - TUMOR_CENTER[0]) * (x - TUMOR_CENTER[0])
            distance += (y - TUMOR_CENTER[1]) * (y - TUMOR_CENTER[1])
            # Check if current point is inside tumor.
            # If yes, set value to tumor specific value
            if distance <= RADIUS*RADIUS:
                values_array[-1, elem_y, elem_x] = TUMOR_VALUE
    # Create netCDF variable.
    nNodes = []
    nNodes.append('time')
    for dim in range(len(values_array.shape), 0, -1):
        nNodes.append('nNodes_' + str(dim-1))
    init_values = nc_file.createVariable(NAME_VARIABLE, 'i', nNodes)
    # Write NumPy Array to file.
    init_values[0,] = values_array

    # Bounding box for trepanation domain.
    rows = np.any(values_array[-1,:,:], axis=1)
    cols = np.any(values_array[-1,:,:], axis=0)
    try:
        params['surface_rmin'], params['surface_rmax'] = np.where(rows)[0][[0, -1]]
    except IndexError:
        params['surface_rmin'], params['surface_rmax'] = 0, dim1-1
    try:
        params['surface_cmin'], params['surface_cmax'] = np.where(cols)[0][[0, -1]]
    except IndexError:
        params['surface_cmin'], params['surface_cmax'] = 0, dim0-1

    return values_array

def create_surface_from_mri(params, nc_file, BRAIN_VALUE, TUMOR_VALUE,
                            NAME_VARIABLE):
    # Get file/grid dimensions.
    dim0, dim1, dim2 = params['N_NODES']
    COORD_NODE_FIRST = params['COORD_NODE_FIRST']
    # Resize array.
    num_elem = dim0 * dim1 * dim2
    values_array = BRAIN_VALUE \
                   * np.ones(num_elem, dtype=int).reshape(dim2, dim1, dim0)

    filepath = params['MRI_DATA_FOLDER']
    iop = read_intra_op_points(filepath)
    t, contains_tumor = read_tumor_point(filepath)
    if contains_tumor == False:
        print('* ERROR: No tumor coordinates in intraoperative data of case', params['MRI_DATA_FOLDER'], 'found.')
        print('  Tumor geometry is not build from MRI data.')
        print('  Tumor center is necessary for building trepanation area correctly.')
        print('  Try other cases or use synthetic data.')
        print('Aborting.')
        exit()
    for point in iop:
        point[...] = switch_space(point)
    t = switch_space(t)
    iop, t = rotate_points(iop, t)
    if params['USE_MRI_FILE'] == False:
        iop, t = move_points(iop, t, t)
        path = get_interpolated_path(iop)
        # Get tumor center location.
        TUMOR_CENTER = params['TUMOR_CENTER']
        # Get points to define open skull.
        pts = interpolation(iop)
        # mm to m.
        pts[:,0] /= 1000
        pts[:,1] /= 1000
        # Transform/move points according to new coordinate system.
        pts[:,0] += TUMOR_CENTER[0]
        pts[:,1] += TUMOR_CENTER[1]
        params['HOLE'] = pts
        # Iterate through array.
        for elem_y in range(0, dim1):
            for elem_x in range(0, dim0):
                # Calculate location of current node.
                x = (elem_x * params['GRIDSIZE'][0]) + COORD_NODE_FIRST[0]
                y = (elem_y * params['GRIDSIZE'][1]) + COORD_NODE_FIRST[1]
                # Transform current node to tumor center system.
                x_trans = (x - TUMOR_CENTER[0])*1000
                y_trans = (y - TUMOR_CENTER[1])*1000
                # Check if current point is inside open skill
                # If yes, set value to tumor specific value.
                if path.contains_point((x_trans,y_trans)) == True:
                    values_array[-1, elem_y, elem_x] = TUMOR_VALUE

    if params['USE_MRI_FILE'] == True:
        path = get_interpolated_path(iop)
        # Get points to define open skull.
        pts = interpolation(iop)
        # mm to m.
        pts[:,0] /= 1000
        pts[:,1] /= 1000
        params['HOLE'] = pts
        # Iterate through array.
        for elem_y in range(0, dim1):
            for elem_x in range(0, dim0):
                # Calculate location of current node.
                x = (elem_x * params['GRIDSIZE'][0]) + COORD_NODE_FIRST[0]
                y = (elem_y * params['GRIDSIZE'][1]) + COORD_NODE_FIRST[1]
                # Transform current node to tumor center system.
                x_trans = x * 1000
                y_trans = y * 1000
                # Check if current point is inside open skill
                # If yes, set value to tumor specific value.
                if path.contains_point((x_trans,y_trans)) == True:
                    values_array[-1, elem_y, elem_x] = TUMOR_VALUE

    # Create netCDF variable.
    nNodes = []
    nNodes.append('time')
    for dim in range(len(values_array.shape), 0, -1):
        nNodes.append('nNodes_' + str(dim-1))
    init_values = nc_file.createVariable(NAME_VARIABLE, 'i', nNodes)
    # Write NumPy Array to file.
    init_values[0,] = values_array

    # Bounding box for trepanation domain.
    rows = np.any(values_array[-1,:,:], axis=1)
    cols = np.any(values_array[-1,:,:], axis=0)
    params['surface_rmin'], params['surface_rmax'] = np.where(rows)[0][[0, -1]]
    params['surface_cmin'], params['surface_cmax'] = np.where(cols)[0][[0, -1]]

    return values_array

def read_vessels_segmentation(params):
    print('Read vessels segmentation.')
    # Load vessels segmentation and save it with the smallest bounding box.
    vessels_seg_path = os.path.join(params['MRI_DATA_FOLDER'],
                                    'vessels_segmentation.csv')
    a = np.genfromtxt(vessels_seg_path, delimiter=',')
    rows = np.any(a, axis=1)
    cols = np.any(a, axis=0)
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]
    b = np.zeros((rmax-rmin+1, cmax-cmin+1))
    b[:,:] = a[rmin:rmax+1,cmin:cmax+1]

    dim0_sparse = params['surface_cmax'] - params['surface_cmin'] + 1
    dim1_sparse = params['surface_rmax'] - params['surface_rmin'] + 1

    dim0_dense = b.shape[1]
    dim1_dense = b.shape[0]

    x_sparse, y_sparse = np.meshgrid(np.linspace(0, dim0_sparse-1, dim0_sparse),
                                     np.linspace(0, dim1_sparse-1, dim1_sparse))
    x_dense, y_dense = np.meshgrid(np.linspace(0, dim0_sparse-1, dim0_dense),
                                   np.linspace(0, dim1_sparse-1, dim1_dense))

    vessels_sparse = griddata(np.array([x_dense.ravel(), y_dense.ravel()]).T,
                              b.ravel(), (x_sparse, y_sparse), method='nearest')

    print('Done.')

    return vessels_sparse

def read_vessels_array(params):
    print('Reading {0}.nc.'.format(params['NAME_VESSELS_FILE']))
    nc_file = nc.Dataset(params['NAME_VESSELS_FILE'] + '.nc')
    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size
    nc_var = nc_file.variables['vessels']
    vessels = np.zeros((dim2, dim1, dim0), dtype=int)
    vessels[:,:,:] = nc_var[-1,:,:,:]

    nc_file.close()

    print('Done.')

    return vessels

def create_init_file(params):
    filepath = params['NAME_INITFILE'] + '.nc'
    SPACE_DIM = params['SPACE_DIM']
    print('Creating {0}.'.format(filepath))

    # Delete old init file.
    if os.path.isfile(filepath) == True:
        os.remove(filepath)

    # Check if region file exists.
    if os.path.isfile(params['NAME_REGION_FILE'] + '.nc') == False:
        print('* ERROR: File for region does not exist.')
        print('Aborting.')
        exit()

    # Open region file.
    nc_file = nc.Dataset(params['NAME_REGION_FILE'] + '.nc')
    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size
    nc_var = nc_file.variables['region']
    region = np.zeros((dim2, dim1, dim0), dtype=int)
    region[:,:,:] = nc_var[-1,:,:,:]

    nc_file.close()

    # Create init file.
    nc_file = nc.Dataset(filepath, 'w', format='NETCDF3_CLASSIC')
    time = nc_file.createDimension('time')
    for dim in range(0, SPACE_DIM):
        nNodes = nc_file.createDimension('nNodes_' + str(dim),
                                         params['N_NODES'][dim])

    if params['MRI_DATA_CASE'] != '':
        surface = create_surface_from_mri(params, nc_file, 0, 1, 'surface')
    else:
        surface = create_surface_array(params, nc_file, 0, 1, 'surface')

    if params['USE_VESSELS_SEGMENTATION'] == True and params['CREATE_VESSELS_FILE'] == True:
        vessels = create_vessels_array(params, surface)
    elif params['USE_VESSELS_SEGMENTATION'] == True and params['CREATE_VESSELS_FILE'] == False:
        vessels = read_vessels_array(params)
    else:
        vessels = 0

    brain = params['BRAIN']
    tumor = params['TUMOR']
    names = {'RHO': 'rho', 'C': 'c', 'LAMBDA': 'lambda',
             'RHO_BLOOD': 'rho_blood', 'C_BLOOD': 'c_blood', 'OMEGA': 'omega',
             'T_BLOOD': 'T_blood', 'Q': 'q', 'T': 'T'}
    for key, value in brain.items():
        create_init_array(params, nc_file, region, brain[key], tumor[key],
                          names[key], vessels, surface)

    nc_file.close()

    print('Done.')

def set_environment_variables(params):
    print('Setting environment variables.')

    # Set all environment from dict.
    os.putenv('SCAFESRUN_N_TIMESTEPS', str(params['N_TIMESTEPS']))
    os.putenv('SCAFESRUN_N_SNAPSHOTS', str(params['N_SNAPSHOTS']))
    os.putenv('SCAFESRUN_START_TIME', str(params['START_TIME']))
    os.putenv('SCAFESRUN_END_TIME', str(params['END_TIME']))
    os.putenv('SCAFESRUN_NAME_CONFIGFILE', params['NAME_CONFIGFILE'])
    os.putenv('SCAFESRUN_SPACE_DIM', str(params['SPACE_DIM']))
    os.putenv('SCAFESRUN_COORD_NODE_FIRST', str(params['COORD_NODE_FIRST_ENV']))
    os.putenv('SCAFESRUN_COORD_NODE_LAST', str(params['COORD_NODE_LAST_ENV']))
    os.putenv('SCAFESRUN_N_NODES', str(params['N_NODES_ENV']))
    os.putenv('SCAFESRUN_NAME_EXECUTABLE', str(params['NAME_EXECUTABLE']))
    os.putenv('SCAFESRUN_THRESHOLD', str(params['THRESHOLD']))
    os.putenv('SCAFESRUN_CHECK_CONV_FIRST_AT_ITER',
              str(params['CHECK_CONV_FIRST_AT_ITER']))
    os.putenv('SCAFESRUN_CHECK_CONV_AT_EVERY_N_ITER',
              str(params['CHECK_CONV_AT_EVERY_N_ITER']))
    # Check if init file should be used and if it exists.
    if params['USE_INITFILE'] == True:
        if os.path.isfile(params['NAME_INITFILE'] + '.nc') == True:
            os.putenv('SCAFESRUN_NAME_INITFILE', params['NAME_INITFILE'])
        else:
            print('* ERROR: USE_INITFILE = True, but',
                  params['NAME_INITFILE'] + '.nc', 'does not exist.')
            print('Aborting.')
            exit()

    print('Done.')

def call_simulation(params, run_script):
    # Check if run sript exists.
    if os.path.isfile(run_script) == False:
        print('* ERROR:', run_script, 'does not exist.')
        print('Aborting.')
        exit()

    # Get time of newest netCDF file BEFORE simulation.
    inifile = params['NAME_CONFIGFILE'].split('.')[0]
    inifile = inifile.split('_')[0]
    searchpath = './' + params['NAME_EXECUTABLE'] + '_' + inifile + '_' + '*.nc'
    files_nc_before = glob.glob(searchpath)
    if len(files_nc_before) > 0:
        latest_file_nc_before = max(files_nc_before, key=os.path.getctime)
        latest_file_nc_before_time = os.path.getctime(latest_file_nc_before)
    else:
        latest_file_nc_before_time = 0

    # Set name of folder as executable with dimension as suffix
    print('Starting {0}.'.format(params['NAME_EXECUTABLE']))
    print()
    # Call bash script to set more environment variables and
    # to start simulation.
    returncode = subprocess.call('./' + run_script)
    # Check if simulation/bash script ran successfully.
    print()
    if returncode == 0:
        print('Done.')
    else:
        print('* ERROR: Simulation returned error code {0}.'.format(returncode))
        print('Aborting.')
        exit()

    # Get time of newest netCDF file AFTER simulation.
    files_nc_after = glob.glob(searchpath)
    if len(files_nc_after) > 0:
        latest_file_nc_after = max(files_nc_after, key=os.path.getctime)
        latest_file_nc_after_time = os.path.getctime(latest_file_nc_after)
    else:
        latest_file_nc_after_time = 0

    # If time of newest file after simulation is newer than time of newest file
    # before simulation, then it is assumend there is a new file written
    # by the simulation.
    if latest_file_nc_after_time > latest_file_nc_before_time:
        params['NAME_RESULTFILE'] = latest_file_nc_after
    else:
        params['NAME_RESULTFILE'] = ''


def main():
    params = {'NAME_CONFIGFILE' : ''}
    params['NAME_RESULTFILE'] = ''
    # Check if path to configfile is provided and if file exists.
    if len(sys.argv) > 1:
        if os.path.isfile(sys.argv[1]) == True:
            if sys.argv[1].endswith('.ini'):
                params['NAME_CONFIGFILE'] = sys.argv[1]
                run_script = 'RUN_HELPER.sh'
            else:
                print('* ERROR:', sys.argv[1], 'may not be a config file.')
                print('Aborting.')
                exit()
        elif os.path.isdir(sys.argv[1]) == True:
            print('* ERROR:', sys.argv[1], 'is not a file.')
        else:
            print('* ERROR:', sys.argv[1], 'does not exist.')
        if len(sys.argv) > 2:
            if os.path.isfile(sys.argv[2]) == True:
                run_script = sys.argv[2]
            else:
                print('* ERROR: Optional run script', sys.argv[2],
                      'does not exist.')
                print('Aborting.')
                exit()
    else:
        print('* ERROR: No command line argument for configfile provided.')

    if params['NAME_CONFIGFILE'] == '':
        print('Usage: python3', sys.argv[0],
              '<PATH/TO/CONFIGFILE> [<PATH/TO/RUN/SCRIPT]')
        print('Aborting.')
        exit()

    parse_config_file(params)
    check_variables(params)
    calc_variables(params)
    check_stability(params)
    if params['USE_MRI_FILE'] == False:
        create_region_file(params)
    create_init_file(params)
    set_environment_variables(params)
    call_simulation(params, run_script)
    if params['NAME_RESULTFILE'] != '' and params['SPACE_DIM'] == 3:
        plot_surface(params['NAME_RESULTFILE'], params)
        open_surface_temperatures(params['NAME_RESULTFILE'],
                                  params['NAME_INITFILE'])
        tumor_temperatures(params['NAME_RESULTFILE'],
                           params['NAME_REGION_FILE'])
        T_tumor = tumor_near_surface_temperatures(params['NAME_RESULTFILE'],
                                                  params['NAME_REGION_FILE'])
        brain_temperatures(params['NAME_RESULTFILE'],
                           params['NAME_REGION_FILE'])
        domain_temperatures(params['NAME_RESULTFILE'])
        if params['USE_VESSELS_SEGMENTATION'] == True:
            T_vessel = vessels_temperatures(params['NAME_RESULTFILE'],
                                            params['NAME_VESSELS_FILE'])
            T_normal = non_vessels_temperatures(params['NAME_RESULTFILE'],
                                               params['NAME_VESSELS_FILE'])
        else:
            T_vessel = -1.0
            T_normal = -1.0
        if params['MRI_DATA_CASE'] != '':
            csv_result_temperatures(params['NAME_RESULTFILE'],
                                    params['MRI_DATA_FOLDER'])
        calc_l2_norm(params['NAME_RESULTFILE'], T_normal, T_tumor, T_vessel,
                     params['T_NORMAL'], params['T_TUMOR'], params['T_VESSEL'])


if __name__ == '__main__':
    main()
