import configparser
import glob
import netCDF4 as nc
import numpy as np
import os
import subprocess
import sys
import matplotlib.path as mpath

from plotSurface import plot_surface
from readMRIData import read_intra_op_points
from readMRIData import read_tumor_point
from readMRIData import rotate_points
from readMRIData import move_points
from readMRIData import interpolation
from readMRIData import get_interpolated_path
from readMRIData import get_path

from postProcessing import mean_surface_temperature
from postProcessing import mean_tumor_temperature

def parse_config_file(params):
    print('Parsing {0}.'.format(params['NAME_CONFIGFILE']))

    # Create configparser and open file.
    config = configparser.ConfigParser()
    config.read(params['NAME_CONFIGFILE'])
    # Get values from section 'Dimension'.
    params['SPACE_DIM'] = config['Dimension'].getint('SPACE_DIM', fallback=3)
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
    params['NAME_INITFILE'] = config['Input'].get('NAME_INITFILE',
                                                  fallback='init')
    params['USE_INITFILE'] = config['Input'].getboolean('USE_INITFILE',
                                                        fallback=False)
    params['CREATE_INITFILE'] = config['Input'].getboolean('CREATE_INITFILE',
                                                           fallback=False)
    params['THRESHOLD'] = config['Input'].getfloat('THRESHOLD',
                                                   fallback=0.00001)
    params['CHECK_CONV_FIRST_AT_ITER'] = config['Input'].getfloat('CHECK_CONV_FIRST_AT_ITER',
                                                                  fallback=1)
    params['CHECK_CONV_AT_EVERY_N_ITER'] = config['Input'].getfloat('CHECK_CONV_AT_EVERY_N_ITER',
                                                                    fallback=1)
    # Get values from section 'MRI'.
    params['MRI_DATA_CASE'] = config['MRI'].get('CASE', fallback='')
    # Get values from section 'Brain'.
    brain = dict(config.items('Brain'))
    for key in brain:
        brain[key] = float(brain[key])
    params['BRAIN'] = brain
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
        parameters['diameter'] = 2.0 * parameters['radius']
    except KeyError:
        pass
    params['PARAMETERS'] = parameters

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
    folder = params['MRI_DATA_CASE']
    if folder != '':
        if os.path.isdir(folder) == True:
            tmp = os.path.join(folder, 'fiducials.csv')
            if os.path.isfile(tmp) != True:
                print('* ERROR:', folder, 'does not contain fiducials.csv.')
                print('Aborting.')
                exit()
        else:
            print('* ERROR:', folder, 'does not exist.')
            print('Aborting.')
            exit()

    print('Done.')

def calc_delta_time(params, material, parameters):
    RHO = material['rho']
    C = material['c']
    LAMBDA = material['lambda']
    RHO_BLOOD = material['rho_blood']
    C_BLOOD = material['c_blood']
    OMEGA = material['omega']
    H = parameters['h']
    GRIDSIZE = params['GRIDSIZE']
    SPACE_DIM = params['SPACE_DIM']
    # Pennes Bioheat Equation.
    tmp = 0
    for dim in range(0, SPACE_DIM):
        tmp += (2.0/(GRIDSIZE[dim]*GRIDSIZE[dim])) * (LAMBDA/(RHO*C))
    # Inner nodes.
    DELTA_TIME = tmp + ((RHO_BLOOD*C_BLOOD)/(RHO*C)) * OMEGA
    DELTA_TIME = 1.0/DELTA_TIME
    # Border with convection.
    DELTA_TIME_CONV = tmp + 2.0*(1.0/GRIDSIZE[SPACE_DIM-1]) \
                                  * (H/(RHO*C))
    DELTA_TIME_CONV += ((RHO_BLOOD*C_BLOOD)/(RHO*C)) * OMEGA
    DELTA_TIME_CONV = 1.0/DELTA_TIME_CONV

    return DELTA_TIME, DELTA_TIME_CONV

def calc_variables(params):
    print('Calculating variables.')

    # Calculate gridsize in each dimension.
    GRIDSIZE = []
    for dim in range(0, params['SPACE_DIM']):
        GRIDSIZE.append((params['COORD_NODE_LAST'][dim] \
                         - params['COORD_NODE_FIRST'][dim])
                        / (params['N_NODES'][dim]-1))
    params['GRIDSIZE'] = GRIDSIZE
    # Calculate delta time.
    if params['N_TIMESTEPS'] < 1:
        print('* WARNING: N_TIMESTEPS not specified.')
        print('  Calculate N_TIMESTEPS from stability criterion.')
        DELTA_TIME_BRAIN, \
        DELTA_TIME_BRAIN_CONV = calc_delta_time(params, params['BRAIN'],
                                                params['PARAMETERS'])
        DELTA_TIME_TUMOR, \
        DELTA_TIME_TUMOR_CONV = calc_delta_time(params, params['TUMOR'],
                                                params['PARAMETERS'])
        # Get minimum for calculation of timesteps.
        DELTA_TIME_MIN = min((DELTA_TIME_BRAIN, DELTA_TIME_BRAIN_CONV,
                              DELTA_TIME_TUMOR, DELTA_TIME_TUMOR_CONV))
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
                         - params['COORD_NODE_FIRST'][0]) / 2.0)
    TUMOR_CENTER.append((params['COORD_NODE_LAST'][1] \
                         - params['COORD_NODE_FIRST'][1]) / 2.0)
    TUMOR_CENTER.append(params['COORD_NODE_LAST'][2]
                        - params['PARAMETERS']['depth'])
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

    DELTA_TIME_BRAIN, \
    DELTA_TIME_BRAIN_CONV = calc_delta_time(params, params['BRAIN'],
                                            params['PARAMETERS'])
    DELTA_TIME_TUMOR, \
    DELTA_TIME_TUMOR_CONV = calc_delta_time(params, params['TUMOR'],
                                            params['PARAMETERS'])
    # Get minimum for calculation of timesteps.
    DELTA_TIME_MIN = min((DELTA_TIME_BRAIN, DELTA_TIME_BRAIN_CONV,
                          DELTA_TIME_TUMOR, DELTA_TIME_TUMOR_CONV))

    DELTA_TIME = params['DELTA_TIME']
    # Abort simulation if stability is not fulfilled.
    if DELTA_TIME > DELTA_TIME_BRAIN:
        print('* ERROR: Stability not fulfilled in healthy brain region.')
        print('  DELTA_TIME = {0}, but has to be DELTA_TIME < {1}.'.format(DELTA_TIME,
                                                                           DELTA_TIME_BRAIN))
        print('Aborting.')
        exit()
    # Abort simulation if stability is not fulfilled.
    if DELTA_TIME > DELTA_TIME_TUMOR:
        print('* ERROR: Stability not fulfilled in tumor region.')
        print('  DELTA_TIME = {0}, but has to be DELTA_TIME < {1}.'.format(DELTA_TIME,
                                                                           DELTA_TIME_TUMOR))
        print('Aborting.')
        exit()
    # Abort simulation if stability is not fulfilled.
    if DELTA_TIME > DELTA_TIME_BRAIN_CONV:
        print('* ERROR: Stability not fulfilled in healty brain region at \
              border with convection.')
        print('  DELTA_TIME = {0}, but has to be DELTA_TIME < {1}.'.format(DELTA_TIME,
                                                                           DELTA_TIME_BRAIN_CONV))
        print('Aborting.')
        exit()
    # Abort simulation if stability is not fulfilled.
    if DELTA_TIME > DELTA_TIME_TUMOR_CONV:
        print('* ERROR: Stability not fulfilled in tumor region at border \
              with convection.')
        print('  DELTA_TIME = {0}, but has to be DELTA_TIME < {1}.'.format(DELTA_TIME,
                                                                           DELTA_TIME_TUMOR_CONV))
        print('Aborting.')
        exit()

    print('Done.')

def create_region_array(params, nc_file, BRAIN_VALUE, TUMOR_VALUE,
                        NAME_VARIABLE):
    RADIUS = params['PARAMETERS']['diameter']/2
    # Get file/grid dimensions.
    dim0 = params['N_NODES'][0]
    dim1 = params['N_NODES'][1]
    dim2 = params['N_NODES'][2]
    # Get tumor center location.
    TUMOR_CENTER = params['TUMOR_CENTER']
    # Resize temperature array.
    num_elem = dim0 * dim1 * dim2
    values_array = BRAIN_VALUE \
                   * np.ones(num_elem, dtype=int).reshape(dim2, dim1, dim0)
    # Iterate through temperature array.
    for elem_z in range(0, values_array.shape[0]):
        for elem_y in range(0, values_array.shape[1]):
            for elem_x in range(0, values_array.shape[2]):
                # Calculate location of current node.
                x = elem_x * params['GRIDSIZE'][0]
                y = elem_y * params['GRIDSIZE'][1]
                z = elem_z * params['GRIDSIZE'][2]
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
    filepath = 'region.nc'
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

def create_init_array(params, nc_file, region, BRAIN_VALUE, TUMOR_VALUE,
                      NAME_VARIABLE):
    # Get file/grid dimensions.
    dim0 = params['N_NODES'][0]
    dim1 = params['N_NODES'][1]
    dim2 = params['N_NODES'][2]
    # Resize temperature array.
    num_elem = dim0 * dim1 * dim2
    values_array = BRAIN_VALUE * np.ones(num_elem).reshape(dim2, dim1, dim0)
    # Iterate through temperature array.
    for elem_z in range(0, values_array.shape[0]):
        for elem_y in range(0, values_array.shape[1]):
            for elem_x in range(0, values_array.shape[2]):
                # Check if current point is inside tumor.
                # If yes, set value to tumor specific value.
                if region[elem_z, elem_y, elem_x] == 1:
                    values_array[elem_z, elem_y, elem_x] = TUMOR_VALUE
    # Write NumPy array to netCDF file.
    write_values_to_file(nc_file, values_array, NAME_VARIABLE)

def create_surface_array(params, nc_file, BRAIN_VALUE, TUMOR_VALUE,
                         NAME_VARIABLE):
    RADIUS = (params['PARAMETERS']['diameter'] \
              * params['PARAMETERS']['hole_factor'])/2
    # Get file/grid dimensions.
    dim0 = params['N_NODES'][0]
    dim1 = params['N_NODES'][1]
    dim2 = params['N_NODES'][2]
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
            x = elem_x * params['GRIDSIZE'][0]
            y = elem_y * params['GRIDSIZE'][1]
            # Calculate distance (squared) to tumor center.
            distance = (x - TUMOR_CENTER[0]) * (x - TUMOR_CENTER[0])
            distance += (y - TUMOR_CENTER[1]) * (y - TUMOR_CENTER[1])
            # Check if current point is inside tumor.
            # If yes, set value to tumor specific value
            if distance <= RADIUS*RADIUS:
                values_array[dim2-1, elem_y, elem_x] = TUMOR_VALUE
    # Create netCDF variable.
    nNodes = []
    nNodes.append('time')
    for dim in range(len(values_array.shape), 0, -1):
        nNodes.append('nNodes_' + str(dim-1))
    init_values = nc_file.createVariable(NAME_VARIABLE, 'i', nNodes)
    # Write NumPy Array to file.
    init_values[0,] = values_array

def create_surface_from_mri(params, nc_file, BRAIN_VALUE, TUMOR_VALUE,
                            NAME_VARIABLE):
    filepath = params['MRI_DATA_CASE']
    iop = read_intra_op_points(filepath)
    t = read_tumor_point(filepath)
    iop, t = rotate_points(iop, t)
    iop, t = move_points(iop, t, t)

    path = get_interpolated_path(iop)
    #path = get_path(iop)

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

    # Get file/grid dimensions.
    dim0 = params['N_NODES'][0]
    dim1 = params['N_NODES'][1]
    dim2 = params['N_NODES'][2]
    # Resize array.
    num_elem = dim0 * dim1 * dim2
    values_array = BRAIN_VALUE \
                   * np.ones(num_elem, dtype=int).reshape(dim2, dim1, dim0)
    # Iterate through array.
    for elem_y in range(0, dim1):
        for elem_x in range(0, dim0):
            # Calculate location of current node.
            x = elem_x * params['GRIDSIZE'][0]
            y = elem_y * params['GRIDSIZE'][1]
            # Transform current node to tumor center system.
            x_trans = (x - TUMOR_CENTER[0])*1000
            y_trans = (y - TUMOR_CENTER[1])*1000
            # Check if current point is inside open skill
            # If yes, set value to tumor specific value.
            if path.contains_point((x_trans,y_trans)) == True:
                values_array[dim2-1, elem_y, elem_x] = TUMOR_VALUE
    # Create netCDF variable.
    nNodes = []
    nNodes.append('time')
    for dim in range(len(values_array.shape), 0, -1):
        nNodes.append('nNodes_' + str(dim-1))
    init_values = nc_file.createVariable(NAME_VARIABLE, 'i', nNodes)
    # Write NumPy Array to file.
    init_values[0,] = values_array

def create_init_file(params):
    filepath = params['NAME_INITFILE'] + '.nc'
    SPACE_DIM = params['SPACE_DIM']
    print('Creating {0}.'.format(filepath))

    # Delete old init file.
    if os.path.isfile(filepath) == True:
        os.remove(filepath)

    # Check if region file exists.
    if os.path.isfile('region.nc') == False:
        print('* ERROR: region.nc does not exist.')
        print('Aborting.')
        exit()

    # Open region file.
    nc_file = nc.Dataset('region.nc')
    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size
    time = nc_file.dimensions['time'].size
    nc_var = nc_file.variables['region']
    region = np.zeros((dim2, dim1, dim0), dtype=int)
    region[:,:,:] = nc_var[(time-1):time,:,:,:]

    nc_file.close()

    # Create init file.
    nc_file = nc.Dataset(filepath, 'w', format='NETCDF3_CLASSIC')
    time = nc_file.createDimension('time')
    for dim in range(0, SPACE_DIM):
        nNodes = nc_file.createDimension('nNodes_' + str(dim),
                                         params['N_NODES'][dim])

    brain = params['BRAIN']
    tumor = params['TUMOR']
    names = {'rho': 'rho', 'c': 'c', 'lambda': 'lambda',
             'rho_blood': 'rho_blood', 'c_blood': 'c_blood', 'omega': 'omega',
             't_blood': 'T_blood', 'q': 'q', 't': 'T'}
    for key, value in brain.items():
        create_init_array(params, nc_file, region, brain[key], tumor[key],
                          names[key])
    if params['MRI_DATA_CASE'] != '':
        create_surface_from_mri(params, nc_file, 0, 1, 'surface')
    else:
        create_surface_array(params, nc_file, 0, 1, 'surface')

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

def start_simulation(params, run_script):
    # Check if run sript exists.
    if os.path.isfile(run_script) == False:
        print('* ERROR:', run_script, 'does not exist.')
        print('Aborting.')
        exit()

    # Get time of newest netCDF file BEFORE simulation.
    files_nc_before = glob.glob('./' + params['NAME_EXECUTABLE'] + '*.nc')
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

    # Get time of newest netCDF file AFTER simulation.
    files_nc_after = glob.glob(params['NAME_EXECUTABLE'] + '*.nc')
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
            params['NAME_CONFIGFILE'] = sys.argv[1]
            run_script = 'RUN_HELPER.sh'
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
    create_region_file(params)
    create_init_file(params)
    set_environment_variables(params)
    start_simulation(params, run_script)
    if params['NAME_RESULTFILE'] != '' and params['SPACE_DIM'] == 3:
        plot_surface(params['NAME_RESULTFILE'], params)
        mean_surface_temperature(params['NAME_RESULTFILE'])
        mean_tumor_temperature(params['NAME_RESULTFILE'])

if __name__ == '__main__':
    main()
