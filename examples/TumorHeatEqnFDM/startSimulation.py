import configparser
import math
import netCDF4 as nc
import numpy as np
import os
import subprocess
import sys

params = {
    'NAME_CONFIGFILE' : '',
    'SPACE_DIM' : 0,
    'GRIDSIZE_GLOBAL' : 0.0,
    'COORD_NODE_FIRST_DIM1' : 0.0,
    'COORD_NODE_FIRST_DIM2' : 0.0,
    'COORD_NODE_FIRST_DIM3' : 0.0,
    'COORD_NODE_LAST_DIM1' : 0.0,
    'COORD_NODE_LAST_DIM2' : 0.0,
    'COORD_NODE_LAST_DIM3' : 0.0,
    'START_TIME' : 0,
    'END_TIME' : 0,
    'DELTA_TIME' : 0.0,
    'N_SNAPSHOTS': 0,
    'N_TIMESTEPS' : 0,
    'N_NODES_DIM1' : 0,
    'N_NODES_DIM2' : 0,
    'N_NODES_DIM3' : 0,
    'NAME_INITFILE' : '',
    'USE_INITFILE' : False,
    'CREATE_INITFILE' : False,
}

def parse_config_file():
    global params
    print('Parsing {0}.'.format(params['NAME_CONFIGFILE']))

    # Create configparser and open file.
    config = configparser.ConfigParser()
    config.read(params['NAME_CONFIGFILE'])
    # Get values from section 'Dimension'.
    params['SPACE_DIM'] = config['Dimension'].getint('SPACE_DIM')
    # Check if dimension makes sense and
    # some functions and variables only work for dimension 1, 2 or 3.
    SPACE_DIM = params['SPACE_DIM']
    if SPACE_DIM < 1 or SPACE_DIM > 3:
        print('SPACE_DIM is {0}.'.format(SPACE_DIM))
        print('SPACE_DIM must be 1, 2 or 3.')
        print('Aborting.')
        exit()
    # Get values from section 'Geometry'.
    params['GRIDSIZE_GLOBAL'] = config['Geometry'].getfloat('GRIDSIZE_GLOBAL')
    params['COORD_NODE_FIRST_DIM1'] = config['Geometry'].getfloat('COORD_NODE_FIRST_DIM1')
    params['COORD_NODE_FIRST_DIM2'] = config['Geometry'].getfloat('COORD_NODE_FIRST_DIM2')
    params['COORD_NODE_FIRST_DIM3'] = config['Geometry'].getfloat('COORD_NODE_FIRST_DIM3')
    params['COORD_NODE_LAST_DIM1'] = config['Geometry'].getfloat('COORD_NODE_LAST_DIM1')
    params['COORD_NODE_LAST_DIM2'] = config['Geometry'].getfloat('COORD_NODE_LAST_DIM2')
    params['COORD_NODE_LAST_DIM3'] = config['Geometry'].getfloat('COORD_NODE_LAST_DIM3')
    # Get values from section 'Time'.
    params['START_TIME'] = config['Time'].getint('START_TIME')
    params['END_TIME']= config['Time'].getint('END_TIME')
    params['DELTA_TIME'] = config['Time'].getfloat('DELTA_TIME')
    # Get values from section 'Output'.
    params['N_SNAPSHOTS'] = config['Output'].getint('N_SNAPSHOTS')
    # Get values from section 'Input'.
    params['NAME_INITFILE'] = config['Input'].get('NAME_INITFILE', fallback='init')
    params['USE_INITFILE'] = config['Input'].getboolean('USE_INITFILE', fallback=False)
    params['CREATE_INITFILE'] = config['Input'].getboolean('CREATE_INITFILE', fallback=False)
    # Get values from section 'Parameters'.
    params['T_INIT'] = config['Parameters'].getfloat('T_I', fallback=0.0)
    params['T_TUMOR'] = config['Parameters'].getfloat('T_TUMOR', fallback=0.0)
    params['DIAMETER'] = config['Parameters'].getfloat('DIAMETER', fallback=0.0)
    params['DEPTH'] = config['Parameters'].getfloat('DEPTH', fallback=0.0)

    print('Done.')

def calc_variables():
    global params
    print('Calculating variables.')

    # Calculate number of nodes in each dimension.
    N_NODES_DIM1 = (params['COORD_NODE_LAST_DIM1'] - params['COORD_NODE_FIRST_DIM1'])/params['GRIDSIZE_GLOBAL']
    params['N_NODES_DIM1'] = int(math.ceil(N_NODES_DIM1)) + 1
    N_NODES_DIM2 = (params['COORD_NODE_LAST_DIM2'] - params['COORD_NODE_FIRST_DIM2'])/params['GRIDSIZE_GLOBAL']
    params['N_NODES_DIM2'] = int(math.ceil(N_NODES_DIM2)) + 1
    N_NODES_DIM3 = (params['COORD_NODE_LAST_DIM3'] - params['COORD_NODE_FIRST_DIM3'])/params['GRIDSIZE_GLOBAL']
    params['N_NODES_DIM3'] = int(math.ceil(N_NODES_DIM3)) + 1
    # Calculate number of timesteps.
    N_TIMESTEPS = (params['END_TIME'] - params['START_TIME'])/params['DELTA_TIME']
    params['N_TIMESTEPS'] = int(math.ceil(N_TIMESTEPS))
    # Check if number of snapshots is possible.
    if params['N_SNAPSHOTS'] > params['N_TIMESTEPS']:
        print('WARNING: N_SNAPSHOTS was bigger than N_TIMESTEPS.')
        params['N_SNAPSHOTS'] = params['N_TIMESTEPS']
        print('N_SNAPSHOTS was set to N_TIMESTEPS')

    print('Done.')

def create_init_file():
    global params
    global NAME_INITFILE
    global T_INIT
    global T
    SPACE_DIM = params['SPACE_DIM']
    NAME_INITFILE = params['NAME_INITFILE']
    T_INIT = params['T_INIT']
    print('Creating {0}.nc.'.format(NAME_INITFILE))

    if SPACE_DIM == 1:
        create_temperature_array_1D()
        create_init_file_1D()
    elif SPACE_DIM == 2:
        create_temperature_array_2D()
        create_init_file_2D()
    else:
        create_temperature_array_3D()
        create_init_file_3D()

    print('Done.')

def create_temperature_array_1D():
    global T
    RADIUS = params['DIAMETER']/2
    T_TUMOR = params['T_TUMOR']
    # Get file/grid dimensions.
    dim0 = params['N_NODES_DIM1']
    # Resize temperature array.
    num_elem = dim0
    T = T_INIT * np.ones(num_elem).reshape(dim0)
    # Calculate location of tumor center.
    TUMOR_CENTER = params['COORD_NODE_LAST_DIM1'] - params['DEPTH']
    # Iterate through temperature array.
    for elem in range(0, T.shape[0]):
        # Calculate location of current node.
        x = elem * params['GRIDSIZE_GLOBAL']
        # Calculate distance (squared) to tumor center.
        distance = (x - TUMOR_CENTER) * (x - TUMOR_CENTER)
        # Check if current point is inside Tumor.
        # If yes, set temperature for this point to tumor temperature.
        if distance <= RADIUS*RADIUS:
            T[elem] = T_TUMOR

def create_init_file_1D():
    # Get file/grid dimensions.
    dim0 = params['N_NODES_DIM1']
    # Create netCDF file.
    nc_file = nc.Dataset(NAME_INITFILE + '.nc', 'w', format='NETCDF3_CLASSIC')
    nNodes_0 = nc_file.createDimension('nNodes_0', dim0)
    time = nc_file.createDimension('time')
    init_values = nc_file.createVariable('T', 'f8', ('time', 'nNodes_0'))
    # Write NumPy Array to file.
    init_values[0,:] = T
    nc_file.close()

def create_temperature_array_2D():
    global T
    RADIUS = params['DIAMETER']/2
    T_TUMOR = params['T_TUMOR']
    # Get file/grid dimensions.
    dim0 = params['N_NODES_DIM1']
    dim1 = params['N_NODES_DIM2']
    # Resize temperature array.
    num_elem = dim0 * dim1
    T = T_INIT * np.ones(num_elem).reshape(dim1, dim0)
    # Calculate location of tumor center.
    TUMOR_CENTER_DIM1 = params['COORD_NODE_LAST_DIM1']/2.0
    TUMOR_CENTER_DIM2 = params['COORD_NODE_LAST_DIM2'] - params['DEPTH']
    # Iterate through temperature array.
    for elem_y in range(0, T.shape[1]):
        for elem_x in range(0, T.shape[0]):
            # Calculate location of current node.
            x = elem_x * params['GRIDSIZE_GLOBAL']
            y = elem_y * params['GRIDSIZE_GLOBAL']
            # Calculate distance (squared) to tumor center.
            distance = ((x - TUMOR_CENTER_DIM1) * (x - TUMOR_CENTER_DIM1)) \
                        + ((y - TUMOR_CENTER_DIM2) * (y - TUMOR_CENTER_DIM2))
            # Check if current point is inside Tumor.
            # If yes, set temperature for this point to tumor temperature.
            if distance <= RADIUS*RADIUS:
                T[elem_y, elem_x] = T_TUMOR

def create_init_file_2D():
    # Get file/grid dimensions.
    dim0 = params['N_NODES_DIM1']
    dim1 = params['N_NODES_DIM2']
    # Create netCDF file.
    nc_file = nc.Dataset(NAME_INITFILE + '.nc', 'w', format='NETCDF3_CLASSIC')
    nNodes_0 = nc_file.createDimension('nNodes_0', dim0)
    nNodes_1 = nc_file.createDimension('nNodes_1', dim1)
    time = nc_file.createDimension('time')
    init_values = nc_file.createVariable('T', 'f8', ('time', 'nNodes_1', 'nNodes_0'))
    # Write NumPy Array to file.
    init_values[0,:,:] = T
    nc_file.close()

def create_temperature_array_3D():
    global T
    RADIUS = params['DIAMETER']/2
    T_TUMOR = params['T_TUMOR']
    # Get file/grid dimensions.
    dim0 = params['N_NODES_DIM1']
    dim1 = params['N_NODES_DIM2']
    dim2 = params['N_NODES_DIM3']
    # Resize temperature array.
    num_elem = dim0 * dim1 * dim2
    T = T_INIT * np.ones(num_elem).reshape(dim2, dim1, dim0)
    # Calculate location of tumor center.
    TUMOR_CENTER_DIM1 = params['COORD_NODE_LAST_DIM1']/2.0
    TUMOR_CENTER_DIM2 = params['COORD_NODE_LAST_DIM2']/2.0
    TUMOR_CENTER_DIM3 = params['COORD_NODE_LAST_DIM3'] - params['DEPTH']
    # Iterate through temperature array.
    for elem_z in range(0, T.shape[2]):
        for elem_y in range(0, T.shape[1]):
            for elem_x in range(0, T.shape[0]):
                # Calculate location of current node.
                x = elem_x * params['GRIDSIZE_GLOBAL']
                y = elem_y * params['GRIDSIZE_GLOBAL']
                z = elem_z * params['GRIDSIZE_GLOBAL']
                # Calculate distance (squared) to tumor center.
                distance = ((x - TUMOR_CENTER_DIM1) * (x - TUMOR_CENTER_DIM1)) \
                            + ((y - TUMOR_CENTER_DIM2) * (y - TUMOR_CENTER_DIM2)) \
                            + ((z - TUMOR_CENTER_DIM3) * (z - TUMOR_CENTER_DIM3))
                # Check if current point is inside Tumor.
                # If yes, set temperature for this point to tumor temperature.
                if distance <= RADIUS*RADIUS:
                    T[elem_z, elem_y, elem_x] = T_TUMOR

def create_init_file_3D():
    # Get file/grid dimensions.
    dim0 = params['N_NODES_DIM1']
    dim1 = params['N_NODES_DIM2']
    dim2 = params['N_NODES_DIM3']
    # Create netCDF file.
    nc_file = nc.Dataset(NAME_INITFILE + '.nc', 'w', format='NETCDF3_CLASSIC')
    nNodes_0 = nc_file.createDimension('nNodes_0', dim0)
    nNodes_1 = nc_file.createDimension('nNodes_1', dim1)
    nNodes_2 = nc_file.createDimension('nNodes_2', dim2)
    time = nc_file.createDimension('time')
    init_values = nc_file.createVariable('T', 'f8', ('time', 'nNodes_2', 'nNodes_1', 'nNodes_0'))
    # Write NumPy Array to file.
    init_values[0,:,:,:] = T
    nc_file.close()

def set_environment_variables():
    global params
    print('Setting environment variables.')

    # Set all environment that can be set directly.
    os.putenv('SCAFESRUN_N_TIMESTEPS', str(params['N_TIMESTEPS']))
    os.putenv('SCAFESRUN_N_SNAPSHOTS', str(params['N_SNAPSHOTS']))
    os.putenv('SCAFESRUN_START_TIME', str(params['START_TIME']))
    os.putenv('SCAFESRUN_END_TIME', str(params['END_TIME']))
    os.putenv('SCAFESRUN_NAME_CONFIGFILE', params['NAME_CONFIGFILE'])
    os.putenv('SCAFESRUN_SPACE_DIM', str(params['SPACE_DIM']))
    # Build environment variables for coordinates.
    COORD_NODE_FIRST = str(params['COORD_NODE_FIRST_DIM1'])
    COORD_NODE_LAST = str(params['COORD_NODE_LAST_DIM1'])
    N_NODES = str(params['N_NODES_DIM1'])
    if params['SPACE_DIM'] > 1:
        COORD_NODE_FIRST += 'x' + str(params['COORD_NODE_FIRST_DIM2'])
        COORD_NODE_LAST += 'x' + str(params['COORD_NODE_LAST_DIM2'])
        N_NODES += 'x' + str(params['N_NODES_DIM2'])
    if params['SPACE_DIM'] > 2:
        COORD_NODE_FIRST += 'x' + str(params['COORD_NODE_FIRST_DIM3'])
        COORD_NODE_LAST += 'x' + str(params['COORD_NODE_LAST_DIM3'])
        N_NODES += 'x' + str(params['N_NODES_DIM3'])
    # Set environment variables for coordinates.
    os.putenv('SCAFESRUN_COORD_NODE_FIRST', COORD_NODE_FIRST)
    os.putenv('SCAFESRUN_COORD_NODE_LAST', COORD_NODE_LAST)
    os.putenv('SCAFESRUN_N_NODES', N_NODES)
    # Set name of folder as executable with dimension as suffix.
    NAME_EXECUTABLE = os.path.basename(os.getcwd()) + str(params['SPACE_DIM']) + 'D'
    # Check if executable exists and set environment variable if it does.
    # Otherwise abort script.
    if os.path.isfile(NAME_EXECUTABLE) == True:
        os.putenv('SCAFESRUN_NAME_EXECUTABLE', NAME_EXECUTABLE)
    else:
        print(NAME_EXECUTABLE, 'does not exist.')
        print('Aborting.')
        exit()
    # Check if init file should be used and if it exists.
    if params['USE_INITFILE'] == True:
        if os.path.isfile(params['NAME_INITFILE'] + '.nc') == True:
            os.putenv('SCAFESRUN_NAME_INITFILE', params['NAME_INITFILE'])
        else:
            print('USE_INITFILE = True, but', params['NAME_INITFILE'] + '.nc', 'does not exist.')
            print('Aborting.')
            exit()
    elif params['CREATE_INITFILE'] == True:
        print('WARNING: CREATE_INITFILE = True, but USE_INITFILE = False.')
    else:
        pass

    print('Done.')

def start_simulation():
    global params

    # Set name of folder as executable with dimension as suffix
    NAME_EXECUTABLE = os.path.basename(os.getcwd()) + str(params['SPACE_DIM']) + 'D'
    print('Starting {0}.'.format(NAME_EXECUTABLE))
    print()
    # Call bash script to set more environment variables and
    # to start simulation.
    returncode = subprocess.call('./RUN_HELPER.sh')
    # Check if simulation/bash script ran successfully.
    print()
    if returncode == 0:
        print('Done.')
    else:
        print('Simulation returned error code {0}.'.format(returncode))


def main():
    global params
    # Check if path to configfile is provided and
    # if file exists.
    if len(sys.argv) > 1:
        if os.path.isfile(sys.argv[1]) == True:
            params['NAME_CONFIGFILE'] = sys.argv[1]
        else:
            print(sys.argv[1], 'does not exist.')
            print('Usage: python3', sys.argv[0], '<PATH/TO/CONFIGFILE>')
            print('Aborting.')
            exit()
    else:
        print('No command line argument for configfile provided.')
        print('Usage: python3', sys.argv[0], '<PATH/TO/CONFIGFILE>')
        print('Aborting.')
        exit()

    parse_config_file()
    calc_variables()
    if params['CREATE_INITFILE'] == True:
        create_init_file()
    set_environment_variables()
    start_simulation()

if __name__ == '__main__':
    main()
