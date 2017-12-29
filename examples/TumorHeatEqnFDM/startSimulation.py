import configparser
import netCDF4 as nc
import numpy as np
import os
import subprocess
import sys

params = {
    'NAME_CONFIGFILE' : ''
}

def parse_config_file():
    global params
    print('Parsing {0}.'.format(params['NAME_CONFIGFILE']))

    # Create configparser and open file.
    config = configparser.ConfigParser()
    config.read(params['NAME_CONFIGFILE'])
    # Get values from section 'Dimension'.
    params['SPACE_DIM'] = config['Dimension'].getint('SPACE_DIM', fallback=1)
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
    params['N_TIMESTEPS'] = config['Time'].getint('N_TIMESTEPS', fallback=1)
    # Get values from section 'Output'.
    params['N_SNAPSHOTS'] = config['Output'].getint('N_SNAPSHOTS')
    # Get values from section 'Input'.
    params['NAME_INITFILE'] = config['Input'].get('NAME_INITFILE', fallback='init')
    params['USE_INITFILE'] = config['Input'].getboolean('USE_INITFILE', fallback=False)
    params['CREATE_INITFILE'] = config['Input'].getboolean('CREATE_INITFILE', fallback=False)
    NAME_VARIABLES = config['Input'].get('NAME_VARIABLES', fallback='U')
    NAME_VARIABLES = list(NAME_VARIABLES.split())
    params['NAME_VARIABLES'] = NAME_VARIABLES
    # Get values from section 'Parameters'.
    params['T_INIT'] = config['Parameters'].getfloat('T_I')
    params['T_TUMOR'] = config['Parameters'].getfloat('T_TUMOR')
    params['DIAMETER'] = config['Parameters'].getfloat('DIAMETER')
    params['DEPTH'] = config['Parameters'].getfloat('DEPTH')
    params['RHO'] = config['Parameters'].getfloat('RHO')
    params['C'] = config['Parameters'].getfloat('C')
    params['K'] = config['Parameters'].getfloat('K')
    params['RHO_B'] = config['Parameters'].getfloat('RHO_B')
    params['C_PB'] = config['Parameters'].getfloat('C_PB')
    params['OMEGA_B_BRAIN'] = config['Parameters'].getfloat('OMEGA_B_BRAIN')
    params['OMEGA_B_TUMOR'] = config['Parameters'].getfloat('OMEGA_B_TUMOR')
    params['H'] = config['Parameters'].getfloat('H')

    print('Done.')

def calc_variables():
    global params
    print('Calculating variables.')

    # Calculate gridsize in each dimension.
    GRIDSIZE = []
    for dim in range(0, params['SPACE_DIM']):
        GRIDSIZE.append((params['COORD_NODE_LAST'][dim] - params['COORD_NODE_FIRST'][dim])/(params['N_NODES'][dim]-1))
    params['GRIDSIZE'] = GRIDSIZE
    # Calculate delta time.
    DELTA_TIME = (params['END_TIME'] - params['START_TIME'])/params['N_TIMESTEPS']
    params['DELTA_TIME'] = DELTA_TIME

    print('Done.')

def check_variables():
    global params
    print('Checking variables.')

    # Check if dimension makes sense and
    # some functions and variables only work for dimension 1, 2 or 3.
    SPACE_DIM = params['SPACE_DIM']
    if SPACE_DIM < 1 or SPACE_DIM > 3:
        print('SPACE_DIM is {0}.'.format(SPACE_DIM))
        print('SPACE_DIM must be 1, 2 or 3.')
        print('Aborting.')
        exit()
    # Check if there are enough coordinates for first node.
    DIM_COORD_NODE_FIRST = len(params['COORD_NODE_FIRST'])
    if DIM_COORD_NODE_FIRST != SPACE_DIM:
        print('Dimension of COORD_NODE_FIRST has to be {0}.'.format(SPACE_DIM))
        print('Dimension of COORD_NODE_FIRST is {0}.'.format(DIM_COORD_NODE_FIRST))
        print('Aborting.')
        exit()
    # Check if there are enough coordinates for last node.
    DIM_COORD_NODE_LAST = len(params['COORD_NODE_LAST'])
    if DIM_COORD_NODE_LAST != SPACE_DIM:
        print('Dimension of COORD_NODE_LAST has to be {0}.'.format(SPACE_DIM))
        print('Dimension of COORD_NODE_LAST is {0}.'.format(DIM_COORD_NODE_LAST))
        print('Aborting.')
        exit()
    # Check if there are enough number of nodes.
    DIM_N_NODES = len(params['N_NODES'])
    if DIM_N_NODES != SPACE_DIM:
        print('Dimension of N_NODES has to be {0}.'.format(SPACE_DIM))
        print('Dimension of N_NODES is {0}.'.format(DIM_N_NODES))
        print('Aborting.')
        exit()
    # Check if END_TIME is after START_TIME.
    START_TIME = params['START_TIME']
    END_TIME = params['END_TIME']
    if END_TIME < START_TIME:
        print('END_TIME is smaller than START_TIME.')
        print('END_TIME must be greater than START_TIME.')
        print('Aborting.')
        exit()
    # Check if number of snapshots is possible.
    if params['N_SNAPSHOTS'] > params['N_TIMESTEPS']:
        print('WARNING: N_SNAPSHOTS was bigger than N_TIMESTEPS.')
        params['N_SNAPSHOTS'] = params['N_TIMESTEPS']
        print('N_SNAPSHOTS was set to N_TIMESTEPS.')
    # Check if combinations of USE_INITFILE and CREATE_INITFILE makes sense.
    if params['USE_INITFILE'] == True and params['CREATE_INITFILE'] == False:
        if os.path.isfile(params['NAME_INITFILE'] + '.nc') == False:
            print('USE_INITFILE = True and CREATE_INITFILE = False, but', params['NAME_INITFILE'] + '.nc', 'does not exist.')
            print('Aborting.')
            exit()
    if params['USE_INITFILE'] == False and params['CREATE_INITFILE'] == True:
        print('WARNING: CREATE_INITFILE = True, but USE_INITFILE = False.')
    else:
        pass
    # Check if executable exists.
    NAME_EXECUTABLE = os.path.basename(os.getcwd()) + str(params['SPACE_DIM']) + 'D'
    if os.path.isfile(NAME_EXECUTABLE) == False:
        print(NAME_EXECUTABLE, 'does not exist.')
        print('Aborting.')
        exit()
    params['NAME_EXECUTABLE'] = NAME_EXECUTABLE

    print('Done.')

def check_stability():
    global params
    print('Checking stability.')

    RHO = params['RHO']
    C = params['C']
    K = params['K']
    RHO_B = params['RHO_B']
    C_PB = params['C_PB']
    OMEGA_B_BRAIN = params['OMEGA_B_BRAIN']
    OMEGA_B_TUMOR = params['OMEGA_B_TUMOR']
    H = params['H']
    GRIDSIZE = params['GRIDSIZE']
    SPACE_DIM = params['SPACE_DIM']
    DELTA_TIME = params['DELTA_TIME']
    # Pennes Bioheat Equation.
    tmp = 0
    for dim in range(0, SPACE_DIM):
        tmp += (2.0/(GRIDSIZE[dim]*GRIDSIZE[dim])) * (K/(RHO*C))
    # Healthy brain region in inner nodes.
    DELTA_TIME_BRAIN = tmp + ((RHO_B*C_PB)/(RHO*C)) * OMEGA_B_BRAIN
    DELTA_TIME_BRAIN = 1.0/DELTA_TIME_BRAIN
    # Tumor region in inner nodes.
    DELTA_TIME_TUMOR = tmp + ((RHO_B*C_PB)/(RHO*C)) * OMEGA_B_TUMOR
    DELTA_TIME_TUMOR = 1.0/DELTA_TIME_TUMOR
    # Healthy brain region at border with convection.
    DELTA_TIME_BRAIN_CONV = tmp + 2.0*(1.0/GRIDSIZE[SPACE_DIM-1])*(H/(RHO*C))
    DELTA_TIME_BRAIN_CONV += ((RHO_B*C_PB)/(RHO*C)) * OMEGA_B_BRAIN
    DELTA_TIME_BRAIN_CONV = 1.0/DELTA_TIME_BRAIN_CONV
    # Tumor brain region at border with convection.
    # Will probably not be the case, but test it anyway.
    DELTA_TIME_TUMOR_CONV = tmp + 2.0*(1.0/GRIDSIZE[SPACE_DIM-1])*(H/(RHO*C))
    DELTA_TIME_TUMOR_CONV += ((RHO_B*C_PB)/(RHO*C)) * OMEGA_B_TUMOR
    DELTA_TIME_TUMOR_CONV = 1.0/DELTA_TIME_TUMOR_CONV

    # Abort simulation if stability is not fulfilled.
    if DELTA_TIME > DELTA_TIME_BRAIN:
        print('Stability not fulfilled in healthy brain region.')
        print('DELTA_TIME = {0}, but has to be DELTA_TIME < {1}.'.format(DELTA_TIME, DELTA_TIME_BRAIN))
        print('Aborting.')
        exit()
    # Abort simulation if stability is not fulfilled.
    if DELTA_TIME > DELTA_TIME_TUMOR:
        print('Stability not fulfilled in tumor region.')
        print('DELTA_TIME = {0}, but has to be DELTA_TIME < {1}.'.format(DELTA_TIME, DELTA_TIME_TUMOR))
        print('Aborting.')
        exit()
    # Abort simulation if stability is not fulfilled.
    if DELTA_TIME > DELTA_TIME_BRAIN_CONV:
        print('Stability not fulfilled in healty brain region at border with convection.')
        print('DELTA_TIME = {0}, but has to be DELTA_TIME < {1}.'.format(DELTA_TIME, DELTA_TIME_BRAIN_CONV))
        print('Aborting.')
        exit()
    # Abort simulation if stability is not fulfilled.
    if DELTA_TIME > DELTA_TIME_TUMOR_CONV:
        print('Stability not fulfilled in tumor region at border with convection.')
        print('DELTA_TIME = {0}, but has to be DELTA_TIME < {1}.'.format(DELTA_TIME, DELTA_TIME_TUMOR_CONV))
        print('Aborting.')
        exit()

    print('Done.')

def create_init_file():
    global params
    SPACE_DIM = params['SPACE_DIM']
    NAME_INITFILE = params['NAME_INITFILE']
    NAME_VARIABLES = params['NAME_VARIABLES']
    print('Creating {0}.nc.'.format(NAME_INITFILE))

    nc_file = nc.Dataset(NAME_INITFILE + '.nc', 'w', format='NETCDF3_CLASSIC')
    time = nc_file.createDimension('time')
    nNodes_0 = nc_file.createDimension('nNodes_0', params['N_NODES'][0])

    if SPACE_DIM == 1:
        create_value_array_1D(nc_file, params['T_INIT'], params['T_TUMOR'], NAME_VARIABLES[0])
    elif SPACE_DIM == 2:
        nNodes_1 = nc_file.createDimension('nNodes_1', params['N_NODES'][1])
        create_value_array_2D(nc_file, params['T_INIT'], params['T_TUMOR'], NAME_VARIABLES[0])
    else:
        nNodes_1 = nc_file.createDimension('nNodes_1', params['N_NODES'][1])
        nNodes_2 = nc_file.createDimension('nNodes_2', params['N_NODES'][2])
        create_value_array_3D(nc_file, params['T_INIT'], params['T_TUMOR'], NAME_VARIABLES[0])

    nc_file.close()
    print('Done.')

def create_value_array_1D(nc_file, BRAIN_VALUE, TUMOR_VALUE, NAME_VARIABLE):
    RADIUS = params['DIAMETER']/2
    TUMOR_CENTER = []
    # Get file/grid dimensions.
    dim0 = params['N_NODES'][0]
    # Resize temperature array.
    num_elem = dim0
    values_array = BRAIN_VALUE * np.ones(num_elem).reshape(dim0)
    # Calculate location of tumor center.
    TUMOR_CENTER.append(params['COORD_NODE_LAST'][0] - params['DEPTH'])
    # Iterate through array.
    for elem in range(0, values_array.shape[0]):
        # Calculate location of current node.
        x = elem * params['GRIDSIZE'][0]
        # Calculate distance (squared) to tumor center.
        distance = (x - TUMOR_CENTER[0]) * (x - TUMOR_CENTER[0])
        # Check if current point is inside tumor.
        # If yes, set value to tumor specific value
        if distance <= RADIUS*RADIUS:
            values_array[elem] = TUMOR_VALUE
    # Write NumPy array to netCDF file.
    write_values_to_file_1D(nc_file, values_array, NAME_VARIABLE)

def write_values_to_file_1D(nc_file, values_array, NAME_VARIABLE):
    # Create netCDF variable.
    init_values = nc_file.createVariable(NAME_VARIABLE, 'f8', ('time', 'nNodes_0'))
    # Write NumPy Array to file.
    init_values[0,:] = values_array

def create_value_array_2D(nc_file, BRAIN_VALUE, TUMOR_VALUE, NAME_VARIABLE):
    RADIUS = params['DIAMETER']/2
    TUMOR_CENTER = []
    # Get file/grid dimensions.
    dim0 = params['N_NODES'][0]
    dim1 = params['N_NODES'][1]
    # Resize temperature array.
    num_elem = dim0 * dim1
    values_array = BRAIN_VALUE * np.ones(num_elem).reshape(dim1, dim0)
    # Calculate location of tumor center.
    TUMOR_CENTER.append(params['COORD_NODE_LAST'][0]/2.0)
    TUMOR_CENTER.append(params['COORD_NODE_LAST'][1] - params['DEPTH'])
    # Iterate through temperature array.
    for elem_y in range(0, values_array.shape[0]):
        for elem_x in range(0, values_array.shape[1]):
            # Calculate location of current node.
            x = elem_x * params['GRIDSIZE'][0]
            y = elem_y * params['GRIDSIZE'][1]
            # Calculate distance (squared) to tumor center.
            distance = (x - TUMOR_CENTER[0]) * (x - TUMOR_CENTER[0])
            distance += (y - TUMOR_CENTER[1]) * (y - TUMOR_CENTER[1])
            # Check if current point is inside tumor.
            # If yes, set value to tumor specific value
            if distance <= RADIUS*RADIUS:
                values_array[elem_y, elem_x] = TUMOR_VALUE
    # Write NumPy array to netCDF file.
    write_values_to_file_2D(nc_file, values_array, NAME_VARIABLE)

def write_values_to_file_2D(nc_file, values_array, NAME_VARIABLE):
    # Create netCDF variable.
    init_values = nc_file.createVariable(NAME_VARIABLE, 'f8', ('time', 'nNodes_1', 'nNodes_0'))
    # Write NumPy Array to file.
    init_values[0,:,:] = values_array

def create_value_array_3D(nc_file, BRAIN_VALUE, TUMOR_VALUE, NAME_VARIABLE):
    RADIUS = params['DIAMETER']/2
    TUMOR_CENTER = []
    # Get file/grid dimensions.
    dim0 = params['N_NODES'][0]
    dim1 = params['N_NODES'][1]
    dim2 = params['N_NODES'][2]
    # Resize temperature array.
    num_elem = dim0 * dim1 * dim2
    values_array = BRAIN_VALUE * np.ones(num_elem).reshape(dim2, dim1, dim0)
    # Calculate location of tumor center.
    TUMOR_CENTER.append(params['COORD_NODE_LAST'][0]/2.0)
    TUMOR_CENTER.append(params['COORD_NODE_LAST'][1]/2.0)
    TUMOR_CENTER.append(params['COORD_NODE_LAST'][2] - params['DEPTH'])
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
    # Write NumPy array to netCDF file.
    write_values_to_file_3D(nc_file, values_array, NAME_VARIABLE)

def write_values_to_file_3D(nc_file, values_array, NAME_VARIABLE):
    # Create netCDF variable.
    init_values = nc_file.createVariable(NAME_VARIABLE, 'f8', ('time', 'nNodes_2', 'nNodes_1', 'nNodes_0'))
    # Write NumPy Array to file.
    init_values[0,:,:,:] = values_array

def set_environment_variables():
    global params
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
    # Check if init file should be used and if it exists.
    if params['USE_INITFILE'] == True:
        if os.path.isfile(params['NAME_INITFILE'] + '.nc') == True:
            os.putenv('SCAFESRUN_NAME_INITFILE', params['NAME_INITFILE'])
        else:
            print('USE_INITFILE = True, but', params['NAME_INITFILE'] + '.nc', 'does not exist.')
            print('Aborting.')
            exit()

    print('Done.')

def start_simulation():
    global params

    # Set name of folder as executable with dimension as suffix
    print('Starting {0}.'.format(params['NAME_EXECUTABLE']))
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
    check_variables()
    calc_variables()
    check_stability()
    if params['CREATE_INITFILE'] == True:
        create_init_file()
    set_environment_variables()
    start_simulation()

if __name__ == '__main__':
    main()
