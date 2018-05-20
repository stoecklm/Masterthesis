import configparser
import time

import netCDF4 as nc
import numpy as np

def temperature_array_from_result(filepath):
    print('Read {}.'.format(filepath))
    nc_file = nc.Dataset(filepath)
    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size

    possible_names = ['T', 'TNewDom', 'TDiff']
    found_name = False
    for name in possible_names:
        try:
            T = nc_file.variables[name]
            found_name = True
            break
        except KeyError:
            pass

    if found_name == False:
        print('* ERROR: No temperature variable found in this file.')
        print('Aborting.')
        exit()

    temp = np.zeros((dim2, dim1, dim0))
    temp[:,:,:] = T[-1,:,:,:]

    nc_file.close()

    return temp

def surface_temperature_array_from_result(filepath):
    temp = temperature_array_from_result(filepath)
    dim2, dim1, dim0 = temp.shape

    temp_surface = np.zeros((dim1, dim0))
    temp_surface[:,:] = temp[-1,:,:]

    return temp_surface

def parse_pymc_from_config_file(params):
    print('Parsing {} for PyMC parameters.'.format(params['NAME_CONFIGFILE_TEMPLATE']))

    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(params['NAME_CONFIGFILE_TEMPLATE'])

    params['ITERATIONS'] = config['PyMC'].getint('ITERATIONS', fallback=5)
    params['BURNS'] = config['PyMC'].getint('BURNS', fallback=1)
    params['T_NORMAL'] = config['PyMC'].getfloat('T_NORMAL', fallback=32.8)
    params['T_TUMOR'] = config['PyMC'].getfloat('T_TUMOR', fallback=30.0)
    params['T_VESSEL'] = config['PyMC'].getfloat('T_VESSEL', fallback=34.5)

    print('Done.')

def create_testcase_name(tested_variables, params):
    case = params['NAME_CONFIGFILE_TEMPLATE'].split('.')[0]
    case = case.split('_')[0]
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(params['NAME_CONFIGFILE_TEMPLATE'])

    n_nodes = config['Geometry'].get('N_NODES')
    current_time = time.strftime("%Y%m%d%H%M%S", time.localtime())
    name = tested_variables + '_' + case + '_' + n_nodes + '_' + current_time

    return name

def create_mcmc_netcdf_file(filepath, dimension):
    print('Save data to {}.'.format(filepath))
    nc_file = nc.Dataset(filepath, 'w', format='NETCDF4')
    nc_file.createDimension('iterations', dimension)
    nc_file.createDimension('scalar', 1)

    return nc_file

def write_ini_file_to_nc_file(nc_file, inifilepath):
    print('Save {} to file.'.format(inifilepath))
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(inifilepath)
    for section in config.sections():
        nc_file.createGroup('ini/' + section)
        items = dict(config.items(section))
        for key, value in items.items():
            var_name = 'ini/' + section + '/' + key
            values = nc_file.createVariable(var_name, str, ('scalar'))
            values[:] = np.array([value], dtype='object')

def save_vector_to_mcmc_file(nc_file, vector_data, vector_name):
    print('Save {} to file.'.format(vector_name))
    values = nc_file.createVariable(vector_name, 'f8', ('iterations'))
    values[:] = vector_data[:]

def close_nc_file(nc_file):
    nc_file.close()

    print('Done.')

if __name__ == '__main__':
    pass
