import configparser
import os
import sys

import netCDF4 as nc
import numpy as np

def surface_temperatures(filepath):
    print('Calc open surface temperatures of {0}.'.format(filepath))

    # Open result file and read data from it.
    nc_file = nc.Dataset(filepath)
    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size
    time = nc_file.dimensions['time'].size

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

    # Create numpy array and save surface data from netCDF file to it.
    temp = np.zeros((dim1, dim0))
    temp[:,:] = T[(time-1):time,(dim2-1),:,:]

    nc_file.close()

    if os.path.isfile('init.nc') == False:
        print('init.nc does not exist.')
        print('Aborting.')
        exit()

    # Open init file and read data from it.
    nc_file = nc.Dataset('init.nc')
    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size
    time = nc_file.dimensions['time'].size

    surface = nc_file.variables['surface']

    # Create numpy array and save surface data from netCDF file to it.
    skull = np.zeros((dim1, dim0))
    skull[:,:] = surface[(time-1):time,(dim2-1),:,:]

    nc_file.close()

    if np.count_nonzero(skull == 1) != 0:
        temp_mean = np.mean(temp[np.where(skull == 1)])
        temp_max = np.max(temp[np.where(skull == 1)])
        temp_min = np.min(temp[np.where(skull == 1)])
        print('Mean temp of open surface: {0}.'.format(temp_mean))
        print('Max temp of open surface: {0}.'.format(temp_max))
        print('Min temp of open surface: {0}.'.format(temp_min))

        config = configparser.ConfigParser()

        config['Open_Surface'] = {}
        config['Open_Surface']['Mean'] = str(temp_mean)
        config['Open_Surface']['Max'] = str(temp_max)
        config['Open_Surface']['Min'] = str(temp_min)

        filepath = os.path.splitext(filepath)[0]
        filepath += '_results.dat'

        print('Write results to {0}.'.format(filepath))

        with open(filepath, 'w') as configfile:
            config.write(configfile)
    else:
        print('No open skull specified.')

    print('Done.')

def tumor_temperatures(filepath):
    print('Calc tumor temperatures of {0}.'.format(filepath))

    # Open result file and read data from it.
    nc_file = nc.Dataset(filepath)
    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size
    time = nc_file.dimensions['time'].size

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

    # Create numpy array and save surface data from netCDF file to it.
    temp = np.zeros((dim2, dim1, dim0))
    temp[:,:] = T[(time-1):time,:,:,:]

    nc_file.close()

    if os.path.isfile('region.nc') == False:
        print('region.nc does not exist.')
        print('Aborting.')
        exit()

    # Open region file and read data from it.
    nc_file = nc.Dataset('region.nc')
    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size
    time = nc_file.dimensions['time'].size

    region = nc_file.variables['region']

    # Create numpy array and save surface data from netCDF file to it.
    tumor = np.zeros((dim2, dim1, dim0))
    tumor[:,:] = region[(time-1):time,:,:,:]

    nc_file.close()

    if np.count_nonzero(tumor == 1) != 0:
        temp_mean = np.mean(temp[np.where(tumor == 1)])
        temp_max = np.max(temp[np.where(tumor == 1)])
        temp_min = np.min(temp[np.where(tumor == 1)])
        print('Mean temp of tumor: {0}.'.format(temp_mean))
        print('Max temp of tumor: {0}.'.format(temp_max))
        print('Min temp of tumor: {0}.'.format(temp_min))

        config = configparser.ConfigParser()

        config['Tumor'] = {}
        config['Tumor']['Mean'] = str(temp_mean)
        config['Tumor']['Max'] = str(temp_max)
        config['Tumor']['Min'] = str(temp_min)

        filepath = os.path.splitext(filepath)[0]
        filepath += '_results.dat'

        print('Write results to {0}.'.format(filepath))

        with open(filepath, 'a') as configfile:
            config.write(configfile)
    else:
        print('No tumor specified.')

    print('Done.')

def brain_temperatures(filepath):
    print('Calc brain temperatures of {0}.'.format(filepath))

    # Open result file and read data from it.
    nc_file = nc.Dataset(filepath)
    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size
    time = nc_file.dimensions['time'].size

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

    # Create numpy array and save surface data from netCDF file to it.
    temp = np.zeros((dim2, dim1, dim0))
    temp[:,:] = T[(time-1):time,:,:,:]

    nc_file.close()

    if os.path.isfile('region.nc') == False:
        print('region.nc does not exist.')
        print('Aborting.')
        exit()

    # Open region file and read data from it.
    nc_file = nc.Dataset('region.nc')
    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size
    time = nc_file.dimensions['time'].size

    region = nc_file.variables['region']

    # Create numpy array and save surface data from netCDF file to it.
    tumor = np.zeros((dim2, dim1, dim0))
    tumor[:,:] = region[(time-1):time,:,:,:]

    nc_file.close()

    if np.count_nonzero(tumor == 1) != 0:
        temp_mean = np.mean(temp[np.where(tumor == 0)])
        temp_max = np.max(temp[np.where(tumor == 0)])
        temp_min = np.min(temp[np.where(tumor == 0)])
    else:
        temp_mean = np.mean(temp)
        temp_max = np.max(temp)
        temp_min = np.min(temp)

    print('Mean temp of brain: {0}.'.format(temp_mean))
    print('Max temp of brain: {0}.'.format(temp_max))
    print('Min temp of brain: {0}.'.format(temp_min))

    config = configparser.ConfigParser()

    config['Brain'] = {}
    config['Brain']['Mean'] = str(temp_mean)
    config['Brain']['Max'] = str(temp_max)
    config['Brain']['Min'] = str(temp_min)

    filepath = os.path.splitext(filepath)[0]
    filepath += '_results.dat'

    print('Write results to {0}.'.format(filepath))

    with open(filepath, 'a') as configfile:
        config.write(configfile)

    print('Done.')

def domain_temperatures(filepath):
    print('Calc domain temperatures of {0}.'.format(filepath))

    # Open result file and read data from it.
    nc_file = nc.Dataset(filepath)
    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size
    time = nc_file.dimensions['time'].size

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

    # Create numpy array and save surface data from netCDF file to it.
    temp = np.zeros((dim2, dim1, dim0))
    temp[:,:] = T[(time-1):time,:,:,:]

    nc_file.close()

    temp_mean = np.mean(temp)
    temp_max = np.max(temp)
    temp_min = np.min(temp)
    print('Mean temp of domain: {0}.'.format(temp_mean))
    print('Max temp of domain: {0}.'.format(temp_max))
    print('Min temp of domain: {0}.'.format(temp_min))

    config = configparser.ConfigParser()

    config['Domain'] = {}
    config['Domain']['Mean'] = str(temp_mean)
    config['Domain']['Max'] = str(temp_max)
    config['Domain']['Min'] = str(temp_min)

    filepath = os.path.splitext(filepath)[0]
    filepath += '_results.dat'

    print('Write results to {0}.'.format(filepath))

    with open(filepath, 'a') as configfile:
        config.write(configfile)

    print('Done.')


def csv_result_temperatures(filepath, csv):
    csv = os.path.join(csv, csv + '.csv')
    print('Calc temperatures of {0}.'.format(csv))

    # Open results file (thermography).
    temp = np.genfromtxt(csv, delimiter=',')
    temp_mean = np.mean(temp[np.where(temp != 0)])
    temp_max = np.max(temp[np.where(temp != 0)])
    temp_min = np.min(temp[np.where(temp != 0)])
    print('Mean temp: {0}.'.format(temp_mean))
    print('Max temp: {0}.'.format(temp_max))
    print('Min temp: {0}.'.format(temp_min))

    config = configparser.ConfigParser()

    section = str(os.path.basename(csv))

    config[section] = {}
    config[section]['Mean'] = str(temp_mean)
    config[section]['Max'] = str(temp_max)
    config[section]['Min'] = str(temp_min)

    filepath = os.path.splitext(filepath)[0]
    filepath += '_results.dat'

    print('Write results to {0}.'.format(filepath))

    with open(filepath, 'a') as configfile:
        config.write(configfile)

    print('Done.')


def main():
    filepath = ''
    # Check if path to netCDF file (i.e. results) is provided,
    # if file exists and if file has .nc extension.
    if len(sys.argv) > 1:
        if os.path.isfile(sys.argv[1]) == True:
            if os.path.splitext(sys.argv[1])[1] == '.nc':
                filepath = sys.argv[1]
            else:
                print(sys.argv[1], 'does not have .nc extension.')
        else:
            print(sys.argv[1], 'does not exist.')
    else:
        print('No command line argument for netCDF file provided.')

    if filepath == '':
        print('Usage: python3', sys.argv[0], '<PATH/TO/FILE>')
        print('Aborting.')
        exit()

    surface_temperatures(filepath)
    tumor_temperatures(filepath)

if __name__ == '__main__':
    main()
