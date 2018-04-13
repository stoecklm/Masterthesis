import configparser
import os
import sys

import netCDF4 as nc
import numpy as np

from helperFunctions import temperature_array_from_result
from helperFunctions import surface_temperature_array_from_result

def print_results(section, temp_mean, temp_max, temp_min, temp_std_dev):
    print('Mean temp of {0}: {1}.'.format(section, temp_mean))
    print('Max temp of {0}: {1}.'.format(section, temp_max))
    print('Min temp of {0}: {1}.'.format(section, temp_min))
    print('Std dev temp of {0}: {1}.'.format(section, temp_std_dev))

def write_results_to_file(section, temp_mean, temp_max, temp_min, temp_std_dev,
                          filepath, file_mode):
    config = configparser.ConfigParser()
    config[section] = {}
    config[section]['Mean'] = str(temp_mean)
    config[section]['Max'] = str(temp_max)
    config[section]['Min'] = str(temp_min)
    config[section]['Std_Dev'] = str(temp_std_dev)

    print('Write results to {0}.'.format(filepath))

    with open(filepath, file_mode) as configfile:
        config.write(configfile)

def region_array_from_file(filepath):
    filepath += '.nc'
    if os.path.isfile(filepath) == False:
        print(filepath, 'does not exist.')
        print('Aborting.')
        exit()

    nc_file = nc.Dataset(filepath)
    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size
    region = nc_file.variables['region']

    tumor = np.zeros((dim2, dim1, dim0))
    tumor[:,:,:] = region[-1,:,:,:]

    nc_file.close()

    return tumor

def surface_temperatures(filepath):
    print('Calc open surface temperatures of {0}.'.format(filepath))

    if os.path.isfile('init.nc') == False:
        print('init.nc does not exist.')
        print('Aborting.')
        exit()

    temp = surface_temperature_array_from_result(filepath)

    nc_file = nc.Dataset('init.nc')
    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size
    surface = nc_file.variables['surface']

    skull = np.zeros((dim1, dim0))
    skull[:,:] = surface[-1,-1,:,:]

    nc_file.close()

    if np.count_nonzero(skull == 1) != 0:
        temp_mean = np.mean(temp[np.where(skull == 1)])
        temp_max = np.max(temp[np.where(skull == 1)])
        temp_min = np.min(temp[np.where(skull == 1)])
        temp_std_dev = np.std(temp[np.where(skull == 1)])

        filepath = os.path.splitext(filepath)[0]
        filepath += '_results.dat'

        print_results('open surface', temp_mean, temp_max, temp_min,
                      temp_std_dev)
        write_results_to_file('Open_Surface', temp_mean, temp_max, temp_min,
                              temp_std_dev, filepath, 'w')

    else:
        print('No open skull specified.')
        temp_mean = -1.0

    print('Done.')

    return temp_mean

def tumor_temperatures(filepath, region_filepath):
    print('Calc tumor temperatures of {0}.'.format(filepath))

    temp = temperature_array_from_result(filepath)
    tumor = region_array_from_file(region_filepath)

    if np.count_nonzero(tumor == 1) != 0:
        temp_mean = np.mean(temp[np.where(tumor == 1)])
        temp_max = np.max(temp[np.where(tumor == 1)])
        temp_min = np.min(temp[np.where(tumor == 1)])
        temp_std_dev = np.std(temp[np.where(tumor == 1)])

        filepath = os.path.splitext(filepath)[0]
        filepath += '_results.dat'

        print_results('tumor', temp_mean, temp_max, temp_min, temp_std_dev)
        write_results_to_file('Tumor', temp_mean, temp_max, temp_min,
                              temp_std_dev, filepath, 'a')
    else:
        print('No tumor specified.')
        temp_mean = -1.0

    print('Done.')

    return temp_mean

def brain_temperatures(filepath, region_filepath):
    print('Calc brain temperatures of {0}.'.format(filepath))

    temp = temperature_array_from_result(filepath)
    tumor = region_array_from_file(region_filepath)

    if np.count_nonzero(tumor == 1) != 0:
        temp_mean = np.mean(temp[np.where(tumor == 0)])
        temp_max = np.max(temp[np.where(tumor == 0)])
        temp_min = np.min(temp[np.where(tumor == 0)])
        temp_std_dev = np.std(temp[np.where(tumor == 0)])
    else:
        temp_mean = np.mean(temp)
        temp_max = np.max(temp)
        temp_min = np.min(temp)
        temp_std_dev = np.std(temp)

    filepath = os.path.splitext(filepath)[0]
    filepath += '_results.dat'

    print_results('brain', temp_mean, temp_max, temp_min, temp_std_dev)
    write_results_to_file('Brain', temp_mean, temp_max, temp_min, temp_std_dev,
                          filepath, 'a')

    print('Done.')

    return temp_mean

def domain_temperatures(filepath):
    print('Calc domain temperatures of {0}.'.format(filepath))

    temp = temperature_array_from_result(filepath)

    temp_mean = np.mean(temp)
    temp_max = np.max(temp)
    temp_min = np.min(temp)
    temp_std_dev = np.std(temp)

    filepath = os.path.splitext(filepath)[0]
    filepath += '_results.dat'

    print_results('domain', temp_mean, temp_max, temp_min, temp_std_dev)
    write_results_to_file('Domain', temp_mean, temp_max, temp_min, temp_std_dev,
                          filepath, 'a')

    print('Done.')

    return temp_mean

def csv_result_temperatures(filepath, csv):
    csv = os.path.join(csv, 'thermo.csv')
    print('Calc temperatures of {0}.'.format(csv))

    # Open results file (thermography).
    temp = np.genfromtxt(csv, delimiter=',')
    temp = np.nan_to_num(temp)
    temp_mean = np.mean(temp[np.where(temp != 0)])
    temp_max = np.max(temp[np.where(temp != 0)])
    temp_min = np.min(temp[np.where(temp != 0)])
    temp_std_dev = np.std(temp[np.where(temp != 0)])

    section = str(os.path.basename(csv))
    filepath = os.path.splitext(filepath)[0]
    filepath += '_results.dat'

    print_results('thermo.csv', temp_mean, temp_max, temp_min, temp_std_dev)
    write_results_to_file(section, temp_mean, temp_max, temp_min, temp_std_dev,
                          filepath, 'a')

    print('Done.')

    return temp_mean

def vessels_temperatures(filepath_nc, vessels):
    print('Calc vessel temperatures of {0}.'.format(filepath_nc))

    temp = surface_temperature_array_from_result(filepath_nc)

    temp_mean = np.mean(temp[np.where(vessels == 1)])
    temp_max = np.max(temp[np.where(vessels == 1)])
    temp_min = np.min(temp[np.where(vessels == 1)])
    temp_std_dev = np.std(temp[np.where(vessels == 1)])
    temp_mean_vessel = temp_mean

    filepath = os.path.splitext(filepath_nc)[0]
    filepath += '_results.dat'

    print_results('vessels', temp_mean, temp_max, temp_min, temp_std_dev)
    write_results_to_file('Vessel', temp_mean, temp_max, temp_min, temp_std_dev,
                          filepath, 'a')

    print('Done.')

    print('Calc non-vessel temperatures of {0}.'.format(filepath_nc))

    temp_mean = np.mean(temp[np.where(vessels == 0)])
    temp_max = np.max(temp[np.where(vessels == 0)])
    temp_min = np.min(temp[np.where(vessels == 0)])
    temp_std_dev = np.std(temp[np.where(vessels == 0)])
    temp_mean_non_vessel = temp_mean

    print_results('non-vessels', temp_mean, temp_max, temp_min, temp_std_dev)
    write_results_to_file('Non_Vessel', temp_mean, temp_max, temp_min,
                          temp_std_dev, filepath, 'a')

    print('Done.')

    return temp_mean_vessel, temp_mean_non_vessel


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
