import configparser
import os
import sys

import netCDF4 as nc
import numpy as np

from helperFunctions import temperature_array_from_result
from helperFunctions import surface_temperature_array_from_result

DEPTH = 10

def print_results(section, temp_mean, temp_max, temp_min, temp_std_dev):
    print('Mean temp of {}: {:02.3f}.'.format(section, temp_mean))
    print('Max temp of {}: {:02.3f}.'.format(section, temp_max))
    print('Min temp of {}: {:02.3f}.'.format(section, temp_min))
    print('Std dev temp of {}: {:02.4f}.'.format(section, temp_std_dev))

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

def array_from_file(filepath, name):
    filepath += '.nc'
    if os.path.isfile(filepath) == False:
        print(filepath, 'does not exist.')
        print('Aborting.')
        exit()

    nc_file = nc.Dataset(filepath)
    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size
    nc_var = nc_file.variables[name]

    np_array = np.zeros((dim2, dim1, dim0))
    np_array[:,:,:] = nc_var[-1,:,:,:]

    nc_file.close()

    return np_array

def region_array_from_file(filepath):
    region = array_from_file(filepath, 'region')

    return region

def surface_array_from_file(filepath):
    surface = array_from_file(filepath, 'surface')

    return surface

def surface_vessels_array_from_file(filepath):
    vessels = array_from_file(filepath, 'vessels')
    vessels = vessels[-1,:,:]

    return vessels

def calc_open_surface_temperatures(temp, surface):
    skull = surface[-1,:,:]
    if np.count_nonzero(skull == 1) != 0:
        temp_mean = np.mean(temp[np.where(skull == 1)])
        temp_max = np.max(temp[np.where(skull == 1)])
        temp_min = np.min(temp[np.where(skull == 1)])
        temp_std_dev = np.std(temp[np.where(skull == 1)])
    else:
        print('No open surface specified.')
        temp_mean = -1.0
        temp_max = -1.0
        temp_min = -1.0
        temp_std_dev = -1.0

    return temp_mean, temp_max, temp_min, temp_std_dev

def open_surface_temperatures(filepath, filepath_init, do_print=True, do_write=True):
    print()
    print('Calc open surface temperatures of {0}.'.format(filepath))

    temp = surface_temperature_array_from_result(filepath)
    surface = surface_array_from_file(filepath_init)
    temp_mean, temp_max, temp_min, temp_std_dev = calc_open_surface_temperatures(temp, surface)

    if do_print == True:
        print_results('open surface', temp_mean, temp_max, temp_min,
                      temp_std_dev)
    if do_write == True:
        filepath = os.path.splitext(filepath)[0] + '_results.dat'
        write_results_to_file('Open_Surface', temp_mean, temp_max, temp_min,
                              temp_std_dev, filepath, 'w')

    print('Done.')

    return temp_mean

def calc_tumor_temperatures(temp, tumor):
    if np.count_nonzero(tumor == 1) != 0:
        temp_mean = np.mean(temp[np.where(tumor == 1)])
        temp_max = np.max(temp[np.where(tumor == 1)])
        temp_min = np.min(temp[np.where(tumor == 1)])
        temp_std_dev = np.std(temp[np.where(tumor == 1)])
    else:
        print('No tumor specified.')
        temp_mean = -1.0
        temp_max = -1.0
        temp_min = -1.0
        temp_std_dev = -1.0

    return temp_mean, temp_max, temp_min, temp_std_dev

def tumor_temperatures(filepath, region_filepath, do_print=True, do_write=True):
    print()
    print('Calc tumor temperatures of {0}.'.format(filepath))

    temp = temperature_array_from_result(filepath)
    tumor = region_array_from_file(region_filepath)
    temp_mean, temp_max, temp_min, temp_std_dev = calc_tumor_temperatures(temp, tumor)

    if do_print == True:
        print_results('tumor', temp_mean, temp_max, temp_min, temp_std_dev)
    if do_write == True:
        filepath = os.path.splitext(filepath)[0] + '_results.dat'
        write_results_to_file('Tumor', temp_mean, temp_max, temp_min,
                              temp_std_dev, filepath, 'a')

    print('Done.')

    return temp_mean

def calc_tumor_near_surface_temperatures(temp, tumor):
    dim2 = np.any(tumor, axis=(1, 2))
    try:
        min2, max2 = np.where(dim2)[0][[0, -1]]
    except IndexError:
        return -1.0, -1.0, -1.0, -1.0

    depth = DEPTH

    tumor = tumor[max2-depth+1:max2+1,:,:]
    temp = temp[max2-depth+1:max2+1,:,:]

    if np.count_nonzero(tumor == 1) != 0:
        temp_mean = np.mean(temp[np.where(tumor == 1)])
        temp_max = np.max(temp[np.where(tumor == 1)])
        temp_min = np.min(temp[np.where(tumor == 1)])
        temp_std_dev = np.std(temp[np.where(tumor == 1)])
    else:
        print('No tumor specified.')
        temp_mean = -1.0
        temp_max = -1.0
        temp_min = -1.0
        temp_std_dev = -1.0

    return temp_mean, temp_max, temp_min, temp_std_dev

def tumor_near_surface_temperatures(filepath, region_filepath, do_print=True, do_write=True):
    print()
    print('Calc tumor temperatures near surface of {0}.'.format(filepath))

    temp = temperature_array_from_result(filepath)
    tumor = region_array_from_file(region_filepath)
    temp_mean, temp_max, temp_min, temp_std_dev = calc_tumor_near_surface_temperatures(temp, tumor)

    depth = DEPTH

    if do_print == True:
        section = 'tumor near surface (first ' + str(depth) + ' nodes)'
        print_results(section, temp_mean, temp_max, temp_min, temp_std_dev)
    if do_write == True:
        filepath = os.path.splitext(filepath)[0] + '_results.dat'
        section = 'Tumor_Near_Surface_' + str(depth) + '_Depth'
        write_results_to_file(section, temp_mean, temp_max, temp_min,
                              temp_std_dev, filepath, 'a')
    print('Done.')

    return temp_mean

def calc_brain_temperatures(temp, tumor):
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

    return temp_mean, temp_max, temp_min, temp_std_dev

def brain_temperatures(filepath, region_filepath, do_print=True, do_write=True):
    print()
    print('Calc brain temperatures of {0}.'.format(filepath))

    temp = temperature_array_from_result(filepath)
    tumor = region_array_from_file(region_filepath)
    temp_mean, temp_max, temp_min, temp_std_dev = calc_brain_temperatures(temp, tumor)

    if do_print == True:
        print_results('brain', temp_mean, temp_max, temp_min, temp_std_dev)
    if do_write == True:
        filepath = os.path.splitext(filepath)[0] + '_results.dat'
        write_results_to_file('Brain', temp_mean, temp_max, temp_min, temp_std_dev,
                              filepath, 'a')

    print('Done.')

    return temp_mean

def calc_domain_temperatures(temp):
    temp_mean = np.mean(temp)
    temp_max = np.max(temp)
    temp_min = np.min(temp)
    temp_std_dev = np.std(temp)

    return temp_mean, temp_max, temp_min, temp_std_dev

def domain_temperatures(filepath, do_print=True, do_write=True):
    print()
    print('Calc domain temperatures of {0}.'.format(filepath))

    temp = temperature_array_from_result(filepath)
    temp_mean, temp_max, temp_min, temp_std_dev = calc_domain_temperatures(temp)

    if do_print == True:
        print_results('domain', temp_mean, temp_max, temp_min, temp_std_dev)
    if do_write == True:
        filepath = os.path.splitext(filepath)[0] + '_results.dat'
        write_results_to_file('Domain', temp_mean, temp_max, temp_min, temp_std_dev,
                              filepath, 'a')

    print('Done.')

    return temp_mean

def calc_csv_result_temperatures(temp):
    temp_mean = np.mean(temp[np.where(temp != 0)])
    temp_max = np.max(temp[np.where(temp != 0)])
    temp_min = np.min(temp[np.where(temp != 0)])
    temp_std_dev = np.std(temp[np.where(temp != 0)])

    return temp_mean, temp_max, temp_min, temp_std_dev

def csv_result_temperatures(filepath, csv, do_print=True, do_write=True):
    csv = os.path.join(csv, 'thermo.csv')
    print()
    print('Calc temperatures of {0}.'.format(csv))

    # Open results file (thermography).
    temp = np.genfromtxt(csv, delimiter=',')
    temp = np.nan_to_num(temp)

    temp_mean, temp_max, temp_min, temp_std_dev = calc_csv_result_temperatures(temp)

    if do_print == True:
        print_results('thermo.csv', temp_mean, temp_max, temp_min, temp_std_dev)
    if do_write == True:
        section = str(os.path.basename(csv))
        filepath = os.path.splitext(filepath)[0] + '_results.dat'
        write_results_to_file(section, temp_mean, temp_max, temp_min, temp_std_dev,
                              filepath, 'a')

    print('Done.')

    return temp_mean

def calc_vessels_temperatures(temp, vessels):
    temp_mean = np.mean(temp[np.where(vessels == 1)])
    temp_max = np.max(temp[np.where(vessels == 1)])
    temp_min = np.min(temp[np.where(vessels == 1)])
    temp_std_dev = np.std(temp[np.where(vessels == 1)])

    return temp_mean, temp_max, temp_min, temp_std_dev

def vessels_temperatures(filepath_nc, filepath_vessels, do_print=True, do_write=True):
    print()
    print('Calc vessel temperatures of {0}.'.format(filepath_nc))

    temp = surface_temperature_array_from_result(filepath_nc)
    vessels = surface_vessels_array_from_file(filepath_vessels)
    temp_mean, temp_max, temp_min, temp_std_dev = calc_vessels_temperatures(temp, vessels)

    if do_print == True:
        print_results('vessels', temp_mean, temp_max, temp_min, temp_std_dev)
    if do_write == True:
        filepath = os.path.splitext(filepath_nc)[0] + '_results.dat'
        write_results_to_file('Vessel', temp_mean, temp_max, temp_min, temp_std_dev,
                              filepath, 'a')

    print('Done.')

    return temp_mean

def calc_non_vessels_temperatures(temp, vessels):
    temp_mean = np.mean(temp[np.where(vessels == 0)])
    temp_max = np.max(temp[np.where(vessels == 0)])
    temp_min = np.min(temp[np.where(vessels == 0)])
    temp_std_dev = np.std(temp[np.where(vessels == 0)])

    return temp_mean, temp_max, temp_min, temp_std_dev

def non_vessels_temperatures(filepath_nc, filepath_vessels, do_print=True, do_write=True):
    print()
    print('Calc non-vessel temperatures of {0}.'.format(filepath_nc))

    temp = surface_temperature_array_from_result(filepath_nc)
    vessels = surface_vessels_array_from_file(filepath_vessels)
    temp_mean, temp_max, temp_min, temp_std_dev = calc_non_vessels_temperatures(temp, vessels)

    if do_print == True:
        print_results('non-vessels', temp_mean, temp_max, temp_min, temp_std_dev)
    if do_write == True:
        filepath = os.path.splitext(filepath_nc)[0]
        filepath += '_results.dat'
        write_results_to_file('Non_Vessel', temp_mean, temp_max, temp_min,
                              temp_std_dev, filepath, 'a')

    print('Done.')

    return temp_mean

def calc_l2_norm(filepath_nc, T_normal, T_tumor, T_vessel,
                 T_normal_thermo, T_tumor_thermo, T_vessel_thermo):
    print()
    print('Calc L2-norm of {0}.'.format(filepath_nc))
    if T_normal_thermo == -1.0 or \
       T_tumor_thermo == -1.0 or \
       T_vessel_thermo == -1.0:
        print('No target values specified.')
    elif T_vessel == -1.0:
        print('No vessels specified.')
    elif T_normal == -1.0:
        print('No normal region specified.')
    elif T_tumor == -1.0:
        print('No tumor specified.')
    else:
        scafes_values = np.asarray([T_normal, T_tumor, T_vessel])
        target_values = np.asarray([T_normal_thermo, T_tumor_thermo, T_vessel_thermo])
        l2_norm = np.linalg.norm(np.subtract(scafes_values, target_values), 2)
        print('Target values:', target_values)
        print('T_normal: {:02.3f}.'.format(T_normal))
        print('T_tumor: {:02.3f}.'.format(T_tumor))
        print('T_vessel: {:02.3f}.'.format(T_vessel))
        print('L2-norm: {:02.3f}.'.format(l2_norm))
        filepath = os.path.splitext(filepath_nc)[0] + '_results.dat'
        config = configparser.ConfigParser()
        config['L2-norm'] = {}
        config['L2-norm']['Target values'] = str(target_values)
        config['L2-norm']['T_normal'] = str(T_normal)
        config['L2-norm']['T_tumor'] = str(T_tumor)
        config['L2-norm']['T_vessel'] = str(T_vessel)
        config['L2-norm']['L2-norm'] = str(l2_norm)

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
