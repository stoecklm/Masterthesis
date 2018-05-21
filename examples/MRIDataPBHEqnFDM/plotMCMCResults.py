import os
import sys

import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.ticker import FormatStrFormatter
import netCDF4 as nc
import numpy as np
import numpy.ma as ma


def read_variable_from_netcdf_file(name_of_variable, filepath):
    print('Read {}.'.format(filepath))
    nc_file = nc.Dataset(filepath)
    iterations = nc_file.dimensions['iterations'].size

    try:
        nc_variable = nc_file.variables[name_of_variable]
    except KeyError:
        print('* ERROR: Variable not found in this file.')
        print('Aborting.')
        exit()

    data_from_file = np.zeros(iterations)
    data_from_file = nc_variable[:]

    nc_file.close()

    return data_from_file

def read_2d_tested_variables(filepath):
    print('Read {}.'.format(filepath))
    nc_file = nc.Dataset(filepath)
    iterations = nc_file.dimensions['iterations'].size

    possible_names = ['omega', 'rho_c', 'q_brain']
    found_name = False
    for name in possible_names:
        try:
            nc_variable_1 = nc_file.variables[name + '_brain']
            nc_variable_2 = nc_file.variables[name + '_tumor']
            found_var = True
            break
        except KeyError:
            pass

    if found_var == False:
        print('Aborting.')
        exit()

    data_from_file_1 = np.zeros(iterations)
    data_from_file_2 = np.zeros(iterations)
    data_from_file_1 = nc_variable_1[:]
    data_from_file_2 = nc_variable_2[:]

    nc_file.close()

    return data_from_file_1, name + '_brain', data_from_file_2, name + '_tumor'

def read_tested_variables_from_netcdf_file(filepath):
    print('Read {}.'.format(filepath))
    nc_file = nc.Dataset(filepath)
    iterations = nc_file.dimensions['iterations'].size

    possible_names = ['lambda_bt', 'T_blood', 'h']
    found_name = False
    for name in possible_names:
        try:
            nc_variable = nc_file.variables[name]
            found_var = True
            found_name = name
            break
        except KeyError:
            pass

    if found_var == True:
        data_from_file_1 = np.zeros(iterations)
        data_from_file_1 = nc_variable[:]
        data_from_file_2 = np.zeros(0)
        name_1 = name
        name_2 = ''
        nc_file.close()
    else:
        nc_file.close()
        data_from_file_1, name_1, data_from_file_2, name_2 =  read_2d_tested_variables(filepath)

    return data_from_file_1, name_1, data_from_file_2, name_2

def main():
    if len(sys.argv) > 1:
        filepath = sys.argv[1]
        if os.path.isfile(filepath) != True:
            print(filepath, 'does not exist.')
            print('Aborting.')
            exit()
        else:
            filepath = sys.argv[1]
    else:
        print('No command line argument for folder provided.')
        print('Aborting.')
        exit()

    results_name = os.path.basename(filepath).split('.')[0]

    l2_norm = read_variable_from_netcdf_file('L2-norm', filepath)
    T_normal = read_variable_from_netcdf_file('T_normal', filepath)
    T_tumor = read_variable_from_netcdf_file('T_tumor', filepath)
    T_vessel = read_variable_from_netcdf_file('T_vessel', filepath)

    var_1, name_1, var_2, name_2 = read_tested_variables_from_netcdf_file(filepath)
    if len(var_2) == 0:
        fig, ax = plt.subplots()
        ax.plot(var_1, l2_norm, 'o')
        ax.set_xlabel(name_1)
        ax.set_ylabel('L2-norm')
        plt.savefig(results_name + '_l2_norm.eps')
        plt.close()
        fig, ax = plt.subplots()
        ax.plot(var_1, T_tumor, marker='s', color='orange', label='T_tumor', linestyle='None')
        ax.plot(var_1, T_normal, marker='o', color='darkorchid', label='T_normal', linestyle='None')
        ax.plot(var_1, T_vessel, marker='^', color='lightseagreen', label='T_vessel', linestyle='None')
        ax.set_xlabel(name_1)
        ax.set_ylabel('Temperature in deg C')
        plt.legend()
        plt.savefig(results_name + '_T.eps')
        plt.close()
    else:
        pass

if __name__ == '__main__':
    main()
