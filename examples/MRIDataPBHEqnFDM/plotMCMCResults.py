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

    possible_names = ['omega', 'rho_c', 'q']
    found_var = False
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
    found_var = False
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

def plot_2d(var_1, name_1, var_2, name_2, data, results_name, name, color):
    print('Save figure to {}.'.format(results_name + '_' + name + '.eps'))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = "3d")
    xpos = var_1
    ypos = var_2
    zpos = np.zeros(var_1.shape) + (data.min() * 0.999)
    dx = (var_1.max() - var_1.min())/20
    dy = (var_2.max() - var_2.min())/20
    dx = np.ones(var_1.shape)*dx
    dy = np.ones(var_1.shape)*dy
    dz = data - (data.min() * 0.999)
    ax.set_xlabel(name_1)
    ax.set_ylabel(name_2)
    ax.set_zlabel(name)
    ax.set_zlim3d(data.min()*0.999, data.max())
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=color, shade=True, zsort='max')
    plt.savefig(results_name + '_' + name + '.eps')

def plot_2d_l2_norm(var_1, name_1, var_2, name_2, l2_norm, results_name):
    plot_2d(var_1, name_1, var_2, name_2, l2_norm, results_name, 'l2_norm', '#1f77b4')

def plot_2d_temperatures(var_1, name_1, var_2, name_2, temp, results_name, name, color):
    plot_2d(var_1, name_1, var_2, name_2, temp, results_name, name, color)

def plot_all_temperatures(temp_1, temp_2, temp_3, var_1, var_2):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = "3d")
    xpos_1 = var_1 - 1000/3
    xpos_2 = var_1 + 1000/3
    xpos_3 = var_1
    ypos = var_2
    dx = np.ones(var_1.shape)*300/3
    dy = np.ones(var_1.shape)*1000
    zpos = np.zeros(var_1.shape)
    dz_1 = temp_1
    dz_2 = temp_2
    dz_3 = temp_3
    ax.bar3d(xpos_1, ypos, zpos, dx, dy, dz_1, color='blue', shade=True)
    ax.bar3d(xpos_2, ypos, zpos, dx, dy, dz_2, color='green', shade=True)
    ax.bar3d(xpos_3, ypos, zpos, dx, dy, dz_3, color='red', shade=True)
    plt.savefig('test.eps')


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
        plot_2d_l2_norm(var_1, name_1, var_2, name_2, l2_norm, results_name)
        plot_2d_temperatures(var_1, name_1, var_2, name_2, T_tumor, results_name, 'T_tumor', 'orange')
        plot_2d_temperatures(var_1, name_1, var_2, name_2, T_normal, results_name, 'T_normal', 'darkorchid')
        plot_2d_temperatures(var_1, name_1, var_2, name_2, T_vessel, results_name, 'T_vessel', 'lightseagreen')
        #plot_all_temperatures(T_tumor, T_normal, T_vessel, var_1, var_2)


if __name__ == '__main__':
    main()
