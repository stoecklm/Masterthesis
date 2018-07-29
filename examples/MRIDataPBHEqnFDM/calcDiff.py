import configparser
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
from scipy.interpolate import RegularGridInterpolator

from helperFunctions import temperature_array_from_result

CMAP = cm.viridis

def save_as_netcdf(data, filepath):
    print('Save data to {}.'.format(filepath))
    nc_file = nc.Dataset(filepath, 'w', format='NETCDF3_CLASSIC')
    nc_file.createDimension('time')
    for dim in range(0, len(data.shape)):
        nc_file.createDimension('nNodes_' + str(dim), data.shape[2-dim])
    nNodes = []
    nNodes.append('time')
    for dim in range(len(data.shape), 0, -1):
        nNodes.append('nNodes_' + str(dim-1))
    init_values = nc_file.createVariable('TDiff', 'f8', nNodes)
    init_values[0,] = data

    nc_file.close()

def plot_diff_on_surface(data, filepath):
    print('Save figure to {}.'.format(filepath))
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x, y = np.meshgrid(np.linspace(0, data.shape[1]-1, data.shape[1]),
                       np.linspace(0, data.shape[0]-1, data.shape[0]))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    fig.suptitle('Surface Temperature Difference in deg C', fontsize=12)
    surf = ax.plot_surface(x, y, data, cmap=CMAP,
                           linewidth=0, antialiased=False)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    ticks = np.linspace(data.min(), data.max(), 11)
    fig.colorbar(surf, ticks=ticks, orientation='vertical', shrink=0.75,
                 aspect=20)
    plt.savefig(filepath)
    plt.gcf().clear()
    plt.close()

def get_regular_grid_interpolator(data):
    x = np.linspace(0, data.shape[0]-1, data.shape[0])
    y = np.linspace(0, data.shape[1]-1, data.shape[1])
    z = np.linspace(0, data.shape[2]-1, data.shape[2])

    reg_grid_interpolator = RegularGridInterpolator((x, y, z), data,
                                                    method='linear',
                                                    bounds_error=False,
                                                    fill_value=0)

    return reg_grid_interpolator

def diff_of_two_equal_files(temp_1, temp_2, filepath):
    temp_diff = np.subtract(temp_1, temp_2)
    temp_diff = np.absolute(temp_diff)

    config = configparser.ConfigParser()
    config.optionxform = str

    mean_diff = np.mean(temp_diff)
    max_diff = np.amax(temp_diff)
    min_diff = np.amin(temp_diff)
    print('Volume:  mean(abs(diff)) = {:02.3e}.'.format(mean_diff))
    print('Volume:   max(abs(diff)) = {:02.3e}.'.format(max_diff))
    print('Volume:   min(abs(diff)) = {:02.3e}.'.format(min_diff))
    config['Volume'] = {}
    config['Volume']['Mean'] = str(mean_diff)
    config['Volume']['Max'] = str(max_diff)
    config['Volume']['Min'] = str(min_diff)

    mean_diff = np.mean(temp_diff[-1,:,:])
    max_diff = np.amax(temp_diff[-1,:,:])
    min_diff = np.amin(temp_diff[-1,:,:])
    print('Surface: mean(abs(diff)) = {:02.3e}.'.format(mean_diff))
    print('Surface:  max(abs(diff)) = {:02.3e}.'.format(max_diff))
    print('Surface:  min(abs(diff)) = {:02.3e}.'.format(min_diff))
    config['Surface'] = {}
    config['Surface']['Mean'] = str(mean_diff)
    config['Surface']['Max'] = str(max_diff)
    config['Surface']['Min'] = str(min_diff)

    print('Write results to {0}.'.format(filepath))

    with open(filepath, 'w') as configfile:
        config.write(configfile)

    return temp_diff

def diff_of_file_dimensions_not_equal(temp_1, temp_2, filepath):
    if temp_1.shape[0] >= temp_2.shape[0] \
        and temp_1.shape[1] >= temp_2.shape[1] \
        and temp_1.shape[2] >= temp_2.shape[2]:
        temp_dense = temp_1
        temp_sparse = temp_2
    elif temp_1.shape[0] <= temp_2.shape[0] \
        and temp_1.shape[1] <= temp_2.shape[1] \
        and temp_1.shape[2] <= temp_2.shape[2]:
        temp_sparse = temp_1
        temp_dense = temp_2
    else:
        temp_diff = diff_of_weird_file_dimensions(temp_1, temp_2, filepath)
        return temp_diff

    dense_interpolator = get_regular_grid_interpolator(temp_dense)

    dim0_sparse, dim1_sparse, dim2_sparse = reversed(temp_sparse.shape)
    dim0_dense, dim1_dense, dim2_dense = reversed(temp_dense.shape)

    dim0_lin = np.linspace(0, dim0_dense-1, dim0_sparse)
    dim1_lin = np.linspace(0, dim1_dense-1, dim1_sparse)
    dim2_lin = np.linspace(0, dim2_dense-1, dim2_sparse)

    dim2_mesh, dim1_mesh, dim0_mesh = np.meshgrid(dim2_lin, dim1_lin, dim0_lin,
                                                  indexing='ij')

    temp_dense_to_sparse = dense_interpolator((dim2_mesh, dim1_mesh, dim0_mesh))

    temp_diff = diff_of_two_equal_files(temp_sparse, temp_dense_to_sparse,
                                        filepath)

    return temp_diff

def diff_of_weird_file_dimensions(temp_1, temp_2, filepath):
    dim0_1, dim1_1, dim2_1 = reversed(temp_1.shape)
    dim0_2, dim1_2, dim2_2 = reversed(temp_2.shape)

    if dim0_1 > dim0_2:
        dim0_sparse = dim0_2
    else:
        dim0_sparse = dim0_1

    if dim1_1 > dim1_2:
        dim1_sparse = dim1_2
    else:
        dim1_sparse = dim1_1

    if dim2_1 > dim2_2:
        dim2_dense = dim2_1
    else:
        dim2_sparse = dim2_1

    temp_1_interpolator = get_regular_grid_interpolator(temp_1)
    temp_2_interpolator = get_regular_grid_interpolator(temp_2)

    dim0_lin = np.linspace(0, dim0_1-1, dim0_sparse)
    dim1_lin = np.linspace(0, dim1_1-1, dim1_sparse)
    dim2_lin = np.linspace(0, dim2_1-1, dim2_sparse)

    dim2_mesh, dim1_mesh, dim0_mesh = np.meshgrid(dim2_lin, dim1_lin, dim0_lin,
                                                  indexing='ij')

    temp_1_new_dim = temp_1_interpolator((dim2_mesh, dim1_mesh, dim0_mesh))

    dim0_lin = np.linspace(0, dim0_2-1, dim0_sparse)
    dim1_lin = np.linspace(0, dim1_2-1, dim1_sparse)
    dim2_lin = np.linspace(0, dim2_2-1, dim2_sparse)

    dim2_mesh, dim1_mesh, dim0_mesh = np.meshgrid(dim2_lin, dim1_lin, dim0_lin,
                                                  indexing='ij')

    temp_2_new_dim = temp_2_interpolator((dim2_mesh, dim1_mesh, dim0_mesh))

    temp_diff = diff_of_two_equal_files(temp_1_new_dim, temp_2_new_dim,
                                        filepath)

    return temp_diff

def calc_diff(filepaths):
    temp_1 = temperature_array_from_result(filepaths[0])
    temp_2 = temperature_array_from_result(filepaths[1])

    filepath = os.path.splitext(filepaths[0])[0] + '_' \
               + os.path.basename(os.path.splitext(filepaths[1])[0]) \
               + '_diff'
    filepath_dat = filepath + '.dat'

    if temp_1.shape == temp_2.shape:
        temp_diff = diff_of_two_equal_files(temp_1, temp_2, filepath_dat)
    else:
        temp_diff = diff_of_file_dimensions_not_equal(temp_1, temp_2,
                                                      filepath_dat)

    filepath_fig = filepath + '_surface.eps'
    plot_diff_on_surface(temp_diff[-1,:,:], filepath_fig)
    filepath_nc = filepath + '.nc'
    save_as_netcdf(temp_diff, filepath_nc)


def main():
    filepaths = []
    # Check if paths to two netCDF files are provided.
    if len(sys.argv) > 2:
        for elem in range(1, 3):
            if os.path.isfile(sys.argv[elem]) == True:
                if os.path.splitext(sys.argv[elem])[1] == '.nc':
                    filepaths.append(sys.argv[elem])
                else:
                    print(sys.argv[elem],
                          'does not have .nc extension.')
            else:
                print(sys.argv[elem], 'does not exist.')
    else:
        print('Not enough command line arguments for netCDF files.')

    if len(filepaths) < 2:
        print('Usage: python3', sys.argv[0],
              '<PATH/TO/FILE_1> <PATH/TO/FILE_2>')
        print('Aborting.')
        exit()

    calc_diff(filepaths)

    print('Done.')

if __name__ == '__main__':
    main()
