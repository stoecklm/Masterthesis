import os
import sys

import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import netCDF4 as nc
import numpy as np

def plot_3d_surface(a, params, title, filepath):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    COORD_NODE_FIRST = params['COORD_NODE_FIRST']
    COORD_NODE_LAST = params['COORD_NODE_LAST']
    DIM = params['N_NODES']
    x, y = np.meshgrid(np.linspace(COORD_NODE_FIRST[0], COORD_NODE_LAST[0],
                                   DIM[0]),
                       np.linspace(COORD_NODE_FIRST[1], COORD_NODE_LAST[1],
                                   DIM[1]))
    # Label for axis.
    ax.set_xlabel('x in m')
    ax.set_ylabel('y in m')
    # Title.
    fig.suptitle('Surface Temperature in deg C for\n' + title, fontsize=12)
    # Plot surface.
    surf = ax.plot_surface(x, y, a, cmap=cm.viridis,
                           linewidth=0, antialiased=False)
    # Customize z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, orientation='vertical', shrink=0.5, aspect=20)
    # Save plot to file.
    print('Save figure to {}.'.format(filepath))
    plt.savefig(filepath)

def plot_heatmap():
    pass

def plot_tumor():
    pass

def plot_surface(filepath, params):
    print('Plotting {0}.'.format(filepath))

    # Open netCDF file and read data from it.
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
    temperature = np.zeros((dim1, dim0))
    temperature[:,:] = T[(time-1):time,(dim2-1),:,:]

    nc_file.close()

    filepath = os.path.splitext(filepath)[0]
    title = filepath
    filepath_heatmap = filepath
    filepath += '_py_surface.eps'
    filepath_heatmap += '_py_heatmap.eps'

    plot_3d_surface(temperature, params, title, filepath)

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

    plot_surface(filepath)

if __name__ == '__main__':
    main()
