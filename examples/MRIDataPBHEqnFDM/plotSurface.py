import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import netCDF4 as nc
import numpy as np
import os
import sys

def plot_surface(filepath):
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
    a = np.zeros((dim1, dim0))
    a[:,:] = T[(time-1):time,(dim2-1),:,:]

    nc_file.close()

    filepath = os.path.splitext(filepath)[0]
    plot_title = filepath
    filepath_heatmap = filepath
    filepath += '_py_surface.eps'
    filepath_heatmap += '_py_heatmap.eps'

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x, y = np.meshgrid(np.linspace(0, dim0-1, dim0), \
                       np.linspace(0, dim1-1, dim1))
    # Label for axis.
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    # Title.
    fig.suptitle('Surface Temperature [deg C] for\n' + \
                 plot_title, fontsize=12)
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
    # Plot heatmap.
    plt.gcf().clear()
    plt.title('Surface Temperature [deg C] for\n' + plot_title)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.imshow(a, cmap=cm.viridis, interpolation='nearest')
    print('Save figure to {}.'.format(filepath_heatmap))
    plt.colorbar()
    plt.savefig(filepath_heatmap)

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
