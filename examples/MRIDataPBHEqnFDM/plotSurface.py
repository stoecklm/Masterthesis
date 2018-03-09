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
    ticks = np.linspace(a.min(), a.max(), 11)
    fig.colorbar(surf, ticks=ticks, orientation='vertical', shrink=0.75,
                 aspect=20)
    # Invert z axis since the relevant part is colder than the other part
    # and therefore hard to see.
    ax.invert_zaxis()
    # Save plot to file.
    print('Save figure to {}.'.format(filepath))
    plt.savefig(filepath)
    plt.gcf().clear()
    plt.close()

def plot_heatmap(a, params, title, filepath):
    COORD_NODE_FIRST = params['COORD_NODE_FIRST']
    COORD_NODE_LAST = params['COORD_NODE_LAST']
    DIM = params['N_NODES']
    TUMOR_CENTER = params['TUMOR_CENTER']
    RADIUS = params['PARAMETERS']['diameter']/2
    x, y = np.meshgrid(np.linspace(COORD_NODE_FIRST[0], COORD_NODE_LAST[0],
                                   DIM[0]),
                       np.linspace(COORD_NODE_FIRST[1], COORD_NODE_LAST[1],
                                   DIM[1]))
    # Plot heatmap with circle around hole.
    fig, ax = plt.subplots()
    heatmap = ax.pcolormesh(x, y, a, cmap=cm.viridis, rasterized=True)
    try:
        pts = params['HOLE']
        ax.plot(pts[:,0], pts[:,1], color='r', linestyle='dashed')
    except KeyError:
        circle = plt.Circle((TUMOR_CENTER[0], TUMOR_CENTER[1]), RADIUS,
                            color='r', fill=False, linestyle='dashed')
        ax.add_artist(circle)
    # Title.
    fig.suptitle('Heatmap in deg C for\n' + title, fontsize=12)
    # Label for axis.
    ax.set_xlabel('x in m')
    ax.set_ylabel('y in m')
    # Add a color bar which maps values to colors.
    ticks = np.linspace(a.min(), a.max(), 11)
    fig.colorbar(heatmap, ticks=ticks, orientation='vertical', shrink=0.75,
                 aspect=20)
    # Equal gridsize.
    plt.gca().set_aspect('equal', adjustable='box')
    # Save plot to file.
    print('Save figure to {}.'.format(filepath))
    plt.savefig(filepath)
    plt.gcf().clear()
    plt.close()

def plot_tumor(a, params, title, filepath):
    COORD_NODE_FIRST = params['COORD_NODE_FIRST']
    COORD_NODE_LAST = params['COORD_NODE_LAST']
    DIM = params['N_NODES']
    TUMOR_CENTER = params['TUMOR_CENTER']
    RADIUS = params['PARAMETERS']['diameter']/2
    x, z = np.meshgrid(np.linspace(COORD_NODE_FIRST[0], COORD_NODE_LAST[0],
                                   DIM[0]),
                       np.linspace(COORD_NODE_FIRST[2], COORD_NODE_LAST[2],
                                   DIM[2]))
    circle = plt.Circle((TUMOR_CENTER[0], TUMOR_CENTER[2]), RADIUS, color='r',
                         fill=False, linestyle='dashed')
    # Plot heatmap with circle around tumor.
    fig, ax = plt.subplots()
    heatmap = ax.pcolormesh(x, z, a, cmap=cm.viridis, rasterized=True)
    ax.add_artist(circle)
    # Title.
    fig.suptitle('Heatmap in deg C for\n' + title, fontsize=12)
    # Customize z axis.
    ax.set_xlabel('x in m')
    ax.set_ylabel('z in m')
    # Add a color bar which maps values to colors.
    ticks = np.linspace(a.min(), a.max(), 11)
    fig.colorbar(heatmap, ticks=ticks, orientation='vertical', shrink=0.75,
                 aspect=20)
    # Equal gridsize.
    plt.gca().set_aspect('equal', adjustable='box')
    # Save plot to file.
    print('Save figure to {}.'.format(filepath))
    plt.savefig(filepath)
    plt.gcf().clear()
    plt.close()

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

    filepath = os.path.splitext(filepath)[0]
    title = filepath
    filepath_heatmap = filepath
    filepath_tumor = filepath
    filepath += '_py_surface.eps'
    filepath_heatmap += '_py_heatmap.eps'
    filepath_tumor += '_py_tumor.eps'

    # Create numpy array and save surface data from netCDF file to it.
    temperature = np.zeros((dim1, dim0))
    temperature[:,:] = T[(time-1):time,(dim2-1),:,:]

    plot_3d_surface(temperature, params, title, filepath)
    plot_heatmap(temperature, params, title, filepath_heatmap)

    # Create numpy array and save tumor (x-z-plane) data from netCDF file to it.
    temperature = np.zeros((dim2, dim0))
    temperature[:,:] = T[(time-1):time,:,int(dim1/2),:]

    plot_tumor(temperature, params, title, filepath_tumor)

    nc_file.close()

    print('Done.')

def main():
    print('Script will be called by startSimulation.py.')
    print('Aborting.')
    exit()

if __name__ == '__main__':
    main()
