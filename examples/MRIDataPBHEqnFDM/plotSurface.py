import os
import sys

import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
import numpy as np
import numpy.ma as ma

from helperFunctions import temperature_array_from_result
from postProcessing import surface_array_from_file

CMAP = cm.viridis

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
    #fig.suptitle('Surface Temperature in deg C for\n' + title, fontsize=12)
    # Plot surface.
    surf = ax.plot_surface(x, y, a, cmap=CMAP,
                           linewidth=0, antialiased=False)
    # Customize z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    # Add a color bar which maps values to colors.
    ticks = np.linspace(a.min(), a.max(), 11)
    cbar = fig.colorbar(surf, ticks=ticks, orientation='vertical', shrink=0.75,
                        aspect=20)
    cbar.set_label('\nTemperature in °C')
    # Invert z axis since the relevant part is colder than the other part
    # and therefore hard to see.
    ax.invert_zaxis()
    # Save plot to file.
    print('Save figure to {}.'.format(filepath))
    plt.savefig(filepath, bbox_inches='tight')
    plt.gcf().clear()
    plt.close()

def plot_heatmap(a, params, title, filepath):
    COORD_NODE_FIRST = params['COORD_NODE_FIRST']
    COORD_NODE_LAST = params['COORD_NODE_LAST']
    DIM = params['N_NODES']
    TUMOR_CENTER = params['TUMOR_CENTER']
    RADIUS = params['PARAMETERS']['DIAMETER']/2
    x, y = np.meshgrid(np.linspace(COORD_NODE_FIRST[0], COORD_NODE_LAST[0],
                                   DIM[0]),
                       np.linspace(COORD_NODE_FIRST[1], COORD_NODE_LAST[1],
                                   DIM[1]))
    # Plot heatmap with circle around hole.
    fig, ax = plt.subplots()
    heatmap = ax.pcolormesh(x, y, a, cmap=CMAP, rasterized=True)
    try:
        pts = params['HOLE']
        ax.plot(pts[:,0], pts[:,1], color='r', linestyle='dashed')
    except KeyError:
        circle = plt.Circle((TUMOR_CENTER[0], TUMOR_CENTER[1]), RADIUS,
                            color='r', fill=False, linestyle='dashed')
        ax.add_artist(circle)
    # Title.
    #fig.suptitle('Heatmap in deg C for\n' + title, fontsize=12)
    # Label for axis.
    ax.set_xlabel('x in m')
    ax.set_ylabel('y in m')
    # Add a color bar which maps values to colors.
    ticks = np.linspace(a.min(), a.max(), 11)
    cbar = fig.colorbar(heatmap, ticks=ticks, orientation='vertical',
                        shrink=0.75, aspect=20)
    cbar.set_label('\nTemperature in °C')
    # Equal gridsize.
    plt.gca().set_aspect('equal', adjustable='box')
    # Save plot to file.
    print('Save figure to {}.'.format(filepath))
    plt.savefig(filepath, bbox_inches='tight')
    plt.gcf().clear()
    plt.close()

def plot_heatmap_scaled(temp, params, title, filepath):
    surface = surface_array_from_file(params['NAME_INITFILE'])
    skull = surface[-1,:,:]
    if np.count_nonzero(skull == 1) != 0:
        temp[np.where(skull == 0)] = float('nan')
        temp = ma.masked_where(np.isnan(temp), temp)
    else:
        print('No open surface specified.')

    COORD_NODE_FIRST = params['COORD_NODE_FIRST']
    COORD_NODE_LAST = params['COORD_NODE_LAST']
    DIM = params['N_NODES']
    TUMOR_CENTER = params['TUMOR_CENTER']
    RADIUS = params['PARAMETERS']['DIAMETER']/2
    x, y = np.meshgrid(np.linspace(COORD_NODE_FIRST[0], COORD_NODE_LAST[0],
                                   DIM[0]),
                       np.linspace(COORD_NODE_FIRST[1], COORD_NODE_LAST[1],
                                   DIM[1]))
    fig, ax = plt.subplots()
    heatmap = ax.pcolormesh(x, y, temp, cmap=CMAP, rasterized=True)
    #fig.suptitle('Heatmap in deg C for\n' + title, fontsize=12)
    ax.set_xlabel('x in m')
    ax.set_ylabel('y in m')
    ticks = np.linspace(temp.min(), temp.max(), 11)
    cbar = fig.colorbar(heatmap, ticks=ticks, orientation='vertical',
                        shrink=0.75, aspect=20)
    cbar.set_label('\nTemperature in °C')
    plt.gca().set_aspect('equal', adjustable='box')
    print('Save figure to {}.'.format(filepath))
    plt.savefig(filepath, bbox_inches='tight')
    plt.gcf().clear()
    plt.close()

def plot_heatmap_thermo_scaled(temp, params, title, filepath, csv_min, csv_max):
    TICKS = np.linspace(csv_min, csv_max, 10)
    LEVELS = MaxNLocator(nbins=1000).tick_values(csv_min, csv_max)
    NORM = BoundaryNorm(LEVELS, ncolors=CMAP.N, clip=True)
    surface = surface_array_from_file(params['NAME_INITFILE'])
    skull = surface[-1,:,:]
    if np.count_nonzero(skull == 1) != 0:
        temp[np.where(skull == 0)] = float('nan')
        temp = ma.masked_where(np.isnan(temp), temp)
    else:
        print('No open surface specified.')

    COORD_NODE_FIRST = params['COORD_NODE_FIRST']
    COORD_NODE_LAST = params['COORD_NODE_LAST']
    DIM = params['N_NODES']
    TUMOR_CENTER = params['TUMOR_CENTER']
    RADIUS = params['PARAMETERS']['DIAMETER']/2
    x, y = np.meshgrid(np.linspace(COORD_NODE_FIRST[0], COORD_NODE_LAST[0],
                                   DIM[0]),
                       np.linspace(COORD_NODE_FIRST[1], COORD_NODE_LAST[1],
                                   DIM[1]))
    fig, ax = plt.subplots()
    heatmap = ax.pcolormesh(x, y, temp, cmap=CMAP, norm=NORM, rasterized=True)
    ax.set_xlabel('x in m')
    ax.set_ylabel('y in m')
    cbar = fig.colorbar(heatmap, ticks=TICKS, orientation='vertical',
                        shrink=0.75, aspect=20)
    cbar.set_label('\nTemperature in °C')
    plt.gca().set_aspect('equal', adjustable='box')
    print('Save figure to {}.'.format(filepath))
    plt.savefig(filepath, bbox_inches='tight')
    plt.gcf().clear()
    plt.close()

def plot_tumor(a, params, title, filepath):
    COORD_NODE_FIRST = params['COORD_NODE_FIRST']
    COORD_NODE_LAST = params['COORD_NODE_LAST']
    DIM = params['N_NODES']
    TUMOR_CENTER = params['TUMOR_CENTER']
    RADIUS = params['PARAMETERS']['DIAMETER']/2
    x, z = np.meshgrid(np.linspace(COORD_NODE_FIRST[0], COORD_NODE_LAST[0],
                                   DIM[0]),
                       np.linspace(COORD_NODE_FIRST[2], COORD_NODE_LAST[2],
                                   DIM[2]))
    # Plot heatmap with circle around tumor.
    fig, ax = plt.subplots()
    heatmap = ax.pcolormesh(x, z, a, cmap=CMAP, rasterized=True)
    if params['USE_MRI_FILE'] == False:
        circle = plt.Circle((TUMOR_CENTER[0], TUMOR_CENTER[2]), RADIUS, color='r',
                             fill=False, linestyle='dashed')
        ax.add_artist(circle)
    # Title.
    #fig.suptitle('Heatmap in deg C for\n' + title, fontsize=12)
    # Customize z axis.
    ax.set_xlabel('x in m')
    ax.set_ylabel('z in m')
    # Add a color bar which maps values to colors.
    ticks = np.linspace(a.min(), a.max(), 11)
    cbar = fig.colorbar(heatmap, ticks=ticks, orientation='vertical',
                        shrink=0.75, aspect=20)
    cbar.set_label('\nTemperature in °C')
    # Equal gridsize.
    plt.gca().set_aspect('equal', adjustable='box')
    # Save plot to file.
    print('Save figure to {}.'.format(filepath))
    plt.savefig(filepath, bbox_inches='tight')
    plt.gcf().clear()
    plt.close()

def plot_thermo(case, folder):
    filepath = case + '_thermo.eps'
    a = np.genfromtxt(os.path.join(folder, 'thermo.csv'), delimiter=',')
    a[np.isnan(a)] = 0
    rows = np.any(a, axis=1)
    cols = np.any(a, axis=0)
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]
    b = np.zeros((rmax-rmin+1, cmax-cmin+1))
    b[:,:] = a[rmin:rmax+1,cmin:cmax+1]
    b[b == 0] = float('nan')
    x, y = np.meshgrid(np.linspace(0, b.shape[1]-1, b.shape[1]), \
                       np.linspace(0, b.shape[0]-1, b.shape[0]))
    fig, ax = plt.subplots()
    c = ma.masked_where(np.isnan(b), b)
    ticks = np.linspace(c.min(), c.max(), 11)
    heatmap = ax.pcolormesh(x, y, c, cmap=CMAP, rasterized=True)
    cbar = fig.colorbar(heatmap, ticks=ticks, orientation='vertical',
                        shrink=0.75, aspect=20)
    cbar.set_label('\nTemperature in °C')
    # Equal gridsize.
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tick_params(axis='both', which='both', bottom='off', top='off',
                    labelbottom='off', right='off', left='off', labelleft='off')
    # Save plot to file.
    print('Save figure to {}.'.format(filepath))
    plt.savefig(filepath, bbox_inches='tight')
    plt.gcf().clear()
    plt.close()

    return c.min(), c.max()

def plot_surface(filepath, params):
    print('Plotting {0}.'.format(filepath))

    T = temperature_array_from_result(filepath)
    dim2, dim1, dim0 = T.shape

    filepath = os.path.splitext(filepath)[0]
    title = filepath
    filepath_surface = filepath + '_surface.eps'
    filepath_heatmap = filepath + '_heatmap.eps'
    filepath_heatmap_scaled = filepath + '_heatmap_scaled.eps'
    filepath_heatmap_thermo_scaled = filepath + '_heatmap_thermo_scaled.eps'
    filepath_tumor = filepath + '_tumor.eps'

    # Surface data.
    temperature = np.zeros((dim1, dim0))
    temperature[:,:] = T[-1,:,:]
    plot_3d_surface(temperature, params, title, filepath_surface)
    plot_heatmap(temperature, params, title, filepath_heatmap)
    plot_heatmap_scaled(temperature, params, title, filepath_heatmap_scaled)
    if params['MRI_DATA_CASE'] != '':
        csv_min, csv_max = plot_thermo(params['MRI_DATA_CASE'],
                                       params['MRI_DATA_FOLDER'])
        plot_heatmap_thermo_scaled(temperature, params, title,
                                   filepath_heatmap_thermo_scaled, csv_min,
                                   csv_max)

    # x-z-plane data.
    temperature = np.zeros((dim2, dim0))
    temperature[:,:] = T[:,int(dim1/2),:]
    plot_tumor(temperature, params, title, filepath_tumor)


    print('Done.')

def main():
    print('Script will be called by startSimulation.py.')
    print('Aborting.')
    exit()

if __name__ == '__main__':
    main()
