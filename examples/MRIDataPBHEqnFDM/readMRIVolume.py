import configparser
import os
from shutil import copyfile
import sys

import netCDF4 as nc
import numpy as np
from numpy import linalg as LA
import nrrd
from scipy.interpolate import RegularGridInterpolator

from readMRIData import plot_lin_plane_fitting_with_bbox
from readMRIData import read_intra_op_points

def read_nrrd_file(filepath):
    data, header = nrrd.read(filepath)

    return data, header

def save_as_netcdf(data, filename):
    nc_file = nc.Dataset(filename, 'w', format='NETCDF3_CLASSIC')
    nc_file.createDimension('nNodes_0', data.shape[0])
    nc_file.createDimension('nNodes_1', data.shape[1])
    nc_file.createDimension('nNodes_2', data.shape[2])
    nc_file.createDimension('time')
    brain = nc_file.createVariable('region', 'i2', ('time', 'nNodes_2',
                                                    'nNodes_1', 'nNodes_0'))
    brain[0,] = np.swapaxes(data, 0, 2)
    nc_file.close()

def data_as_binary_data(data):
    binary_data = np.zeros(data.shape)
    binary_data = np.where(data == 255, 1, 0)

    return binary_data

def get_ijk_to_lps(header):
    space_origin = list(map(float, header['space origin']))
    space_origin.append(1)

    row0 = list(map(float, header['space directions'][0]))
    row0.append(0)
    row1 = list(map(float, header['space directions'][1]))
    row1.append(0)
    row2 = list(map(float, header['space directions'][2]))
    row2.append(0)
    ijk_to_lps = np.array([row0, row1, row2, space_origin]).T

    return ijk_to_lps

def ijk_to_ras(ijk, header):
    if len(ijk) == 3:
        ijk = np.append(ijk, 1)

    ijk_to_lps = get_ijk_to_lps(header)
    lps_to_ras = np.diag([-1, -1, 1, 1])

    lps = np.dot(ijk_to_lps, ijk)
    ras = np.matmul(lps, lps_to_ras)

    return ras[0:3]

def ras_to_ijk(ras, header):
    if len(ras) == 3:
        ras = np.append(ras, 1)

    ijk_to_lps = get_ijk_to_lps(header)
    lps_to_ijk = LA.inv(ijk_to_lps)
    ras_to_lps = np.diag([-1, -1, 1, 1])

    lps = np.matmul(ras, ras_to_lps)
    ijk = np.dot(lps_to_ijk, lps)

    return ijk[0:3]

def volume_dimensions(header):
    start = np.asarray([0, 0, 0])
    end = np.asarray(header['sizes'])

    dim0_length = np.subtract(ijk_to_ras(start, header),
                              ijk_to_ras([end[0],start[1],start[2]], header))
    dim0_length = LA.norm(dim0_length)

    dim1_length = np.subtract(ijk_to_ras(start, header),
                              ijk_to_ras([start[0],end[1],start[2]], header))
    dim1_length = LA.norm(dim1_length)

    dim2_length = np.subtract(ijk_to_ras(start, header),
                              ijk_to_ras([start[0],start[1],end[2]], header))
    dim2_length = LA.norm(dim2_length)

    return dim0_length, dim1_length, dim2_length

def interpolate_data(data, header):
    dim0_length, dim1_length, dim2_length = volume_dimensions(header)

    x_mri = np.linspace(0, dim0_length, data.shape[0])
    y_mri = np.linspace(0, dim1_length, data.shape[1])
    z_mri = np.linspace(0, dim2_length, data.shape[2])

    dim0_gridsize = 120/(120-1)
    dim1_gridsize = 120/(120-1)
    dim2_gridsize = 60/(50-1)

    dim0 = int(dim0_length/dim0_gridsize)+1
    dim1 = int(dim1_length/dim1_gridsize)+1
    dim2 = int(dim2_length/dim2_gridsize)+1

    x_scafes = np.linspace(0, dim0_length, dim0)
    y_scafes = np.linspace(0, dim1_length, dim1)
    z_scafes = np.linspace(0, dim2_length, dim2)

    my_interpolating_function = RegularGridInterpolator((x_mri, y_mri, z_mri),
                                                        data, method='nearest',
                                                        bounds_error=False,
                                                        fill_value=0)

    new_data = np.zeros(dim0*dim1*dim2).reshape((dim0, dim1, dim2))

    for elem_x in range(0, new_data.shape[0]):
        for elem_y in range(0, new_data.shape[1]):
            for elem_z in range(0, new_data.shape[2]):
                tmp = my_interpolating_function([x_scafes[elem_x],
                                                 y_scafes[elem_y],
                                                 z_scafes[elem_z]])
                new_data[elem_x, elem_y, elem_z] = tmp

    return new_data

def return_grid(data, header, case):
    dim0, dim1, dim2 = data.shape

    new_dim0 = 120
    new_dim1 = 120
    new_dim2 = 50

    num_elem = new_dim0 * new_dim1 * new_dim2
    new_data = np.zeros(num_elem).reshape((new_dim0, new_dim1, new_dim2))

    i = int((new_dim0-dim0)/2)
    j = int((new_dim1-dim1)/2)
    #k = int((new_dim2-dim2)/2)
    k = 5

    new_data[i:i+dim0, j:j+dim1, k:k+dim2] = data
    #new_data = np.flip(new_data, 2)
    new_data = new_data[:,:,::-1]

    return new_data

def bounding_box(start, end):
    bbox = np.asarray([[start[0], start[1], start[2]],
                       [end[0], start[1], start[2]],
                       [end[0], end[1], start[2]],
                       [start[0], end[1], start[2]],
                       [start[0], start[1], end[2]],
                       [end[0], start[1], end[2]],
                       [end[0], end[1], end[2]],
                       [start[0], end[1], end[2]]])

    return bbox

def write_ini_file(case):
    print('Write {}.ini.'.format(case))

    config = configparser.ConfigParser()
    config.optionxform = str

    if os.path.isfile('Parameters.ini') == True:
        config.read('Parameters.ini')
    else:
        print('Parameters.ini does not exist.')
        print('Aborting.')
        exit()

    config['MRI']['CASE'] = case
    config['Input']['NAME_REGION_FILE'] = case + '_region'
    config['Input']['USE_MRI_FILE'] = 'True'
    config['Input']['USE_INITFILE'] = 'True'
    config['Input']['CREATE_INITFILE'] = 'True'

    with open(case + '.ini', 'w') as configfile:
        config.write(configfile)

    print('Done.')


def main():
    if len(sys.argv) > 1:
        if os.path.isdir(sys.argv[1]) == True:
            tmp = os.path.join(os.sys.argv[1], '_1.nrrd')
            if os.path.isfile(tmp) != True:
                print(sys.argv[1], 'does not contain _1.nrrd.')
                print('Aborting.')
                exit()
            else:
                folderpath = sys.argv[1]
                filepath = tmp
        else:
            print(sys.argv[1], 'does not exist.')
            print('Aborting.')
            exit()
    else:
        print('No command line argument for folder provided.')
        print('Aborting.')
        exit()

    print('Create region file.')

    case = filepath.split('_')[0]

    data, header = read_nrrd_file(filepath)
    save_as_netcdf(data, case + '_seg_tumor.nc')

    binary_data = data_as_binary_data(data)
    save_as_netcdf(binary_data, case + '_binary_tumor.nc')

    interpolated_data = interpolate_data(binary_data, header)
    save_as_netcdf(interpolated_data, case + '_interpolated_data.nc')

    region = return_grid(interpolated_data, header, case)
    save_as_netcdf(region, case + '_region.nc')

    iop = read_intra_op_points(folderpath)
    start = np.asarray([0, 0, 0])
    end = np.asarray(header['sizes'])
    bbox = bounding_box(ijk_to_ras(start, header), ijk_to_ras(end, header))
    plot_lin_plane_fitting_with_bbox(iop, bbox, case + '_domain')

    write_ini_file(case)

    print('Done.')

if __name__ == '__main__':
    main()
