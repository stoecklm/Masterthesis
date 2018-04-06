import configparser
import os
import sys

import netCDF4 as nc
import numpy as np
from numpy import linalg as LA
import nrrd
from scipy.interpolate import RegularGridInterpolator

def interpolate_3d(data, start, end, header):
    dim0_length = np.subtract(ijk_to_ras([0,0,0], header),
                              ijk_to_ras([data.shape[0]-1,0,0], header))
    dim0_length = LA.norm(dim0_length)

    dim1_length = np.subtract(ijk_to_ras([0,0,0], header),
                              ijk_to_ras([0,data.shape[1]-1,0], header))
    dim1_length = LA.norm(dim1_length)

    dim2_length = np.subtract(ijk_to_ras([0,0,0], header),
                              ijk_to_ras([0,0,data.shape[2]-1], header))
    dim2_length = LA.norm(dim2_length)
    print(dim0_length)
    print(dim1_length)
    print(dim2_length)

    x_mri = np.linspace(0, dim0_length, data.shape[0])
    y_mri = np.linspace(0, dim1_length, data.shape[1])
    z_mri = np.linspace(0, dim2_length, data.shape[2])

    dim0 = 120
    dim1 = 120
    dim2 = 50

    x_scafes = np.linspace(0, 120, dim0)
    y_scafes = np.linspace(0, 120, dim1)
    z_scafes = np.linspace(0, 60, dim2)

    my_interpolating_function = RegularGridInterpolator((x_mri, y_mri, z_mri),
                                                        data,
                                                        bounds_error=False,
                                                        fill_value=0)

    x, y, z = np.meshgrid(x_scafes, y_scafes, z_scafes, sparse=True,
                          indexing='ij')

    new_data = np.zeros(dim0*dim1*dim2).reshape((dim0, dim1, dim2))
    for elem_x in range(0, new_data.shape[0]):
        for elem_y in range(0, new_data.shape[1]):
            for elem_z in range(0, new_data.shape[2]):
                tmp = my_interpolating_function([x_scafes[elem_x],
                                                 y_scafes[elem_y],
                                                 z_scafes[elem_z]])
                new_data[elem_x, elem_y, elem_z] = tmp

    return new_data

def stats_from_tumor(data):
    print('Mean: {}.'.format(np.mean(data)))
    print('Min: {}.'.format(np.min(data)))
    print('Max: {}.'.format(np.max(data)))
    print('StdDev: {}.'.format(np.std(data)))

def binary_data(data):
    bin_data = np.zeros(data.shape)
    mean = np.mean(data)
    std_dev = np.std(data)

    bin_data = np.where(data > mean + std_dev, 1, 0)

    return bin_data

def extract_tumor_from_volume(data, start, end):
    new_data = data[start[0]:end[0]+1,
                    start[1]:end[1]+1,
                    start[2]:end[2]+1]

    return new_data

def save_volume_as_netcdf(data, filename):
    nc_file = nc.Dataset(filename, 'w', format='NETCDF3_CLASSIC')
    nNodes = nc_file.createDimension('nNodes_0', data.shape[0])
    nNodes = nc_file.createDimension('nNodes_1', data.shape[1])
    nNodes = nc_file.createDimension('nNodes_2', data.shape[2])
    time = nc_file.createDimension('time')
    brain = nc_file.createVariable('region', 'i2', ('time', 'nNodes_2',
                                                  'nNodes_1', 'nNodes_0'))
    brain[0,:,:,:] = data.T

    nc_file.close()

def read_volume(filepath):
    data, header = nrrd.read(filepath)

    return data, header

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

def bounding_box_tumor(start, end):
    bbox = np.asarray([[start[0], start[1], start[2]],
                       [end[0], start[1], start[2]],
                       [end[0], end[1], start[2]],
                       [start[0], end[1], start[2]],
                       [start[0], start[1], end[2]],
                       [end[0], start[1], end[2]],
                       [end[0], end[1], end[2]],
                       [start[0], end[1], end[2]]])

    return bbox

def main():
    filepath = ''
    if len(sys.argv) > 1:
        if os.path.isdir(sys.argv[1]) == True:
            tmp = os.path.join(os.sys.argv[1], 'volume.nrrd')
            if os.path.isfile(tmp) != True:
                print(sys.argv[1], 'does not contain volume.nrrd.')
                print('Aborting.')
                exit()
            else:
                folderpath = sys.argv[1]
                filepath = tmp
        else:
            print(sys.argv[1], 'does not exist.')
            print('Aborting.')
            exit()
        tmp = os.path.join(folderpath, 'bounding-box.ini')
        if os.path.isfile(tmp) != True:
            print(sys.argv[1], 'does not contain bounding-box.ini.')
            print('Aborting.')
            exit()
        else:
            bbox = tmp
    else:
        print('No command line argument for folder provided.')
        print('Aborting.')
        exit()

    data, header = read_volume(filepath)

    save_volume_as_netcdf(data, 'whole_mri_volume.nc')

    #point = (nNodes_0, nNodes_1, nNodes_2)
    space_origin = [0,0,0]
    space_end = list(map(float, header['space origin']))
    config = configparser.ConfigParser()
    config.read(bbox)

    start = config['Start'].get('START')
    start = np.asarray(list(map(int, start.split('x'))))
    end = config['End'].get('END')
    end = np.asarray(list(map(int, end.split('x'))))

    tumor = extract_tumor_from_volume(data, start, end)
    save_volume_as_netcdf(tumor, 'tumor_mri.nc')

    binary_tumor = binary_data(tumor)
    save_volume_as_netcdf(binary_tumor, 'bin_data.nc')

    #interpolated_tumor = interpolate_3d(binary_tumor, start, end, header)
    #save_volume_as_netcdf(interpolated_tumor, 'region.nc')

    print('Done.')

if __name__ == '__main__':
    main()
