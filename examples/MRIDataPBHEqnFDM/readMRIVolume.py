import os
import sys

import netCDF4 as nc
import numpy as np
import nrrd
from scipy.interpolate import RegularGridInterpolator

def interpolate_3d(data):
    x_mri = np.linspace(0, 1, data.shape[0])
    y_mri = np.linspace(0, 1, data.shape[1])
    z_mri = np.linspace(0, 1, data.shape[2])

    dim0 = 120
    dim1 = 120
    dim2 = 50

    x_scafes = np.linspace(0, 1, dim0)
    y_scafes = np.linspace(0, 1, dim1)
    z_scafes = np.linspace(0, 1, dim2)

    my_interpolating_function = RegularGridInterpolator((x_mri, y_mri, z_mri), data)

    x, y, z = np.meshgrid(x_scafes, y_scafes, z_scafes, sparse=True, indexing='ij')

    new_data = np.zeros(dim0*dim1*dim2).reshape((dim0, dim1, dim2))

    print(new_data.shape)

    for elem_x in range(0, new_data.shape[0]):
        for elem_y in range(0, new_data.shape[1]):
            for elem_z in range(0, new_data.shape[2]):
                new_data[elem_x, elem_y, elem_z] = my_interpolating_function([x_scafes[elem_x],
                                                                              y_scafes[elem_y],
                                                                              z_scafes[elem_z]])

    save_volume_as_netcdf(new_data, 'region.nc')


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

    save_volume_as_netcdf(bin_data, 'bin_data.nc')

    return bin_data

def extract_tumor_from_volume(data, start, end):
    length = (end[0]-start[0]+1, end[1]-start[1]+1, end[2]-start[2]+1)
    num_elem = length[0]*length[1]*length[2]
    print(length)
    new_data = np.arange(num_elem).reshape(length)
    print(new_data.shape)
    new_data = data[start[0]:end[0]+1,
                    start[1]:end[1]+1,
                    start[2]:end[2]+1]
    print(new_data.shape)
    print('....')

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

    print(data.shape)
    print(header)

    return data

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
                filepath = tmp
        else:
            print(sys.argv[1], 'does not exist.')
            print('Aborting.')
            exit()
    else:
        print('No command line argument for folder provided.')
        print('Aborting.')
        exit()

    data = read_volume(filepath)

    save_volume_as_netcdf(data, 'whole_mri_volume.nc')

    #start = (nNodes_0, nNodes_1, nNodes_2)
    start = (90,45,20)
    end = (150,100,80)
    # end = (224, 256, 176)
    tumor = extract_tumor_from_volume(data, start, end)

    save_volume_as_netcdf(tumor, 'tumor_mri.nc')

    stats_from_tumor(tumor)

    binary_tumor = binary_data(tumor)

    interpolate_3d(binary_tumor)

    print('Done.')

if __name__ == '__main__':
    main()
