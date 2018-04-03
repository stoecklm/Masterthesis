import os
import sys

import netCDF4 as nc
import numpy as np
import nrrd

def extract_tumor_from_volume(data, start, end):
    length = (end[0]-start[0], end[1]-start[1], end[2]-start[2])
    num_elem = length[0]*length[1]*length[2]
    print(length)
    new_data = np.arange(num_elem).reshape(length)
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
    brain = nc_file.createVariable('brain', 'f8', ('time', 'nNodes_2',
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

    print('Done.')

if __name__ == '__main__':
    main()
