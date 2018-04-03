import os
import sys

import netCDF4 as nc
import numpy as np
import nrrd

def save_volume_as_netcdf(data):
    nc_file = nc.Dataset('volume.nc', 'w', format='NETCDF3_CLASSIC')
    nNodes = nc_file.createDimension('nNodes_0', data.shape[0])
    nNodes = nc_file.createDimension('nNodes_1', data.shape[1])
    nNodes = nc_file.createDimension('nNodes_2', data.shape[2])
    time = nc_file.createDimension('time')
    brain = nc_file.createVariable('brain', 'f8', ('time', 'nNodes_2',
                                                  'nNodes_1', 'nNodes_0'))

    brain[0] = data.T

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

    save_volume_as_netcdf(data)

    print('Done.')

if __name__ == '__main__':
    main()
