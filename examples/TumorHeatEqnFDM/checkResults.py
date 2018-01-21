import netCDF4 as nc
import numpy as np
import os
import sys

def main():
    filepath = []
    # Check if paths to netCDF files (i.e. results) are provided,
    # if files exist and if files have .nc extension.
    if len(sys.argv) > 2:
        for elem in range(1, 3):
            if os.path.isfile(sys.argv[elem]) == True:
                if os.path.splitext(sys.argv[elem])[1] == '.nc':
                    filepath.append(sys.argv[elem])
                else:
                    print(sys.argv[elem], 'does not have .nc extension.')
            else:
                print(sys.argv[elem], 'does not exist.')
    else:
        print('Not enough command line arguments for netCDF files provided.')

    if len(filepath) < 2:
        print('Usage: python3', sys.argv[0], '<PATH/TO/FILE_1> <PATH/TO/FILE_2>')
        print('Aborting.')
        exit()

    print('Read data from {0}.'.format(filepath[0]))

    nc_file_1 = nc.Dataset(filepath[0])
    dim0 = nc_file_1.dimensions['nNodes_0'].size
    dim1 = nc_file_1.dimensions['nNodes_1'].size
    dim2 = nc_file_1.dimensions['nNodes_2'].size
    time = nc_file_1.dimensions['time'].size
    TNewDom = nc_file_1.variables['TNewDom']

    a_1 = np.zeros((dim2, dim1, dim0))
    a_1[:,:,:] = TNewDom[(time-1):time,:,:,]

    print('Read data from {0}.'.format(filepath[1]))

    nc_file_2 = nc.Dataset(filepath[1])
    dim0 = nc_file_2.dimensions['nNodes_0'].size
    dim1 = nc_file_2.dimensions['nNodes_1'].size
    dim2 = nc_file_2.dimensions['nNodes_2'].size
    time = nc_file_2.dimensions['time'].size
    TNewDom = nc_file_2.variables['TNewDom']

    a_2 = np.zeros((dim2, dim1, dim0))
    a_2[:,:,:] = TNewDom[(time-1):time,:,:,]

    print()
    if np.array_equal(a_1, a_2) == True:
        print('SUCCESS: Last timestep is equal.')
    else:
        if a_1.shape[0] != a_2.shape[0] or \
           a_1.shape[1] != a_2.shape[1] or \
           a_1.shape[2] != a_2.shape[2]:
            print('Shape of files is not equal.')
        else:
            print('WARNING: Last timestep is NOT equal.')

    nc_file_1.close()
    nc_file_2.close()

    print('Done.')

if __name__ == '__main__':
    main()
