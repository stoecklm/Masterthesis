import netCDF4 as nc
import numpy as np
import os
import sys

def main():
    if len(sys.argv) > 2:
        if os.path.isfile(sys.argv[1]) == True and os.path.isfile(sys.argv[2]) == True:
            filepath_1 = sys.argv[1]
            filepath_2 = sys.argv[2]
        elif os.path.isfile(sys.argv[1]) == False and os.path.isfile(sys.argv[2]) == False:
            print(sys.argv[1], 'does not exist.')
            print(sys.argv[2], 'does not exist.')
            print('Usage: python3', sys.argv[0], '<PATH/TO/FILE_1> <PATH/TO/FILE_2>')
            print('Aborting.')
            exit()
        elif os.path.isfile(sys.argv[1]) == False:
            print(sys.argv[1], 'does not exist.')
            print('Usage: python3', sys.argv[0], '<PATH/TO/FILE_1> <PATH/TO/FILE_2>')
            print('Aborting.')
            exit()
        else: #os.path.isfile(sys.argv[2]) == False:
            print(sys.argv[2], 'does not exist.')
            print('Usage: python3', sys.argv[0], '<PATH/TO/FILE_1> <PATH/TO/FILE_2>')
            print('Aborting.')
            exit()
    elif len(sys.argv) == 2:
        print('Only one command line argument for files provided.')
        print('Usage: python3', sys.argv[0], '<PATH/TO/FILE_1> <PATH/TO/FILE_2>')
        print('Aborting.')
        exit()
    else:
        print('No command line arguments for files provided.')
        print('Usage: python3', sys.argv[0], '<PATH/TO/FILE_1> <PATH/TO/FILE_2>')
        print('Aborting.')
        exit()

    nc_file_1 = nc.Dataset(filepath_1)
    nc_file_2 = nc.Dataset(filepath_2)

    dim0 = nc_file_1.dimensions['nNodes_0'].size
    dim1 = nc_file_1.dimensions['nNodes_1'].size
    dim2 = nc_file_1.dimensions['nNodes_2'].size
    time = nc_file_1.dimensions['time'].size
    TNewDom = nc_file_1.variables['TNewDom']

    a_1 = np.zeros((dim2, dim1, dim0))
    a_1[:,:,:] = TNewDom[(time-1):time,:,:,]

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

if __name__ == '__main__':
    main()
