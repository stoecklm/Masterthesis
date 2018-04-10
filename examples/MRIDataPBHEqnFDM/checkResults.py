import netCDF4 as nc
import numpy as np
import os
import sys

def temperature_array_from_result(filepath):
    print('Read data from {0}.'.format(filepath))

    nc_file = nc.Dataset(filepath)
    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size

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
        print('* ERROR: No temperature variable found in {0}.'.format(filepath))
        print('Aborting.')
        exit()

    temp = np.zeros((dim2, dim1, dim0))
    temp[:,:,:] = T[-1,:,:,:]

    nc_file.close()

    return temp

def check_results(filepath):
    a_1 = temperature_array_from_result(filepath[0])
    a_2 = temperature_array_from_result(filepath[1])
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

    print('Done.')


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

    check_results(filepath)

if __name__ == '__main__':
    main()
