import netCDF4 as nc
import numpy as np
import os
import subprocess
import sys

def main():
    global params
    # Check if path to configfile is provided and
    # if file exists.
    if len(sys.argv) > 1:
        if os.path.isfile(sys.argv[1]) == True:
            file_name = sys.argv[1]
        else:
            print(sys.argv[1], 'does not exist.')
            print('Usage: python3', sys.argv[0], '<PATH/TO/FILE>')
            print('Aborting.')
            exit()
    else:
        print('No command line argument for netCDF file provided.')
        print('Usage: python3', sys.argv[0], '<PATH/TO/FILE>')
        print('Aborting.')
        exit()

    nc_file = nc.Dataset(file_name)
    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size
    time = nc_file.dimensions['time'].size
    T = nc_file.variables['TNewDom']

    a = np.zeros((dim1, dim0))
    a[:,:] = T[(time-1):time,(dim2-1),:,:]

    nc_file.close()

    text_file = open('surface.dat', 'w')

    for elem_y in range(0, a.shape[0]):
        for elem_x in range(0, a.shape[1]):
            text_file.write('{0} {1} {2}\n'.format(str(elem_x), str(elem_y), str(a[elem_y, elem_x])))
        text_file.write('\n')

if __name__ == '__main__':
    main()
