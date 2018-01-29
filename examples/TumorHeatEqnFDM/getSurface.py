import netCDF4 as nc
import numpy as np
import os
import sys

def main():
    filepath = ''
    # Check if path to netCDF file (i.e. results) is provided,
    # if file exists and if file has .nc extension.
    if len(sys.argv) > 1:
        if os.path.isfile(sys.argv[1]) == True:
            if os.path.splitext(sys.argv[1])[1] == '.nc':
                filepath = sys.argv[1]
            else:
                print(sys.argv[1], 'does not have .nc extension.')
        else:
            print(sys.argv[1], 'does not exist.')
    else:
        print('No command line argument for netCDF file provided.')

    if filepath == '':
        print('Usage: python3', sys.argv[0], '<PATH/TO/FILE>')
        print('Aborting.')
        exit()

    print('Read data from {0}.'.format(filepath))
    # Open netCDF file and read data from it.
    nc_file = nc.Dataset(filepath)
    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size
    time = nc_file.dimensions['time'].size
    T = nc_file.variables['TNewDom']
    # Create numpy array and save surface data from netCDF file to it.
    a = np.zeros((dim1, dim0))
    a[:,:] = T[(time-1):time,(dim2-1),:,:]

    nc_file.close()

    filepath = os.path.splitext(filepath)[0]
    filepath += '_surface.dat'

    text_file = open(filepath, 'w')

    print('Write surface data to {0}.'.format(filepath))

    # Iterate through all elements of the surface.
    for elem_y in range(0, a.shape[0]):
        for elem_x in range(0, a.shape[1]):
            text_file.write('{0} {1} {2}\n'.format(str(elem_x), str(elem_y), str(a[elem_y, elem_x])))
        text_file.write('\n')

    text_file.close()

    print('Done.')

if __name__ == '__main__':
    main()
