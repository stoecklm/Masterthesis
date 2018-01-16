import netCDF4 as nc
import numpy as np
import os
import sys

def main():
    if len(sys.argv) > 1:
        if os.path.isfile(sys.argv[1]) == True:
            filepath = sys.argv[1]
        else: #os.path.isfile(sys.argv[1]) == False:
            print(sys.argv[1], 'does not exist.')
            print('Usage: python3', sys.argv[0], '<PATH/TO/FILE_1>')
            print('Aborting.')
            exit()
    else:
        print('No command line arguments for files provided.')
        print('Usage: python3', sys.argv[0], '<PATH/TO/FILE_1>')
        print('Aborting.')
        exit()

    nc_file = nc.Dataset(filepath)

    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size
    time = nc_file.dimensions['time'].size
    TNewDom = nc_file.variables['TNewDom']

    if (time < 3):
        print('Not enough timesteps in file.')
        print('Aborting.')
        exit()

    a_old = np.zeros((dim2, dim1, dim0))
    a_new = np.zeros((dim2, dim1, dim0))
    a_diff = np.zeros((dim2, dim1, dim0))

    eps_start = 1.0
    eps_end = 1.0e-15
    eps_curr = eps_start

    for step in range(time-1):
        a_old[:,:,:] = TNewDom[step:(step+1),:,:,]
        a_new[:,:,:] = TNewDom[(step+1):(step+2),:,:,]
        max_diff = np.amax(np.absolute(np.subtract(a_new, a_old)))
        while (max_diff < eps_curr):
            print('diff < {0} after {1} steps.'.format(eps_curr, step+1))
            eps_curr /= 10

    nc_file.close()

    print('Done.')

if __name__ == '__main__':
    main()
