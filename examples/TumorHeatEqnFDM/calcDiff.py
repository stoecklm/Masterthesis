import numpy as np
import os
from scipy.interpolate import griddata
import sys

def two_netcdf_files_surface(filepath, a_1, a_2):
    import netCDF4 as nc
    filepath = os.path.splitext(filepath[0])[0] + '_' \
               + os.path.basename(os.path.splitext(filepath[1])[0]) \
               + '_diff_surface.nc'
    print('Write surface data to {0}.'.format(filepath))
    if a_1.shape[2] == a_2.shape[2] and \
       a_1.shape[1] == a_2.shape[1]:
        dim0_sparse = a_1.shape[2]
        dim1_sparse = a_1.shape[1]
        a_3 = np.subtract(a_1[a_1.shape[0]-1,:,:], a_2[a_2.shape[0]-1,:,:])
    elif a_1.shape[2] != a_2.shape[2] and \
         a_1.shape[1] != a_2.shape[1]:
        print('* WARNING: Shape in x dim and y dim is not equal.')
        print('  Using downsampling and cubic interpolation to calc diff',
              'between results.')
        if a_1.shape[1] < a_2.shape[1]:
            # a_1 = sparse
            # a_2 = dense
            dim0_sparse = a_1.shape[2]
            dim1_sparse = a_1.shape[1]
            dim0_dense = a_2.shape[2]
            dim1_dense = a_2.shape[1]
        else:
            # a_1 = dense
            # a_2 = sparse
            dim0_sparse = a_2.shape[2]
            dim1_sparse = a_2.shape[1]
            dim0_dense = a_1.shape[2]
            dim1_dense = a_1.shape[1]

        x_sparse, y_sparse = np.meshgrid(np.linspace(0, dim0_sparse-1,
                                                     dim0_sparse),
                                         np.linspace(0, dim1_sparse-1,
                                                     dim1_sparse))
        x_dense, y_dense = np.meshgrid(np.linspace(0, dim0_sparse-1,
                                                   dim0_dense),
                                       np.linspace(0, dim1_sparse-1,
                                                   dim1_dense))
        if a_1.shape[1] < a_2.shape[1]:
            # a_1 = sparse
            # a_2 = dense
            a_sparse = griddata(np.array([x_dense.ravel(),
                                          y_dense.ravel()]).T,
                                a_2[a_2.shape[0]-1,:,:].ravel(),
                                (x_sparse, y_sparse), method='cubic')
            a_3 = np.subtract(a_1[a_1.shape[0]-1,:,:], a_sparse[:,:])
        else:
            # a_1 = dense
            # a_2 = sparse
            a_sparse = griddata(np.array([x_dense.ravel(), y_dense.ravel()]).T,
                                a_1[a_1.shape[0]-1,:,:].ravel(),
                                (x_sparse, y_sparse), method='cubic')
            a_3 = np.subtract(a_sparse[:,:], a_2[a_2.shape[0]-1,:,:])

    nc_file = nc.Dataset(filepath, 'w', format='NETCDF3_CLASSIC')
    nNodes = nc_file.createDimension('nNodes_0', dim0_sparse)
    nNodes = nc_file.createDimension('nNodes_1', dim1_sparse)
    time = nc_file.createDimension('time')
    diff = nc_file.createVariable('TDiff', 'f8', ('time', 'nNodes_1',
                                                  'nNodes_0'))

    diff[0,] = a_3

    nc_file.close()

    max_diff = np.amax(np.absolute(a_3))
    sum_diff = np.mean(np.absolute(a_3))
    print('Surface:  max(abs(diff)) = {0}.'.format(max_diff))
    print('Surface: mean(abs(diff)) = {0}.'.format(sum_diff))

    filepath = os.path.splitext(filepath)[0]
    filepath += '.dat'
    text_file = open(filepath, 'w')

    print('Write surface data to {0}.'.format(filepath))

    for elem_y in range(0, a_3.shape[0]):
        for elem_x in range(0, a_3.shape[1]):
            text_file.write('{0} {1} {2}\n'.format(str(elem_x), str(elem_y),
                                                   str(a_3[elem_y, elem_x])))
        text_file.write('\n')

    text_file.close()

def two_netcdf_files_volume(filepath, a_1, a_2, dim0, dim1, dim2):
    import netCDF4 as nc
    filepath = os.path.splitext(filepath[0])[0] + '_' \
               + os.path.basename(os.path.splitext(filepath[1])[0]) \
               + '_diff_volume.nc'
    print('Write volume data to {0}.'.format(filepath))

    nc_file = nc.Dataset(filepath, 'w', format='NETCDF3_CLASSIC')
    nNodes = nc_file.createDimension('nNodes_0', dim0)
    nNodes = nc_file.createDimension('nNodes_1', dim1)
    nNodes = nc_file.createDimension('nNodes_2', dim2)
    time = nc_file.createDimension('time')
    diff = nc_file.createVariable('TDiff', 'f8', ('time', 'nNodes_2',
                                                  'nNodes_1', 'nNodes_0'))
    a_3 = np.subtract(a_1, a_2)
    diff[0,] = a_3

    nc_file.close()

    max_diff = np.amax(np.absolute(a_3))
    max_diff_index = np.unravel_index(np.argmax(np.absolute(a_3)), a_3.shape)
    max_diff_index = max_diff_index[::-1]
    sum_diff = np.mean(np.absolute(a_3))
    print('Volume:   max(abs(diff)) = {0} at {1}.'.format(max_diff,
                                                          max_diff_index))
    print('Volume:  mean(abs(diff)) = {0}.'.format(sum_diff))
    max_diff = np.amax(np.absolute(a_3[dim2-1,:,:]))
    max_diff_index = np.unravel_index(np.argmax(np.absolute(a_3[dim2-1,:,:])),
                                      a_3.shape)
    max_diff_index = max_diff_index[::-1]
    max_diff_index = max_diff_index[0], max_diff_index[1], dim2-1
    sum_diff = np.mean(np.absolute(a_3[dim2-1,:,:]))
    print('Surface:  max(abs(diff)) = {0} at {1}.'.format(max_diff,
                                                          max_diff_index))
    print('Surface: mean(abs(diff)) = {0}.'.format(sum_diff))

    filepath = os.path.splitext(filepath)[0]
    filepath += '_surface.dat'
    text_file = open(filepath, 'w')

    print('Write surface data to {0}.'.format(filepath))

    a = np.zeros((dim1, dim0))
    a[:,:] = a_3[dim2-1,:,:]

    for elem_y in range(0, a.shape[0]):
        for elem_x in range(0, a.shape[1]):
            text_file.write('{0} {1} {2}\n'.format(str(elem_x), str(elem_y),
                                                   str(a[elem_y, elem_x])))
        text_file.write('\n')

    text_file.close()

def two_netcdf_files(filepath):
    import netCDF4 as nc
    print('Read data from {0}.'.format(filepath[0]))

    nc_file_1 = nc.Dataset(filepath[0])
    dim0 = nc_file_1.dimensions['nNodes_0'].size
    dim1 = nc_file_1.dimensions['nNodes_1'].size
    dim2 = nc_file_1.dimensions['nNodes_2'].size
    time = nc_file_1.dimensions['time'].size
    TNewDom = nc_file_1.variables['TNewDom']

    a_1 = np.zeros((dim2, dim1, dim0))
    a_1[:,:,:] = TNewDom[(time-1):time,:,:,]

    nc_file_1.close()

    print('Read data from {0}.'.format(filepath[1]))

    nc_file_2 = nc.Dataset(filepath[1])
    dim0 = nc_file_2.dimensions['nNodes_0'].size
    dim1 = nc_file_2.dimensions['nNodes_1'].size
    dim2 = nc_file_2.dimensions['nNodes_2'].size
    time = nc_file_2.dimensions['time'].size
    TNewDom = nc_file_2.variables['TNewDom']

    a_2 = np.zeros((dim2, dim1, dim0))
    a_2[:,:,:] = TNewDom[(time-1):time,:,:,]

    nc_file_2.close()

    if a_1.shape[0] == a_2.shape[0] and \
       a_1.shape[1] == a_2.shape[1] and \
       a_1.shape[2] == a_2.shape[2]:
        print('All shapes of files are equal.')
        print('Calc diff of volume and surface.')
        two_netcdf_files_volume(filepath, a_1, a_2, dim0, dim1, dim2)
    elif a_1.shape[1] <= a_2.shape[1] and \
         a_1.shape[2] <= a_2.shape[2]:
        print('Not all shapes of files are equal.')
        print('Calc diff only of surface.')
        two_netcdf_files_surface(filepath, a_1, a_2)
    elif a_1.shape[1] >= a_2.shape[1] and \
         a_1.shape[2] >= a_2.shape[2]:
        print('Not all shapes of files are equal.')
        print('Calc diff only of surface.')
        two_netcdf_files_surface(filepath, a_1, a_2)
    else:
        print('One dim of surface of file 1 is bigger than same dim of',
              'surface of file 2.')
        print('The other dim of surface of file 1 is smaller than same dim of',
              'surface of file 2.')
        print('Cannot work with shape of files.')
        print('Aborting.')
        exit()

def two_dat_files(filepath):
    print('Read data from {0}.'.format(filepath[0]))
    a_1 = np.loadtxt(filepath[0])
    print('Read data from {0}.'.format(filepath[1]))
    a_2 = np.loadtxt(filepath[1])

    filepath_dat = os.path.splitext(filepath[0])[0] + '_' \
               + os.path.basename(os.path.splitext(filepath[1])[0]) \
               + '_diff'

    print('Write surface data to {0}.'.format(filepath_dat + '_abs.dat'))
    diff_file_abs = open(filepath_dat + '_abs.dat', 'w')
    print('Write surface data to {0}.'.format(filepath_dat + '_rel.dat'))
    diff_file_rel = open(filepath_dat + '_rel.dat', 'w')

    diff_max = 0.0
    diff_sum = 0.0

    if a_1.shape[0] != a_2.shape[0]:
        print('Number of lines do not match')
        print('Aborting.')
        exit()
    for num in range(0, a_1.shape[0]):
        if a_1.shape[1] != 3:
            print('Not three values in file {0} in line {1}.'.format(filepath_1,
                                                                     num))
            print('Aborting.')
            exit()
        if a_2.shape[1] != 3:
            print('Not three values in file {0} in line {1}.'.format(filepath_2,
                                                                     num))
            print('Aborting.')
            exit()
        if int(a_1[num,0]) != int(a_2[num,0]):
            print('X coordinates do not match in line {0}.')
            print('Aborting.')
            exit()
        if int(a_1[num,1]) != int(a_2[num,1]):
            print('Y coordinates do not match in line {0}.')
            print('Aborting.')
            exit()
        diff = a_1[num,2] - a_2[num,2]
        diff_sum += abs(diff)
        if abs(diff) > abs(diff_max):
            diff_max = diff
        string = str(int(a_1[num,0])) + ' ' + str(int(a_1[num,1])) + ' '
        string_abs = string + str(diff) + '\n'
        string_rel = string + str((diff/a_1[num,2])*100) + '\n'
        diff_file_abs.write(string_abs)
        diff_file_rel.write(string_rel)

    diff_mean = diff_sum/a_1.shape[0]
    print('Surface:  max(abs(diff)) = {0}.'.format(abs(diff_max)))
    print('Surface: mean(abs(diff)) = {0}.'.format(diff_mean))

    diff_file_abs.close()
    diff_file_rel.close()

def netcdf_and_dat_file(filepath):
    import netCDF4 as nc

    if os.path.splitext(filepath[0])[1] == '.nc':
        print('Read data from {0}.'.format(filepath[0]))
        nc_file = nc.Dataset(filepath[0])
        filepath[0] = os.path.splitext(filepath[0])[0]
        filepath[0] += '_surface.dat'
        print('Write surface data to {0}.'.format(filepath[0]))
        text_file = open(filepath[0], 'w')
    else:
        print('Read data from {0}.'.format(filepath[1]))
        nc_file = nc.Dataset(filepath[1])
        filepath[1] = os.path.splitext(filepath[1])[0]
        filepath[1] += '_surface.dat'
        print('Write surface data to {0}.'.format(filepath[1]))
        text_file = open(filepath[1], 'w')

    dim0 = nc_file.dimensions['nNodes_0'].size
    dim1 = nc_file.dimensions['nNodes_1'].size
    dim2 = nc_file.dimensions['nNodes_2'].size
    time = nc_file.dimensions['time'].size
    T = nc_file.variables['TNewDom']

    a = np.zeros((dim1, dim0))
    a[:,:] = T[(time-1):time,(dim2-1),:,:]

    nc_file.close()

    for elem_y in range(0, a.shape[0]):
        for elem_x in range(0, a.shape[1]):
            text_file.write('{0} {1} {2}\n'.format(str(elem_x), str(elem_y),
                                                   str(a[elem_y, elem_x])))
        text_file.write('\n')

    text_file.close()

    two_dat_files(filepath)

def main():
    filepath = []
    # Check if paths to two files (netCDF or .dat file, i.e. resuls) are provided,
    # if files exist and if files either have .nc or .dat extension.
    if len(sys.argv) > 2:
        for elem in range(1, 3):
            if os.path.isfile(sys.argv[elem]) == True:
                if os.path.splitext(sys.argv[elem])[1] == '.nc' or \
                   os.path.splitext(sys.argv[elem])[1] == '.dat':
                    filepath.append(sys.argv[elem])
                else:
                    print(sys.argv[elem],
                          'does not have .nc or .dat extension.')
            else:
                print(sys.argv[elem], 'does not exist.')
    else:
        print('Not enough command line arguments for netCDF files or .dat',
              'files provided.')

    if len(filepath) < 2:
        print('Usage: python3', sys.argv[0],
              '<PATH/TO/FILE_1> <PATH/TO/FILE_2>')
        print('Aborting.')
        exit()

    if os.path.splitext(filepath[0])[1] == '.nc' and \
       os.path.splitext(filepath[1])[1] == '.nc':
        two_netcdf_files(filepath)
    elif os.path.splitext(filepath[0])[1] == '.dat' and \
         os.path.splitext(filepath[1])[1] == '.dat':
        two_dat_files(filepath)
    else:
        netcdf_and_dat_file(filepath)

    print('Done.')

if __name__ == '__main__':
    main()
