import netCDF4 as nc
import numpy as np

def temperature_array_from_result(filepath):
    print('Read {}.'.format(filepath))
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
        print('* ERROR: No temperature variable found in this file.')
        print('Aborting.')
        exit()

    temp = np.zeros((dim2, dim1, dim0))
    temp[:,:,:] = T[-1,:,:,:]

    nc_file.close()

    return temp

def surface_temperature_array_from_result(filepath):
    temp = temperature_array_from_result(filepath)
    dim2, dim1, dim0 = temp.shape

    temp_surface = np.zeros((dim1, dim0))
    temp_surface[:,:] = temp[-1,:,:]

    return temp_surface

if __name__ == '__main__':
    pass
