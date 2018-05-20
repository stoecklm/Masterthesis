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

def save_1d_mcmc_fit_results_as_netcdf(l2_norm, variable, filename):
    print('Save data to {}.'.format(filename))
    nc_file = nc.Dataset(filename, 'w', format='NETCDF3_CLASSIC')
    nc_file.createDimension('iterations', l2_norm.shape[0])
    values = nc_file.createVariable('L2_Norm', 'f8', ('iterations'))
    values[:] = l2_norm[:]
    values = nc_file.createVariable('tested_variable', 'f8', ('iterations'))
    values[:] = variable[:]
    nc_file.close()

    print('Done.')

if __name__ == '__main__':
    pass
