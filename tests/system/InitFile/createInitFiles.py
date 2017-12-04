import netCDF4 as nc
import numpy as np

dim = 8
file_name = "init"

def create_nc_file_1D():
    global dim
    global file_name
    nc_file_name = file_name + "1D.nc"
    nc_file = nc.Dataset(nc_file_name, "w", format="NETCDF3_CLASSIC")
    nNodes_0 = nc_file.createDimension("nNodes_0", dim)
    time = nc_file.createDimension("time")
    init_values = nc_file.createVariable("U", "f8", ("time", "nNodes_0"))
    num_elem = dim
    a = np.ones(num_elem).reshape(dim)
    init_values[0,:] = a
    init_values[0,0] = 100.0
    init_values[0,(dim-1)] = 100.0
    nc_file.close()

def create_nc_file_2D():
    global dim
    global file_name
    nc_file_name = file_name + "2D.nc"
    nc_file = nc.Dataset(nc_file_name, "w", format="NETCDF3_CLASSIC")
    nNodes_0 = nc_file.createDimension("nNodes_0", dim)
    nNodes_1 = nc_file.createDimension("nNodes_1", dim)
    time = nc_file.createDimension("time")
    init_values = nc_file.createVariable("U", "f8", ("time", "nNodes_0", "nNodes_1"))
    num_elem = dim * dim
    a = np.ones(num_elem).reshape(dim, dim)
    init_values[0,:,:] = a
    init_values[0,0,:] = 100.0
    init_values[0,:,0] = 100.0
    init_values[0,(dim-1),:] = 100.0
    init_values[0,:,(dim-1)] = 100.0
    nc_file.close()

def create_nc_file_3D():
    global dim
    global file_name
    nc_file_name = file_name + "3D.nc"
    nc_file = nc.Dataset(nc_file_name, "w", format="NETCDF3_CLASSIC")
    nNodes_0 = nc_file.createDimension("nNodes_0", dim)
    nNodes_1 = nc_file.createDimension("nNodes_1", dim)
    nNodes_2 = nc_file.createDimension("nNodes_2", dim)
    time = nc_file.createDimension("time")
    init_values = nc_file.createVariable("U", "f8", ("time", "nNodes_0", "nNodes_1", "nNodes_2"))
    num_elem = dim * dim * dim
    a = np.ones(num_elem).reshape(dim, dim, dim)
    init_values[0,:,:,:] = a
    init_values[0,0,:,:] = 100.0
    init_values[0,:,0,:] = 100.0
    init_values[0,:,:,0] = 100.0
    init_values[0,(dim-1),:,:] = 100.0
    init_values[0,:,(dim-1),:] = 100.0
    init_values[0,:,:,(dim-1)] = 100.0
    nc_file.close()

def main():
    create_nc_file_1D()
    create_nc_file_2D()
    create_nc_file_3D()

if __name__ == '__main__':
    main()
