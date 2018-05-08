import configparser
import os
from shutil import copyfile
import sys

import netCDF4 as nc
import numpy as np
from numpy import linalg as LA
import nrrd
from scipy.interpolate import RegularGridInterpolator

from readMRIData import get_rotation_matrix
from readMRIData import plot_lin_plane_fitting_with_bbox
from readMRIData import read_intra_op_points
from readMRIData import rotate_point_by_rotation_matrix

def read_nrrd_file(filepath):
    data, header = nrrd.read(filepath)

    return data, header

def save_as_netcdf(data, filename):
    print('Save data to {}.'.format(filename))
    nc_file = nc.Dataset(filename, 'w', format='NETCDF3_CLASSIC')
    nc_file.createDimension('nNodes_0', data.shape[0])
    nc_file.createDimension('nNodes_1', data.shape[1])
    nc_file.createDimension('nNodes_2', data.shape[2])
    nc_file.createDimension('time')
    brain = nc_file.createVariable('region', 'i2', ('time', 'nNodes_2',
                                                    'nNodes_1', 'nNodes_0'))
    brain[0,] = np.swapaxes(data, 0, 2)
    nc_file.close()

    print('Done.')

def data_as_binary_data(data):
    binary_data = np.zeros(data.shape)
    binary_data = np.where(data > 1.0, 1, 0)

    return binary_data

def get_ijk_to_lps(header):
    space_origin = return_string_list_as_float_numpy_array(header['space origin'])
    space_origin = np.append(space_origin, 1)

    row0 = return_string_list_as_float_numpy_array(header['space directions'][0])
    row0 = np.append(row0, 0)
    row1 = return_string_list_as_float_numpy_array(header['space directions'][1])
    row1 = np.append(row1, 0)
    row2 = return_string_list_as_float_numpy_array(header['space directions'][2])
    row2 = np.append(row2, 0)
    ijk_to_lps = np.array([row0, row1, row2, space_origin]).T

    return ijk_to_lps

def ijk_to_lps(ijk, header):
    if len(ijk) == 3:
        ijk = np.append(ijk, 1)

    ijk_to_lps = get_ijk_to_lps(header)

    lps = np.dot(ijk_to_lps, ijk)

    return lps[0:3]

def ijk_to_ras(ijk, header):
    if len(ijk) == 3:
        ijk = np.append(ijk, 1)

    ijk_to_lps = get_ijk_to_lps(header)
    lps_to_ras = np.diag([-1, -1, 1, 1])

    lps = np.dot(ijk_to_lps, ijk)
    ras = np.matmul(lps, lps_to_ras)

    return ras[0:3]

def lps_to_ijk(lps, header):
    if len(lps) == 3:
        lps = np.append(lps, 1)

    ijk_to_lps = get_ijk_to_lps(header)
    lps_to_ijk = LA.inv(ijk_to_lps)

    ijk = np.dot(lps_to_ijk, lps)

    return ijk[0:3]

def ras_to_ijk(ras, header):
    if len(ras) == 3:
        ras = np.append(ras, 1)

    ijk_to_lps = get_ijk_to_lps(header)
    lps_to_ijk = LA.inv(ijk_to_lps)
    ras_to_lps = np.diag([-1, -1, 1, 1])

    lps = np.matmul(ras, ras_to_lps)
    ijk = np.dot(lps_to_ijk, lps)

    return ijk[0:3]

def switch_space(point):
    if len(point) == 3:
        point = np.append(point, 1)

    switch_mat = np.diag([-1, -1, 1, 1])

    switched_point = np.matmul(point, switch_mat)

    return switched_point[0:3]

def get_regular_grid_interpolator(data):
    x = np.linspace(0, data.shape[0]-1, data.shape[0])
    y = np.linspace(0, data.shape[1]-1, data.shape[1])
    z = np.linspace(0, data.shape[2]-1, data.shape[2])

    reg_grid_interpolator = RegularGridInterpolator((x, y, z), data,
                                                    method='linear',
                                                    bounds_error=False,
                                                    fill_value=0)

    return reg_grid_interpolator

def read_default_ini_file(params):

    config = configparser.ConfigParser()
    config.optionxform = str

    if os.path.isfile('Parameters.ini') == True:
        config.read('Parameters.ini')
    else:
        print('Parameters.ini does not exist.')
        print('Aborting.')
        exit()

    COORD_NODE_FIRST = config['Geometry'].get('COORD_NODE_FIRST')
    COORD_NODE_FIRST = list(map(float, COORD_NODE_FIRST.split('x')))
    params['COORD_NODE_FIRST'] = COORD_NODE_FIRST

    COORD_NODE_LAST = config['Geometry'].get('COORD_NODE_LAST')
    COORD_NODE_LAST = list(map(float, COORD_NODE_LAST.split('x')))
    params['COORD_NODE_LAST'] = COORD_NODE_LAST

    N_NODES = config['Geometry'].get('N_NODES')
    N_NODES = list(map(int, N_NODES.split('x')))
    params['N_NODES'] = N_NODES

    return params

def write_ini_file(params, filename):
    case = params['CASE']
    print('Write {}.ini.'.format(case))

    folderpath = params['FOLDERPATH']
    start = params['START']
    end = params['END']

    start_as_string = str(start[0]/1000) + 'x' + str(start[1]/1000) + 'x' \
                      + str(start[2]/1000)
    end_as_string = str(end[0]/1000) + 'x' + str(end[1]/1000) + 'x' \
                    + str(end[2]/1000)

    config = configparser.ConfigParser()
    config.optionxform = str

    if os.path.isfile('Parameters.ini') == True:
        config.read('Parameters.ini')
    else:
        print('Parameters.ini does not exist.')
        print('Aborting.')
        exit()

    config['MRI']['CASE'] = case
    config['Input']['NAME_REGION_FILE'] = os.path.join(folderpath, 'region')
    config['Input']['NAME_VESSELS_FILE'] = os.path.join(folderpath, 'vessels')
    config['Input']['NAME_INITFILE'] = os.path.join(folderpath, 'init')
    config['Input']['USE_MRI_FILE'] = 'True'
    config['Input']['USE_INITFILE'] = 'True'
    config['Input']['CREATE_INITFILE'] = 'True'
    config['Geometry']['COORD_NODE_FIRST'] = start_as_string
    config['Geometry']['COORD_NODE_LAST'] = end_as_string

    with open(case + '.ini', 'w') as configfile:
        config.write(configfile)

    print('Done.')

def return_string_list_as_float_numpy_array(list_to_array):
    numpy_array = np.asarray(list(map(float, list_to_array)))

    return numpy_array

def return_float_numpy_array_as_string_list(array_to_list):
    string_list = list(map(str, array_to_list))

    return string_list

def rotate_tumor_data(params, data, iop, header_rot):
    dim0_length = (params['COORD_NODE_LAST'][0] \
                   - params['COORD_NODE_FIRST'][0])*1000
    dim1_length = (params['COORD_NODE_LAST'][1] \
                   - params['COORD_NODE_FIRST'][1])*1000
    dim2_length = (params['COORD_NODE_LAST'][2] \
                   - params['COORD_NODE_FIRST'][2])*1000

    x_min = np.min(iop[:,0])
    x_max = np.max(iop[:,0])
    y_min = np.min(iop[:,1])
    y_max = np.max(iop[:,1])

    #z = np.polyfit(iop[:,0], iop[:,2], 0)
    z = np.max(iop[:,2])
    #z = np.min(iop[:,2])
    z = np.ceil(z)

    x_size = x_max - x_min
    y_size = y_max - y_min

    missing_x = dim0_length - x_size
    missing_y = dim1_length - y_size

    x_start = x_min - missing_x/2.0
    x_start = np.floor(x_start)
    x_end = dim0_length + x_start

    y_start = y_min - missing_y/2.0
    y_start = np.floor(y_start)
    y_end = y_start + dim1_length

    z_start = float(z) - dim2_length
    z_end = float(z)

    start = [x_start, y_start, z_start]
    end = [x_end, y_end, z_end]

    params['START'] = start

    params['END'] = end

    dim0 = params['N_NODES'][0]
    dim1 = params['N_NODES'][1]
    dim2 = params['N_NODES'][2]

    delta_x = dim0_length/(dim0-1)
    delta_y = dim1_length/(dim1-1)
    delta_z = dim2_length/(dim2-1)

    print('Create LPS coordinates array.')

    coord_lps = np.zeros(dim0*dim1*dim2*3).reshape((dim0, dim1, dim2, 3))

    for elem_x in range(0, coord_lps.shape[0]):
        x = x_start + elem_x * delta_x
        for elem_y in range(0, coord_lps.shape[1]):
            y = y_start + elem_y * delta_y
            for elem_z in range(0, coord_lps.shape[2]):
                z = z_start + elem_z * delta_z
                coord_lps[elem_x, elem_y, elem_z] = np.asarray([x, y, z])

    print('Done.')
    print('Create IJK coordinates array.')

    coord_ijk = np.zeros(dim0*dim1*dim2*3).reshape((dim0, dim1, dim2, 3))

    for elem_x in range(0, coord_ijk.shape[0]):
        for elem_y in range(0, coord_ijk.shape[1]):
            for elem_z in range(0, coord_ijk.shape[2]):
                coord_ijk[elem_x, elem_y, elem_z] = lps_to_ijk(coord_lps[elem_x,
                                                                         elem_y,
                                                                         elem_z],
                                                               header_rot)

    print('Done.')
    print('Interpolate data.')

    reg_grid_inter = get_regular_grid_interpolator(data)

    final_data = np.zeros(dim0*dim1*dim2).reshape((dim0, dim1, dim2))

    for elem_x in range(0, final_data.shape[0]):
        for elem_y in range(0, final_data.shape[1]):
            for elem_z in range(0, final_data.shape[2]):
                tmp = coord_ijk[elem_x, elem_y, elem_z]
                data_pt = reg_grid_inter(tmp)
                final_data[elem_x, elem_y, elem_z] = data_pt

    print('Done.')

    return final_data

def build_rotated_header(header, rmat):
    space_dir_0 = return_string_list_as_float_numpy_array(header['space directions'][0])
    space_dir_1 = return_string_list_as_float_numpy_array(header['space directions'][1])
    space_dir_2 = return_string_list_as_float_numpy_array(header['space directions'][2])
    space_origin = return_string_list_as_float_numpy_array(header['space origin'])

    space_dir_0_rot = rotate_point_by_rotation_matrix(space_dir_0, rmat)
    space_dir_1_rot = rotate_point_by_rotation_matrix(space_dir_1, rmat)
    space_dir_2_rot = rotate_point_by_rotation_matrix(space_dir_2, rmat)
    space_origin_rot = rotate_point_by_rotation_matrix(space_origin, rmat)

    header_rot = {'keyvaluepairs': ''}
    header_rot['space directions'] = [return_float_numpy_array_as_string_list(space_dir_0_rot),
                                      return_float_numpy_array_as_string_list(space_dir_1_rot),
                                      return_float_numpy_array_as_string_list(space_dir_2_rot)]
    header_rot['space origin'] = return_float_numpy_array_as_string_list(space_origin_rot)

    return header_rot

def ijk_bounding_box(start, end):
    bbox = np.asarray([[start[0], start[1], start[2]],
                       [end[0], start[1], start[2]],
                       [end[0], end[1], start[2]],
                       [start[0], end[1], start[2]],
                       [start[0], start[1], end[2]],
                       [end[0], start[1], end[2]],
                       [end[0], end[1], end[2]],
                       [start[0], end[1], end[2]]])

    return bbox

def ijk_bounding_box_to_lps(bbox, header):
    for point in bbox:
        point[...] = ijk_to_lps(point, header)

    return bbox


def main():
    params = {'keyvaluepairs' : ''}
    if len(sys.argv) > 1:
        if os.path.isdir(sys.argv[1]) == True:
            tmp = os.path.join(os.sys.argv[1], '_1.nrrd')
            if os.path.isfile(tmp) != True:
                print(sys.argv[1], 'does not contain _1.nrrd.')
                print('Aborting.')
                exit()
            else:
                folderpath = sys.argv[1]
                filepath = tmp
        else:
            print(sys.argv[1], 'does not exist.')
            print('Aborting.')
            exit()
    else:
        print('No command line argument for folder provided.')
        print('Aborting.')
        exit()

    print('Create region file.')

    case = filepath.split('_')[0]
    params['CASE'] = case
    params['FOLDERPATH'] = folderpath

    data, header = read_nrrd_file(filepath)
    save_as_netcdf(data, os.path.join(folderpath, 'tumor.nc'))
    read_default_ini_file(params)

    iop = read_intra_op_points(folderpath)

    for point in iop:
        point[...] = switch_space(point)

    rmat = get_rotation_matrix(iop)

    start_ijk = np.asarray([0, 0, 0])
    end_ijk = np.asarray(header['sizes']) - np.asarray([1.0, 1.0, 1.0])
    bbox_ijk = ijk_bounding_box(start_ijk, end_ijk)
    bbox_lps = ijk_bounding_box_to_lps(bbox_ijk, header)
    plot_lin_plane_fitting_with_bbox(iop, bbox_lps,
                                     os.path.join(folderpath, 'iop_and_bbox'))

    for point in iop:
        point[...] = rotate_point_by_rotation_matrix(point, rmat)

    header_rot = build_rotated_header(header, rmat)

    final_data = rotate_tumor_data(params, data, iop, header_rot)
    final_data = data_as_binary_data(final_data)
    filename = os.path.join(folderpath, 'region.nc')
    save_as_netcdf(final_data, filename)
    write_ini_file(params, filename)

    print('Done.')

if __name__ == '__main__':
    main()
