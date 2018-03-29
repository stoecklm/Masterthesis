import csv
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import linalg as LA
import os
import re
from scipy.interpolate import splprep, splev
import scipy.linalg
import sys
import warnings
warnings.filterwarnings(action="ignore", module="scipy",
                        message="^internal gelsd")
warnings.filterwarnings(action="ignore", module="scipy",
                        message="^Setting x")

def plot_points(points, filename):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    for point in points:
        ax.plot([point[0]], [point[1]], [point[2]],
                markeredgecolor='k', marker='o', markersize=5, alpha=0.6)

    ax.legend([i for i in range(1, points.shape[0]+1)])
    print('Save figure to {}.eps.'.format(filename))
    plt.savefig(filename + '.eps')
    plt.close()
    print('Done.')

def print_length(points):
    for point in points:
        print('Length of Points:')
        print((point[0]**2 + point[1]**2 + point[2]**2)**0.5)

def move_points(points, tumor, move):
    print('Move points.')
    for point in points:
        point[...] = point - move

    tumor = tumor - move

    print('Done.')
    return points, tumor

# https://stackoverflow.com/a/31466013
def interpolation(points):
    print('Doing interpolation.')
    # Add last point, since this value will be overwritten by splrep.
    pts = np.zeros((points.shape[0]+1,2))
    pts[0:points.shape[0],:] = points[:,0:2]
    pts[-1,:] = points[-1,0:2]
    # Interpolation
    tck, u = splprep(pts.T, u=None, s=0.0, per=1)
    u_new = np.linspace(u.min(), u.max(), 1000)
    x_new, y_new = splev(u_new, tck, der=0)

    a = np.zeros((1000,2))
    a[:,0] = x_new
    a[:,1] = y_new

    print('Done.')
    return a

# https://stackoverflow.com/a/31466013
def plot_interpolation(points, filename):
    print('Plot interpolation.')
    # Add last point, since this value will be overwritten by splrep.
    pts = np.zeros((points.shape[0]+1,2))
    pts[0:points.shape[0],:] = points[:,0:2]
    pts[-1,:] = points[-1,0:2]
    # Interpolation
    tck, u = splprep(pts.T, u=None, s=0.0, per=1)
    u_new = np.linspace(u.min(), u.max(), 1000)
    x_new, y_new = splev(u_new, tck, der=0)

    plt.plot(pts[:,0], pts[:,1], 'ro')
    plt.plot(x_new, y_new, 'b--')
    print('Save figure to {}.eps'.format(filename))
    plt.savefig(filename + '.eps')
    plt.close()

    print('Done.')

def read_intra_op_points(folderpath):
    print('Read IntraOp points.')
    if os.path.isfile(os.path.join(folderpath, 'fiducials.csv')):
        filepath = os.path.join(folderpath, 'fiducials.csv')
    else:
        filepath = os.path.join(folderpath, 'OpenIGTLink.fcsv')
    points = list()
    with open(filepath, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            # Check if this line is a IntraOp Point or Tumor.
            try:
                tmp = re.search('IntraOp', str(row[-3]))
                try:
                    # If it does not fail, it is a IntraOp Point
                    tmp.group(0)
                    x = float(row[1])
                    y = float(row[2])
                    z = float(row[3])
                    xyz = [x, y, z]
                    points.append(xyz)
                except AttributeError:
                    pass
            except IndexError:
                pass

    print('Done.')
    points = np.asarray(points)
    return points

def read_tumor_point(folderpath):
    print('Read tumor point.')
    if os.path.isfile(os.path.join(folderpath, 'fiducials.csv')):
        filepath = os.path.join(folderpath, 'fiducials.csv')
    else:
        filepath = os.path.join(folderpath, 'OpenIGTLink.fcsv')
    xyz = [0, 0, 0]
    with open(filepath, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            # Check if this line is a IntraOp Point or Tumor.
            try:
                tmp = re.search('Tumor', str(row[-3]))
                try:
                    # If it does not fail, it is the Tumor Center.
                    tmp.group(0)
                    x = float(row[1])
                    y = float(row[2])
                    z = float(row[3])
                    xyz = [x, y, z]
                except AttributeError:
                    pass
            except IndexError:
                pass

    print('Done.')
    tumor = np.asarray(xyz)
    return tumor

# https://gist.github.com/amroamroamro/1db8d69b4b65e8bc66a6
def plot_lin_plane_fitting(points, filename):
    print('Plot plane fitting.')
    data = np.c_[points[:,0], points[:,1], points[:,2]]
    # Regular grid covering the domain of the data.
    mn = np.min(data, axis=0)
    mx = np.max(data, axis=0)
    X,Y = np.meshgrid(np.linspace(mn[0], mx[0], 20),
                      np.linspace(mn[1], mx[1], 20))
    # Fit linear plane.
    A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
    C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])    # coefficients

    # Evaluate plane.
    Z = C[0]*X + C[1]*Y + C[2]

    # Plot points and fitted surface
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_wireframe(X, Y, Z, rstride=2, cstride=2)
    ax.scatter(data[:,0], data[:,1], data[:,2], c='r', s=50)
    plt.xlabel('X')
    plt.ylabel('Y')
    ax.set_zlabel('Z')
    ax.axis('equal')
    ax.axis('tight')
    print('Save figure to {}.eps'.format(filename))
    plt.savefig(filename + '.eps')
    plt.close()

    print('Done.')

# https://gist.github.com/amroamroamro/1db8d69b4b65e8bc66a6
def lin_plane_fitting(points):
    print('Plane fitting.')
    X = points[:,0]
    Y = points[:,1]
    Z = points[:,2]
    data = np.c_[X,Y,Z]
    # regular grid covering the domain of the data
    mn = np.min(data, axis=0)
    mx = np.max(data, axis=0)
    X,Y = np.meshgrid(np.linspace(mn[0], mx[0], 20),
                      np.linspace(mn[1], mx[1], 20))
    # best-fit linear plane
    A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
    C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])    # coefficients

    print('Done.')
    return C

# https://stackoverflow.com/a/9423864
def rotate_points(points, tumor):
    print('Rotate IntraOp points and tumor.')
    C = lin_plane_fitting(points)
    M = np.array([-1.0*C[0], -1.0*C[1], 1])
    N = np.array([0, 0, 1])

    costheta = M.dot(N) / (LA.norm(M)*LA.norm(N))

    axis = np.cross(M, N) / LA.norm(np.cross(M, N))

    c = costheta
    s = np.sqrt(1-c*c)
    C = 1-c
    x = axis[0]
    y = axis[1]
    z = axis[2]

    rmat = np.array([[x*x*C+c, x*y*C-z*s, x*z*C+y*s],
                     [y*x*C+z*s, y*y*C+c, y*z*C-x*s],
                     [z*x*C-y*s, z*y*C+x*s, z*z*C+c]])

    for point in points:
        point[...] = np.dot(rmat, point)

    tumor = np.dot(rmat, tumor)

    print('Done.')
    return points, tumor

def get_interpolated_path(points):
    return mpath.Path(interpolation(points))

def get_path(points):
    return mpath.Path(points[:,0:2])

def main():
    filepath = ''
    # Check if path to csv file (i.e. results) is provided,
    # if file exists and if file has .csv extension.
    if len(sys.argv) > 1:
        if os.path.isdir(sys.argv[1]) == True:
            tmp1 = os.path.join(os.sys.argv[1], 'fiducials.csv')
            tmp2 = os.path.join(os.sys.argv[1], 'OpenIGTLink.fcsv')
            if os.path.isfile(tmp1) != True and os.path.isfile(tmp2) != True:
                print(sys.argv[1], 'does not contain fiducials.csv',
                      'or OpenIGTLink.fcsv.')
                print('Aborting.')
                exit()
            else:
                filepath = sys.argv[1]
                filepath = os.path.normpath(filepath)
        else:
            print(sys.argv[1], 'does not exist.')
            print('Aborting.')
            exit()
    else:
        print('No command line argument for folder provided.')
        print('Running test.')
        filepath = 'test'

    # Original data.
    iop = read_intra_op_points(filepath)
    print('Set of IntraOp points:')
    print(iop)
    t = read_tumor_point(filepath)
    print('Tumor point:')
    print(t)
    plot_points(iop, filepath + '_org_points')
    plot_lin_plane_fitting(iop, filepath + '_org_plane')
    plot_interpolation(iop, filepath + '_org_inter')
    print()
    # Rotation of points.
    iop, t = rotate_points(iop, t)
    plot_points(iop, filepath + '_rot_points')
    plot_lin_plane_fitting(iop, filepath + '_rot_plane')
    plot_interpolation(iop, filepath + '_rot_inter')
    print('Set of IntraOp points after rotation:')
    print(iop)
    print('Tumor point after rotation:')
    print(t)
    print()
    # Move tumor to origin.
    iop, t = move_points(iop, t, t)
    plot_points(iop, filepath + '_move_points')
    plot_lin_plane_fitting(iop, filepath + '_move_plane')
    plot_interpolation(iop, filepath + '_move_inter')
    print('Set of IntraOp points after moving:')
    print(iop)
    print('Tumor point after moving:')
    print(t)
    print()

if __name__ == '__main__':
    main()
