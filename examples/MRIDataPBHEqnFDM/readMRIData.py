import csv
import numpy as np
import os
import sys
from scipy.interpolate import splprep, splev
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_points(points, filepath):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    for point in points:
        ax.plot([point[0]], [point[1]], [point[2]],
                markeredgecolor='k', marker='o', markersize=5, alpha=0.6)

    ax.legend([i for i in range(1, points.shape[0]+1)])
    print('Save figure to {}.'.format(filepath))
    plt.savefig(filepath)
    plt.close()

def get_R_x(angle):
    theta_x = np.radians(angle)

    R_x = np.array([[1, 0, 0],
                    [0, np.cos(theta_x), -1.0*np.sin(theta_x)],
                    [0, np.sin(theta_x), np.cos(theta_x)]])

    return R_x

def get_R_y(angle):
    theta_y = np.radians(angle)

    R_y = np.array([[np.cos(theta_y), 0, np.sin(theta_y)],
                    [0, 1, 0],
                    [-1.0*np.sin(theta_y), 0, np.cos(theta_y)]])

    return R_y

def get_R_z(angle):
    theta_z = np.radians(angle)

    R_z = np.array([[np.cos(theta_z), -1.0*np.sin(theta_z), 0],
                    [np.sin(theta_z), np.cos(theta_z), 0],
                    [0, 0, 1]])

    return R_z

def print_length(points):
    for point in points:
        print((point[0]**2 + point[1]**2 + point[2]**2)**0.5)

def move_points(points, move):
    for point in points:
        point[...] = point - move

    return points

def interpolation(pp, filepath):
    # Add last point, since this value will be overwritten by splrep.
    pts = np.zeros((pp.shape[0]+1,2))
    pts[0:pp.shape[0],:] = pp[:,0:2]
    pts[-1,:] = pp[-1,0:2]
    # Interpolation
    tck, u = splprep(pts.T, u=None, s=0.0, per=1)
    u_new = np.linspace(u.min(), u.max(), 1000)
    x_new, y_new = splev(u_new, tck, der=0)

    plt.plot(pts[:,0], pts[:,1], 'ro')
    plt.plot(x_new, y_new, 'b--')
    print('Save figure to {}.'.format(filepath))
    plt.savefig(filepath)
    plt.close()

def readMRIData(filepath):
    case = os.path.dirname(filepath)
    points = list()
    with open(filepath, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            x = float(row[1])
            y = float(row[2])
            z = float(row[3])
            xyz = [x, y, z]
            points.append(xyz)

    p = np.array([np.asarray(points[0]),
                  np.asarray(points[1]),
                  np.asarray(points[2]),
                  np.asarray(points[3]),
                  np.asarray(points[4]),
                  np.asarray(points[5]),
                  np.asarray(points[6]),
                  np.asarray(points[7])])

    t = np.asarray(points[8])

    print('Original set of Points:')
    print(p)
    print('Original tumor location:')
    print(t)
    print()
    plot_points(p, case + '_org.eps')
    interpolation(p, case + '_org_interpolate.eps')

    deg_x_min = 0
    deg_y_min = 0
    z_diff_min = 1000000000
    for deg_x in range(-180, 0):
        for deg_y in range(-180, 0):
            R_x = get_R_x(deg_x)
            R_y = get_R_y(deg_y)
            R_z = get_R_z(0)

            pp = p

            pp = pp.dot(R_x)
            pp = pp.dot(R_y)
            pp = pp.dot(R_z)

            z_min = np.amin(pp[:,2])
            z_max = np.amax(pp[:,2])
            z_diff = abs(z_max - z_min)
            if abs(z_diff) < abs(z_diff_min):
                z_diff_min = z_diff
                deg_x_min = deg_x
                deg_y_min = deg_y

    print()
    print('x angle:', deg_x_min)
    print('y angle:', deg_y_min)
    print('Diff between z_min and z_max:', z_diff_min)
    print()

    pp = p
    pp = pp.dot(get_R_x(deg_x_min))
    pp = pp.dot(get_R_y(deg_y_min))
    pp = pp.dot(get_R_z(0))
    print('Set of Points ater rotation:')
    print(pp)
    t = t.dot(get_R_x(deg_x_min))
    t = t.dot(get_R_y(deg_y_min))
    t = t.dot(get_R_z(0))
    print('Tumor location after rotation:')
    print(t)
    print()

    plot_points(pp, case + '_rotated_points.eps')
    interpolation(pp, case + '_rotated_interpolate.eps')
    print()

    pp = move_points(pp, t)
    print('Set of Points ater move:')
    print(pp)
    t = t - t
    print('Tumor location after move:')
    print(t)
    print()

    plot_points(pp, case + '_moved_points.eps')
    interpolation(pp, case + '_moved_interpolate.eps')
    print()

    print('Done.')


def main():
    filepath = ''
    # Check if path to csv file (i.e. results) is provided,
    # if file exists and if file has .csv extension.
    if len(sys.argv) > 1:
        if os.path.isfile(sys.argv[1]) == True:
            if os.path.splitext(sys.argv[1])[1] == '.csv':
                filepath = sys.argv[1]
            else:
                print(sys.argv[1], 'does not have .csv extension.')
        else:
            print(sys.argv[1], 'does not exist.')
    else:
        print('No command line argument for csv file provided.')

    if filepath == '':
        print('Usage: python3', sys.argv[0], '<PATH/TO/FILE>')
        print('Aborting.')
        exit()

    readMRIData(filepath)


if __name__ == '__main__':
    main()
