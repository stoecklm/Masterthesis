import numpy as np
import os
import sys

def main():
    if len(sys.argv) > 2:
        if os.path.isfile(sys.argv[1]) == True and os.path.isfile(sys.argv[2]) == True:
            filepath_1 = sys.argv[1]
            filepath_2 = sys.argv[2]
        elif os.path.isfile(sys.argv[1]) == False and os.path.isfile(sys.argv[2]) == False:
            print(sys.argv[1], 'does not exist.')
            print(sys.argv[2], 'does not exist.')
            print('Usage: python3', sys.argv[0], '<PATH/TO/FILE_1> <PATH/TO/FILE_2>')
            print('Aborting.')
            exit()
        elif os.path.isfile(sys.argv[1]) == False:
            print(sys.argv[1], 'does not exist.')
            print('Usage: python3', sys.argv[0], '<PATH/TO/FILE_1> <PATH/TO/FILE_2>')
            print('Aborting.')
            exit()
        else: #os.path.isfile(sys.argv[2]) == False:
            print(sys.argv[2], 'does not exist.')
            print('Usage: python3', sys.argv[0], '<PATH/TO/FILE_1> <PATH/TO/FILE_2>')
            print('Aborting.')
            exit()
    elif len(sys.argv) == 2:
        print('Only one command line argument for files provided.')
        print('Usage: python3', sys.argv[0], '<PATH/TO/FILE_1> <PATH/TO/FILE_2>')
        print('Aborting.')
        exit()
    else:
        print('No command line arguments for files provided.')
        print('Usage: python3', sys.argv[0], '<PATH/TO/FILE_1> <PATH/TO/FILE_2>')
        print('Aborting.')
        exit()

    diff_file_abs = open('diff_abs.dat', 'w')
    diff_file_rel = open('diff_rel.dat', 'w')

    a_1 = np.loadtxt(filepath_1)
    a_2 = np.loadtxt(filepath_2)

    diff_max = 0.0
    diff_sum = 0.0

    if a_1.shape[0] != a_2.shape[0]:
        print('Number of lines do not match')
        print('Aborting.')
        exit()
    for num in range(0, a_1.shape[0]):
        if a_1.shape[1] != 3:
            print('Not three values in file {0} in line {1}.'.format(filepath_1, num))
            print('Aborting.')
            exit()
        if a_2.shape[1] != 3:
            print('Not three values in file {0} in line {1}.'.format(filepath_2, num))
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
    print('max(abs(diff)) = {0}.'.format(abs(diff_max)))
    print('mean(abs(diff)) = {0}.'.format(diff_mean))

if __name__ == '__main__':
    main()
