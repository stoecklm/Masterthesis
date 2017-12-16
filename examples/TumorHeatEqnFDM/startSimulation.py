import configparser
import math
import os
import subprocess
import sys

def main():
    if len(sys.argv) > 1:
        if os.path.isfile(sys.argv[1]) == True:
            NAME_CONFIGFILE = sys.argv[1]
        else:
            print(sys.argv[1], "does not exist.")
            print("Usage: python3", sys.argv[0], "<PATH/TO/CONFIGFILE>")
            print("Aborting.")
            exit()
    else:
        print("No command line argument for configfile provided.")
        print("Usage: python3", sys.argv[0], "<PATH/TO/CONFIGFILE>")
        print("Aborting.")
        exit()

    config = configparser.ConfigParser()
    config.read(NAME_CONFIGFILE)

    SPACE_DIM = config['Dimension'].getint('SPACE_DIM')

    GRIDSIZE_GLOBAL = config['Geometry'].getfloat('GRIDSIZE_GLOBAL')
    COORD_NODE_FIRST_DIM1 = config['Geometry'].getfloat('COORD_NODE_FIRST_DIM1')
    COORD_NODE_FIRST_DIM2 = config['Geometry'].getfloat('COORD_NODE_FIRST_DIM2')
    COORD_NODE_FIRST_DIM3 = config['Geometry'].getfloat('COORD_NODE_FIRST_DIM3')
    COORD_NODE_LAST_DIM1 = config['Geometry'].getfloat('COORD_NODE_LAST_DIM1')
    COORD_NODE_LAST_DIM2 = config['Geometry'].getfloat('COORD_NODE_LAST_DIM2')
    COORD_NODE_LAST_DIM3 = config['Geometry'].getfloat('COORD_NODE_LAST_DIM3')

    N_NODES_DIM1 = (COORD_NODE_LAST_DIM1 - COORD_NODE_FIRST_DIM1)/GRIDSIZE_GLOBAL
    N_NODES_DIM1 = int(math.ceil(N_NODES_DIM1))
    N_NODES_DIM2 = (COORD_NODE_LAST_DIM2 - COORD_NODE_FIRST_DIM2)/GRIDSIZE_GLOBAL
    N_NODES_DIM2 = int(math.ceil(N_NODES_DIM2))
    N_NODES_DIM3 = (COORD_NODE_LAST_DIM3 - COORD_NODE_FIRST_DIM3)/GRIDSIZE_GLOBAL
    N_NODES_DIM3 = int(math.ceil(N_NODES_DIM3))

    START_TIME = config['Time'].getint('START_TIME')
    END_TIME = config['Time'].getint('END_TIME')
    DELTA_TIME = config['Time'].getfloat('DELTA_TIME')

    N_TIMESTEPS = (END_TIME - START_TIME)/DELTA_TIME
    N_TIMESTEPS = int(math.ceil(N_TIMESTEPS))

    N_SNAPSHOTS = config['Output'].getint('N_SNAPSHOTS')
    if N_SNAPSHOTS > N_TIMESTEPS:
        N_SNAPSHOTS = N_TIMESTEPS

    os.putenv("SCAFESRUN_N_TIMESTEPS", str(N_TIMESTEPS))
    os.putenv("SCAFESRUN_N_SNAPSHOTS", str(N_SNAPSHOTS))
    os.putenv("SCAFESRUN_START_TIME", str(START_TIME))
    os.putenv("SCAFESRUN_END_TIME", str(END_TIME))

    if SPACE_DIM == 2:
        tmp = str(COORD_NODE_FIRST_DIM1) + "x" + str(COORD_NODE_FIRST_DIM2)
        os.putenv("SCAFESRUN_COORD_NODE_FIRST", tmp)
        tmp = str(COORD_NODE_LAST_DIM1) + "x" + str(COORD_NODE_LAST_DIM2)
        os.putenv("SCAFESRUN_COORD_NODE_LAST", tmp)
        tmp = str(N_NODES_DIM1) + "x" + str(N_NODES_DIM2)
        os.putenv("SCAFESRUN_N_NODES", tmp)
    if SPACE_DIM == 3:
        tmp = str(COORD_NODE_FIRST_DIM1) + "x" + str(COORD_NODE_FIRST_DIM2) + "x" + str(COORD_NODE_FIRST_DIM3)
        os.putenv("SCAFESRUN_COORD_NODE_FIRST", tmp)
        tmp = str(COORD_NODE_LAST_DIM1) + "x" + str(COORD_NODE_LAST_DIM2) + "x" + str(COORD_NODE_LAST_DIM3)
        os.putenv("SCAFESRUN_COORD_NODE_LAST", tmp)
        tmp = str(N_NODES_DIM1) + "x" + str(N_NODES_DIM2) + "x" + str(N_NODES_DIM3)
        os.putenv("SCAFESRUN_N_NODES", tmp)
    os.putenv("SCAFESRUN_NAME_CONFIGFILE", NAME_CONFIGFILE)
    os.putenv("SCAFESRUN_SPACE_DIM", str(SPACE_DIM))

    subprocess.call("./RUN_HELPER.sh")

if __name__ == '__main__':
    main()
