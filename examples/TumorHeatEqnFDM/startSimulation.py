import configparser
import math
import netCDF4 as nc
import numpy as np
import os
import subprocess
import sys

params = {
    'NAME_CONFIGFILE' : '',
    'SPACE_DIM' : 0,
    'GRIDSIZE_GLOBAL' : 0.0,
    'COORD_NODE_FIRST_DIM1' : 0.0,
    'COORD_NODE_FIRST_DIM2' : 0.0,
    'COORD_NODE_FIRST_DIM3' : 0.0,
    'COORD_NODE_LAST_DIM1' : 0.0,
    'COORD_NODE_LAST_DIM2' : 0.0,
    'COORD_NODE_LAST_DIM3' : 0.0,
    'START_TIME' : 0,
    'END_TIME' : 0,
    'DELTA_TIME' : 0.0,
    'N_SNAPSHOTS': 0,
    'N_TIMESTEPS' : 0,
    'N_NODES_DIM1' : 0,
    'N_NODES_DIM2' : 0,
    'N_NODES_DIM3' : 0,
}

def parseConfigFile():
    global params

    config = configparser.ConfigParser()
    config.read(params['NAME_CONFIGFILE'])

    params['SPACE_DIM'] = config['Dimension'].getint('SPACE_DIM')

    params['GRIDSIZE_GLOBAL'] = config['Geometry'].getfloat('GRIDSIZE_GLOBAL')
    params['COORD_NODE_FIRST_DIM1'] = config['Geometry'].getfloat('COORD_NODE_FIRST_DIM1')
    params['COORD_NODE_FIRST_DIM2'] = config['Geometry'].getfloat('COORD_NODE_FIRST_DIM2')
    params['COORD_NODE_FIRST_DIM3'] = config['Geometry'].getfloat('COORD_NODE_FIRST_DIM3')
    params['COORD_NODE_LAST_DIM1'] = config['Geometry'].getfloat('COORD_NODE_LAST_DIM1')
    params['COORD_NODE_LAST_DIM2'] = config['Geometry'].getfloat('COORD_NODE_LAST_DIM2')
    params['COORD_NODE_LAST_DIM3'] = config['Geometry'].getfloat('COORD_NODE_LAST_DIM3')

    params['START_TIME'] = config['Time'].getint('START_TIME')
    params['END_TIME']= config['Time'].getint('END_TIME')
    params['DELTA_TIME'] = config['Time'].getfloat('DELTA_TIME')

    params['N_SNAPSHOTS'] = config['Output'].getint('N_SNAPSHOTS')

def calcVariables():
    global params

    N_NODES_DIM1 = (params['COORD_NODE_LAST_DIM1'] - params['COORD_NODE_FIRST_DIM1'])/params['GRIDSIZE_GLOBAL']
    params['N_NODES_DIM1'] = int(math.ceil(N_NODES_DIM1))
    N_NODES_DIM2 = (params['COORD_NODE_LAST_DIM2'] - params['COORD_NODE_FIRST_DIM2'])/params['GRIDSIZE_GLOBAL']
    params['N_NODES_DIM2'] = int(math.ceil(N_NODES_DIM2))
    N_NODES_DIM3 = (params['COORD_NODE_LAST_DIM3'] - params['COORD_NODE_FIRST_DIM3'])/params['GRIDSIZE_GLOBAL']
    params['N_NODES_DIM3'] = int(math.ceil(N_NODES_DIM3))

    N_TIMESTEPS = (params['END_TIME'] - params['START_TIME'])/params['DELTA_TIME']
    params['N_TIMESTEPS'] = int(math.ceil(N_TIMESTEPS))

    if params['N_SNAPSHOTS'] > params['N_TIMESTEPS']:
        params['N_SNAPSHOTS'] = params['N_TIMESTEPS']

def setEnvironmentVariables():
    global params

    os.putenv("SCAFESRUN_N_TIMESTEPS", str(params['N_TIMESTEPS']))
    os.putenv("SCAFESRUN_N_SNAPSHOTS", str(params['N_SNAPSHOTS']))
    os.putenv("SCAFESRUN_START_TIME", str(params['START_TIME']))
    os.putenv("SCAFESRUN_END_TIME", str(params['END_TIME']))

    os.putenv("SCAFESRUN_NAME_CONFIGFILE", params['NAME_CONFIGFILE'])
    os.putenv("SCAFESRUN_SPACE_DIM", str(params['SPACE_DIM']))

    COORD_NODE_FIRST = str(params['COORD_NODE_FIRST_DIM1'])
    COORD_NODE_LAST = str(params['COORD_NODE_LAST_DIM1'])
    N_NODES = str(params['N_NODES_DIM1'])

    if params['SPACE_DIM'] > 1:
        COORD_NODE_FIRST += "x" + str(params['COORD_NODE_FIRST_DIM2'])
        COORD_NODE_LAST += "x" + str(params['COORD_NODE_LAST_DIM2'])
        N_NODES += "x" + str(params['N_NODES_DIM2'])
    if params['SPACE_DIM'] > 2:
        COORD_NODE_FIRST += "x" + str(params['COORD_NODE_FIRST_DIM3'])
        COORD_NODE_LAST += "x" + str(params['COORD_NODE_LAST_DIM3'])
        N_NODES += "x" + str(params['N_NODES_DIM3'])

    os.putenv("SCAFESRUN_COORD_NODE_FIRST", COORD_NODE_FIRST)
    os.putenv("SCAFESRUN_COORD_NODE_LAST", COORD_NODE_LAST)
    os.putenv("SCAFESRUN_N_NODES", N_NODES)

    NAME_EXECUTABLE = os.path.basename(os.getcwd()) + str(params['SPACE_DIM']) + 'D'

    os.putenv('SCAFESRUN_NAME_EXECUTABLE', NAME_EXECUTABLE)

def main():
    global params
    # Check if path to configfile is provided.
    if len(sys.argv) > 1:
        if os.path.isfile(sys.argv[1]) == True:
            params['NAME_CONFIGFILE'] = sys.argv[1]
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

    parseConfigFile()
    calcVariables()
    setEnvironmentVariables()
    subprocess.call("./RUN_HELPER.sh")

if __name__ == '__main__':
    main()
