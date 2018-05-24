import configparser
import glob
import os
import sys

import matplotlib
matplotlib.use('Agg')
import pymc
import numpy as np # use numpy 1.11.3.. newer version break pymc

from startSimulation import parse_config_file
from startSimulation import check_variables
from startSimulation import calc_variables
from startSimulation import check_stability
from startSimulation import create_region_file
from startSimulation import create_init_file
from startSimulation import set_environment_variables
from startSimulation import call_simulation

from postProcessing import calc_tumor_near_surface_temperatures
from postProcessing import calc_vessels_temperatures
from postProcessing import calc_non_vessels_temperatures
from postProcessing import region_array_from_file
from postProcessing import surface_vessels_array_from_file

from helperFunctions import close_nc_file
from helperFunctions import create_mcmc_netcdf_file
from helperFunctions import create_testcase_name
from helperFunctions import parse_pymc_from_config_file
from helperFunctions import save_vector_to_mcmc_file
from helperFunctions import temperature_array_from_result
from helperFunctions import write_ini_file_to_nc_file


## Einflussgroessen (unabhaengig von der Simulation)
# Tumortiefe vs. Temperatur an der Oberflaeche

count = 0
params = {'NAME_CONFIGFILE_TEMPLATE' : ''}
TESTED_VARIABLES = 'vessels'

def fitSimulation(targetValues):
    #omega_brain = pymc.Uniform('omega_brain', 0.0014, 0.014, value=0.004)
    omega_brain = pymc.Uniform('omega_brain', 0.001, 0.017, value=0.004)
    #omega_tumor = pymc.Uniform('omega_tumor', 0.0005, 0.017, value=0.00975)
    omega_tumor = pymc.Uniform('omega_tumor', 0.0003, 0.020, value=0.00975)
    #omega_vessel = pymc.Uniform('omega_vessel', 0.0014, 0.014, value=0.004)
    omega_vessel = pymc.Uniform('omega_vessel', 0.001, 0.017, value=0.001)

    @pymc.deterministic(plot=False)
    def callScaFES(omega_brain=omega_brain,
                   omega_tumor=omega_tumor,
                   omega_vessel=omega_vessel):

        global count
        global params
        count += 1

        print()
        print('##### ScaFES iteration: {} #####'.format(count))

        case = params['NAME_CONFIGFILE_TEMPLATE'].split('.')[0]
        case = case.split('_')[0]
        params['NAME_CONFIGFILE'] = case + '-pymc-' + TESTED_VARIABLES + '.ini'
        params['NAME_RESULTFILE'] = ''
        config = configparser.ConfigParser()
        config.optionxform = str
        config.read(params['NAME_CONFIGFILE_TEMPLATE'])
        # Omega.
        config['Brain']['OMEGA'] = str(omega_brain)
        config['Tumor']['OMEGA'] = str(omega_tumor)
        config['MRI']['USE_VESSELS_SEGMENTATION'] = 'True'
        config['MRI']['VARIABLES_VESSELS'] = 'omega'
        config['MRI']['VALUES_VESSELS'] = str(omega_vessel)
        config['MRI']['VALUES_NON_VESSELS'] = str(omega_brain)

        with open(params['NAME_CONFIGFILE'], 'w') as configfile:
            config.write(configfile)

        # Call simulation.
        parse_config_file(params)
        check_variables(params)
        calc_variables(params)
        check_stability(params)
        if params['USE_MRI_FILE'] == False:
            create_region_file(params)
        create_init_file(params)
        set_environment_variables(params)
        call_simulation(params, params['RUN_SCRIPT'])
        if params['NAME_RESULTFILE'] != '':
            # Compute temperatures of normal, tumor, vessel tisue.
            temp = temperature_array_from_result(params['NAME_RESULTFILE'])
            tumor = region_array_from_file(params['NAME_REGION_FILE'])
            vessels = surface_vessels_array_from_file(params['NAME_VESSELS_FILE'])

            T_tumor,_,_,_ = calc_tumor_near_surface_temperatures(temp, tumor)
            if params['USE_VESSELS_SEGMENTATION'] == True:
                temp = temp[-1,:,:]
                T_vessel,_,_,_ = calc_vessels_temperatures(temp, vessels)
                T_normal,_,_,_ = calc_non_vessels_temperatures(temp, vessels)
            else:
                T_vessel = -1.0
                T_normal = -1.0
        else:
            print('* ERROR: No result file written.')
            print('Aborting.')
            exit()

        global dat_file_name
        with open(dat_file_name, 'a') as backupfile:
            backupfile.write(str(omega_brain) + '\t' \
                             + str(omega_tumor) + '\t' \
                             + str(omega_vessel) + '\t' \
                             + str(T_normal) + '\t' \
                             + str(T_tumor) + '\t' \
                             + str(T_vessel) + '\n')

        if (count % 25 == 0):
            searchpath = './' + params['NAME_EXECUTABLE'] + '-' + case + '-*.*'
            log_files = glob.glob(searchpath)
            for f in log_files:
                os.remove(f)

        return [T_normal, T_tumor, T_vessel]

    y = pymc.Normal('simulated temperatures', mu=callScaFES, tau=1,
                    value=targetValues, observed=True)

    return locals()


def main():
    print('Starting MCMC simulation and fit.')
    global params
    # Check if path to configfile is provided and if file exists.
    if len(sys.argv) > 1:
        if os.path.isfile(sys.argv[1]) == True:
            params['NAME_CONFIGFILE_TEMPLATE'] = sys.argv[1]
            params['RUN_SCRIPT'] = 'RUN_HELPER.sh'
        else:
            print('* ERROR:', sys.argv[1], 'does not exist.')
        if len(sys.argv) > 2:
            if os.path.isfile(sys.argv[2]) == True:
                params['RUN_SCRIPT'] = sys.argv[2]
            else:
                print('* ERROR: Optional run script', sys.argv[2],
                      'does not exist.')
                print('Aborting.')
                exit()
    else:
        print('* ERROR: No command line argument for configfile provided.')

    if params['NAME_CONFIGFILE_TEMPLATE'] == '':
        print('Usage: python3', sys.argv[0],
              '<PATH/TO/CONFIGFILE> [<PATH/TO/RUN/SCRIPT]')
        print('Aborting.')
        exit()

    name = create_testcase_name(TESTED_VARIABLES, params)
    db_name = name + '.pickle'
    parse_pymc_from_config_file(params)

    sample_iterations = params['ITERATIONS']
    sample_burns = params['BURNS']
    print('Number of sample iterations: {}.'.format(sample_iterations))
    print('Number of sample burns: {}.'.format(sample_burns))
    targetValues = [params['T_NORMAL'], params['T_TUMOR'], params['T_VESSEL']]
    print('Target values: [T_normal, T_tumor, T_vessel]')
    print('Target values for this dataset: {}.'.format(targetValues))

    # Create backup file.
    global dat_file_name
    dat_file_name = name + '.dat'
    with open(dat_file_name, 'w') as backupfile:
        backupfile.write('#omega_brain\t' + \
                         'omega_tumor\t' + \
                         'omega_vessel\t' + \
                         'T_normal\t' + \
                         'T_tumor\t' + \
                         'T_vessel\n')

    # Apply MCMC sampler.
    MDL = pymc.MCMC(fitSimulation(targetValues), db='pickle',
                    dbname=db_name)
    MDL.sample(iter=sample_iterations, burn=sample_burns)
    print()

    # Extract and plot results.
    temperatures = MDL.stats()['callScaFES']['mean']
    omega_brain = MDL.stats()['omega_brain']['mean']
    omega_tumor = MDL.stats()['omega_tumor']['mean']
    omega_vessel = MDL.stats()['omega_vessel']['mean']
    pymc.Matplot.plot(MDL)
    print()
    print('T_final: [T_normal, T_tumor, T_vessel]')
    print('T_final:', temperatures)
    print()
    print('omega_brain:', omega_brain)
    print('omega_tumor:', omega_tumor)
    print('omega_vessel:', omega_vessel)

    graph = pymc.graph.graph(MDL)
    graph.write_png('graph.png')

    print()
    print('Number of ScaFES calls:', count)
    print()

    T_normal = MDL.trace('callScaFES')[:,0]
    T_tumor = MDL.trace('callScaFES')[:,1]
    T_vessel = MDL.trace('callScaFES')[:,2]
    omega_brain = MDL.trace('omega_brain')[:]
    omega_tumor = MDL.trace('omega_tumor')[:]
    omega_vessel = MDL.trace('omega_vessel')[:]

    l2_norm = np.linalg.norm(np.subtract(MDL.trace('callScaFES')[:],
                                         targetValues), 2, axis=1)
    iterations = l2_norm.shape[0]
    nc_file = create_mcmc_netcdf_file('pymc_' + name + '.nc', iterations)
    save_vector_to_mcmc_file(nc_file, omega_brain, 'omega_brain')
    save_vector_to_mcmc_file(nc_file, omega_tumor, 'omega_tumor')
    save_vector_to_mcmc_file(nc_file, omega_vessel, 'omega_vessel')
    save_vector_to_mcmc_file(nc_file, l2_norm, 'L2-norm')
    save_vector_to_mcmc_file(nc_file, T_normal, 'T_normal')
    save_vector_to_mcmc_file(nc_file, T_tumor, 'T_tumor')
    save_vector_to_mcmc_file(nc_file, T_vessel, 'T_vessel')
    write_ini_file_to_nc_file(nc_file, params['NAME_CONFIGFILE_TEMPLATE'])
    close_nc_file(nc_file)

    print('Backup file written to {}.'.format(dat_file_name))

    MDL.db.close()
    print()
    print('Done.')

if __name__ == "__main__":
    main()
