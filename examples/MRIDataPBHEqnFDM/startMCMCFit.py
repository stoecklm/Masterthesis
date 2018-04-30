import configparser

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

from helperFunctions import temperature_array_from_result

## Variablen
# Durchblutungsrate (Normal, Tumor, Vessel)
# Waermeuebergangskoeff. (Randbedingung)
# metabolische Waermeproduktion (Normal, Tumor)
# T_A (Normal, Tumor)
# T_Blut?
# Lambda Temperaturleitfaehigkeit
# (\rho c) zusammen gefasst

## Einflussgroessen (unabhaengig von der Simulation)
# Tumortiefe vs. Temperatur an der Oberflaeche

def fitSimulation(targetValues):
    normal = pymc.Uniform('w_normal', 0.0014, 0.014, value= 0.004)
    tumor = pymc.Uniform('w_tumor', 0.0005, 0.017, value= 0.00975)
    vessel = pymc.Uniform('w_vessel', 0.0014, 0.014, value= 0.004)

    @pymc.deterministic(plot=False)
    def callScaFES(normal=normal, tumor=tumor, vessel=vessel):

        # Set normal, tumor, vessel,  perfusion to respective values.
        params = {'NAME_CONFIGFILE' : 'Parameters.ini'}
        params['NAME_RESULTFILE'] = ''
        config = configparser.ConfigParser()
        config.optionxform = str
        config.read(params['NAME_CONFIGFILE'])
        config['Tumor']['OMEGA'] = str(tumor)
        config['Brain']['OMEGA'] = str(normal)
        config['MRI']['USE_VESSELS_SEGMENTATION'] = 'True'
        config['MRI']['VARIABLES_VESSELS'] = 'omega'
        config['MRI']['VALUES_VESSELS'] = str(vessel)
        config['MRI']['VALUES_NON_VESSELS'] = str(normal)

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
        call_simulation(params, 'RUN_HELPER.sh')

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

        return [T_normal, T_tumor, T_vessel]

    y = pymc.Normal('simulated temperatures', mu=callScaFES, tau=1,
                    value=targetValues, observed=True)

    return locals()


def main():
    # Target values for this dataset.
    # [T_normal, T_tumor, T_vessel]
    targetValues = [32.8, 30.0, 34.5]

    # Apply MCMC sampler.
    MDL = pymc.MCMC(fitSimulation(targetValues))
    MDL.sample(iter=500, burn=50)
    print()

    # Extract and plot results.
    temperatures = MDL.stats()['callScaFES']['mean']
    normal = MDL.stats()['w_normal']['mean']
    tumor = MDL.stats()['w_tumor']['mean']
    vessel = MDL.stats()['w_vessel']['mean']
    pymc.Matplot.plot(MDL)
    print('T_final: [T_normal, T_tumor, T_vessel]')
    print('T_final: ', temperatures)
    print('Perfusion rates: Normal:', normal, 'Tumor:', tumor, 'Vessel:', vessel)
    graph = pymc.graph.graph(MDL)
    graph.write_png('graph.png')

    print('Done.')

if __name__ == "__main__":
    print(np.__version__)
    main()