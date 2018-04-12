import pymc
import numpy as np # use numpy 1.11.3.. newer version break pymc
#import os
#os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'

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
    tumor = pymc.Uniform('w_tumor', 0.05, 0.4, value= 0.15)
    vessel = pymc.Uniform('w_vessel', -2.5, 0, value= -0.5)
    normal = pymc.Uniform('w_normal', 0.13, 40, value= 0.15)

    @pymc.deterministic(plot=False)
    def callScaFES(tumor=tumor, vessel=vessel, normal=normal):
        # set tumor, vessel, normal perfusion to respective values

        # call simulation

        # compute temperatures of normal, tumor, vessel tisue

        T_normal = normal
        T_tumor = tumor
        T_vessel = vessel
        return [T_normal, T_tumor, T_vessel]

    y = pymc.Normal('simulated temperatures', mu=callScaFES, tau=1, value=targetValues, observed=True)

    return locals()


def main():
    # target values for this dataset
    # [T_normal, T_tumor, T_vessel]
    targetValues = [10,0.3,-1]

    # apply MCMC sampler
    MDL = pymc.MCMC(fitSimulation(targetValues))
    MDL.sample(iter=1e4, burn=1000)
    print("\n")
    # extract and plot results
    temperatures = MDL.stats()['callScaFES']['mean']
    tumor = MDL.stats()['w_tumor']['mean']
    vessel = MDL.stats()['w_vessel']['mean']
    normal = MDL.stats()['w_normal']['mean']
    pymc.Matplot.plot(MDL)
    print("T_final: "+str(temperatures))
    print("Perfusion rates: "+str(normal)+"-"+str(tumor)+"-"+str(vessel))
    graph = pymc.graph.graph(MDL)
    graph.write_png("graph.png")

if __name__ == "__main__":
    print(np.__version__)
    main()