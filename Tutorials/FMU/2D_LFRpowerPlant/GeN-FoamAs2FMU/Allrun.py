"""
Author: Thomas Guilbaud, EPFL/Transmutex, 2023/05
"""

#=============================================================================*
# Imports

from pyfmi import load_fmu
from pyfmi.master import Master
import pandas as pd
import numpy as np
import sys

#=============================================================================*

# Load the FMUs
logLevel = 2
ReactorFMU = load_fmu("Reactor.fmu", log_level=logLevel)
SecondaryCircuitFMU = load_fmu("modelica/SecondaryCircuitWithControl.fmu", log_level=logLevel)

# Create the list of FMUs model
models = [ReactorFMU, SecondaryCircuitFMU]

# Create the link map (output to input)
connections = [ 
    (SecondaryCircuitFMU, "steamGenTemperature1", ReactorFMU, "gfSteamGenTemperature1"),
    (SecondaryCircuitFMU, "steamGenTemperature2", ReactorFMU, "gfSteamGenTemperature2"),
    (SecondaryCircuitFMU, "steamGenTemperature3", ReactorFMU, "gfSteamGenTemperature3"),
    (SecondaryCircuitFMU, "steamGenTemperature4", ReactorFMU, "gfSteamGenTemperature4"),
    (SecondaryCircuitFMU, "steamGenTemperature5", ReactorFMU, "gfSteamGenTemperature5"),
    (SecondaryCircuitFMU, "steamGenTemperature6", ReactorFMU, "gfSteamGenTemperature6"),
    (SecondaryCircuitFMU, "steamGenTemperature7", ReactorFMU, "gfSteamGenTemperature7"),
    (SecondaryCircuitFMU, "steamGenTemperature8", ReactorFMU, "gfSteamGenTemperature8"),
    (SecondaryCircuitFMU, "externalReactivity", ReactorFMU, "gfExternalReactivity"),
    (SecondaryCircuitFMU, "externalSourceModulation", ReactorFMU, "gfExternalSourceMod"),
    (ReactorFMU, "gfInnerCorePowerOut", SecondaryCircuitFMU, "innerCorePower"),
    (ReactorFMU, "gfOuterCorePowerOut", SecondaryCircuitFMU, "outerCorePower"),
    (ReactorFMU, "gfSteamGenPowerOut1", SecondaryCircuitFMU, "steamGenPower1"),
    (ReactorFMU, "gfSteamGenPowerOut2", SecondaryCircuitFMU, "steamGenPower2"),
    (ReactorFMU, "gfSteamGenPowerOut3", SecondaryCircuitFMU, "steamGenPower3"),
    (ReactorFMU, "gfSteamGenPowerOut4", SecondaryCircuitFMU, "steamGenPower4"),
    (ReactorFMU, "gfSteamGenPowerOut5", SecondaryCircuitFMU, "steamGenPower5"),
    (ReactorFMU, "gfSteamGenPowerOut6", SecondaryCircuitFMU, "steamGenPower6"),
    (ReactorFMU, "gfSteamGenPowerOut7", SecondaryCircuitFMU, "steamGenPower7"),
    (ReactorFMU, "gfSteamGenPowerOut8", SecondaryCircuitFMU, "steamGenPower8")
]

# Set variables for GeN-Foam
deltaT = 0.05
ReactorFMU.set("deltaT", deltaT)
ReactorFMU.set("port", 6053)
ReactorFMU.set("outputPath", "Reactor")

# Set variables for BoP
SecondaryCircuitFMU.set("rampLoad.startTime", 100)
SecondaryCircuitFMU.set("rampLoad.height", -20e6)
SecondaryCircuitFMU.set("rampLoad2.startTime", 600)
SecondaryCircuitFMU.set("rampLoad2.height", 20e6)
SecondaryCircuitFMU.set("stepActivateExtReactivity.startTime", 0)

# Initialize the master simulator
master_simulator = Master(models, connections)
opts = master_simulator.simulate_options()
opts["step_size"] = deltaT
opts["result_handling"] = "file" # csv not supported
# opts["initialize"] = True

#=============================================================================*

# Run the simulation
res = master_simulator.simulate(
    start_time = 0,
    final_time = 1000,
    options = opts
)
print("Simulation done")

#=============================================================================*

# Extract the variables into a list

def getListVars(filename: str) -> list:
    listVar = []
    with open(filename, 'r') as file:
        isExtract = False
        for line in file.readlines():
            if (len(line) == 1 and isExtract):
                break
            if (isExtract):
                listVar.append(line.split()[0])
            if ("char name(" in line):
                isExtract = True
    return(listVar)

listOfVariablesSecCir = getListVars("SecondaryCircuitWithControl_1_result.txt")

# Create the result file
resultsSecCir = pd.DataFrame()
for var in listOfVariablesSecCir:
    resultsSecCir[var] = res[SecondaryCircuitFMU][var]
resultsSecCir.to_csv("Reactor/secondaryCircuit.csv", index=False)

print("End Reactor/secondaryCircuit.csv writting")


#=============================================================================*
