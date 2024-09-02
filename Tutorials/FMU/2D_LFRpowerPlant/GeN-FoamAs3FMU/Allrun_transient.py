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
SecondaryCircuitFMU = load_fmu("modelica/SecondaryCircuitWithoutControl.fmu", log_level=logLevel)
PIDsFMU = load_fmu("simulink/pidControllers.fmu", log_level=logLevel)

# Create the list of FMUs model
models = [ReactorFMU, SecondaryCircuitFMU, PIDsFMU]

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
    (SecondaryCircuitFMU, "electricalPowerOutput", PIDsFMU, "measElectricalPower"),
    (SecondaryCircuitFMU, "frequencyOutput", PIDsFMU, "measFrequency"),
    (PIDsFMU, "pidFrequencyResponse", SecondaryCircuitFMU, "valveControlInput"),
    (PIDsFMU, "externalReactivity", ReactorFMU, "gfExternalReactivity"),
    (PIDsFMU, "externalSourceModulation", ReactorFMU, "gfExternalSourceMod"),
    (ReactorFMU, "gfInnerCorePowerOut", PIDsFMU, "innerCorePower"),
    (ReactorFMU, "gfOuterCorePowerOut", PIDsFMU, "outerCorePower"),
    (ReactorFMU, "gfSteamGenPowerOut1", SecondaryCircuitFMU, "steamGenPower1"),
    (ReactorFMU, "gfSteamGenPowerOut2", SecondaryCircuitFMU, "steamGenPower2"),
    (ReactorFMU, "gfSteamGenPowerOut3", SecondaryCircuitFMU, "steamGenPower3"),
    (ReactorFMU, "gfSteamGenPowerOut4", SecondaryCircuitFMU, "steamGenPower4"),
    (ReactorFMU, "gfSteamGenPowerOut5", SecondaryCircuitFMU, "steamGenPower5"),
    (ReactorFMU, "gfSteamGenPowerOut6", SecondaryCircuitFMU, "steamGenPower6"),
    (ReactorFMU, "gfSteamGenPowerOut7", SecondaryCircuitFMU, "steamGenPower7"),
    (ReactorFMU, "gfSteamGenPowerOut8", SecondaryCircuitFMU, "steamGenPower8")
]

#-----------------------------------------------------------------------------*
# Set variables for GeN-Foam
deltaT = 0.05
startTime = 100
endTime = 200

#SecondaryCircuitFMU.initialize()
#PIDsFMU.initialize()

#-----------------------------------------------------------------------------*
# Reactor
ReactorFMU.set("deltaT", deltaT)
ReactorFMU.set("port", 6059)
ReactorFMU.set("outputPath", "Reactor")

print("End set reactor")

#-----------------------------------------------------------------------------*
# Secondary circuit
# SecondaryCircuitFMU.set("time", startTime)

dataSecCir = pd.read_csv("Reactor/secondaryCircuit.csv")
lastIdx = dataSecCir.index.max()
# Set all the variables for the restart based on the previous simulation

for d in dataSecCir.columns:
    val = np.float64(dataSecCir[d][lastIdx])
    if (d != "time"):
        SecondaryCircuitFMU.set(d, val)

SecondaryCircuitFMU.set("rampLoad.startTime", 110)
SecondaryCircuitFMU.set("rampLoad.height", -20e6)
SecondaryCircuitFMU.set("rampLoad2.startTime", 160)
SecondaryCircuitFMU.set("rampLoad2.height", 20e6)

print("End set secondary circuit")

#-----------------------------------------------------------------------------*
# PIDs
PIDsFMU.set("time", startTime)

dataPIDs = pd.read_csv("Reactor/pidControllers.csv")
lastIdx = dataPIDs.index.max()
PIDsFMU.set("measFrequency", dataPIDs["measFrequency"][lastIdx])
PIDsFMU.set("measElectricalPower", dataPIDs["measElectricalPower"][lastIdx])
PIDsFMU.set("innerCorePower", dataPIDs["innerCorePower"][lastIdx])
PIDsFMU.set("outerCorePower", dataPIDs["outerCorePower"][lastIdx])
PIDsFMU.set("pidFrequencyResponse", dataPIDs["pidFrequencyResponse"][lastIdx])

# Force reset of PID frequency (https://ch.mathworks.com/help/simulink/slref/pidcontroller.html)
PIDsFMU.set("integralInitValue", dataPIDs["integralOutput"][lastIdx])
PIDsFMU.set("resetPidFrequency", 1)


print("End set PID controllers")

#-----------------------------------------------------------------------------*
# Initialize the master simulator
master_simulator = Master(models, connections)
opts = master_simulator.simulate_options()
opts["step_size"] = deltaT
opts["result_handling"] = "file" # csv not supported
# opts["initialize"] = False

#=============================================================================*

# Run the simulation
res = master_simulator.simulate(
    start_time = startTime,
    final_time = endTime,
    options = opts#, input = input_object
)
print("Simulation done")

#=============================================================================*

# Extract the variables into a list
listeOfVariables = [
    'time', 
    'frequencySensor.f', 
    'corePowerMeasure.y', 'injectedPower', 'powerSensorMech.power', 'powerSensorElectric.P',
    'steamGenTemperature1', 'steamGenTemperature2', 'steamGenTemperature3', 'steamGenTemperature4', 
    'steamGenTemperature5', 'steamGenTemperature6', 'steamGenTemperature7', 'steamGenTemperature8',
    'heatSourceSG1.power', 'heatSourceSG2.power', 'heatSourceSG3.power', 'heatSourceSG4.power',
    'heatSourceSG5.power', 'heatSourceSG6.power', 'heatSourceSG7.power', 'heatSourceSG8.power',
    'SG1.h[2]', 'SG2.h[2]', 'SG3.h[2]', 'SG4.h[2]',
    'SG5.h[2]', 'SG6.h[2]', 'SG7.h[2]', 'SG8.h[2]',
    'SG1.p', 'SG2.p', 'SG3.p', 'SG4.p', 'SG5.p', 'SG6.p', 'SG7.p', 'SG8.p',
    'SG1.x[2]', 'SG2.x[2]', 'SG3.x[2]', 'SG4.x[2]',
    'SG5.x[2]', 'SG6.x[2]', 'SG7.x[2]', 'SG8.x[2]',
    'SG1.h[1]', 'SG2.h[1]', 'SG3.h[1]', 'SG4.h[1]',
    'SG5.h[1]', 'SG6.h[1]', 'SG7.h[1]', 'SG8.h[1]',
    'SG1.x[1]', 'SG2.x[1]', 'SG3.x[1]', 'SG4.x[1]',
    'SG5.x[1]', 'SG6.x[1]', 'SG7.x[1]', 'SG8.x[1]',
    'SG1.T[1]', 'SG2.T[1]', 'SG3.T[1]', 'SG4.T[1]',
    'SG5.T[1]', 'SG6.T[1]', 'SG7.T[1]', 'SG8.T[1]',
    'valveVapAdmissionTurb.theta', 'externalReactivity', 'sensW.w'
]
# with open("SecondaryCircuitWithoutControl_1_result.txt") as file:
#     isExtract = False
#     for line in file.readlines():
#         if (len(line) == 1 and isExtract):
#             break
#         if (isExtract):
#             listeOfVariables.append(line.split()[0])
#         if ("char name(" in line):
#             isExtract = True

# PIDsFMU.get_fmu_state() # not supported

# Create the result file
results = pd.DataFrame()
for var in listeOfVariables:
    results[var] = res[SecondaryCircuitFMU][var]
results.to_csv("Reactor/resultsTransient.csv", index=False)

print("End results.csv writting")

#=============================================================================*
