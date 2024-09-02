"""
Script to analyse the csv output from the FMU coupled simulation with GeN-Foam.
Able to render in real time the coupled GeN-Foam/FMU simulation. 

Author: Thomas Guilbaud, EPFL/Transmutex SA, 04/2023
"""

#=============================================================================*
# Imports
#=============================================================================*

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
import os
import re

#=============================================================================*

def listModifier(l: list, scale: float=1, offset: float=0) -> list:
    return([scale*e+offset for e in l])

def extractFromLogPK(filename: str, pattern: str, idx: int=2, tmin: float=0, tmax: float=1e9):
    time = 0
    times, values = [], []
    with open(filename, 'r') as file:
        for line in file.readlines():
            if ("Time = " in line and "ExecutionTime" not in line):
                time = float(line.split()[2])
            elif (pattern in line):
                if (tmin <= time and time <= tmax):
                    times.append(time)
                    values.append(float(line.split()[idx]))
    return(times, values)

def extractFromLog(filename: str, regex: str, offset: float=0, scale: float=1) -> list:
    res = []
    with open(filename, "r") as file:
        res = re.findall(regex, file.read())
    return([scale*float(x)-offset for x in res])

#=============================================================================*

# Test that the user gives an output file in the correct format to be analysed
if (len(sys.argv) != 2 or sys.argv[1].split(".")[-1] != "csv"):
    print('Usage: python Allplot_FMUanimate.py path/to/fmuResults.csv')
    sys.exit(1)

# Create the figure and adjust it
fig, axes = plt.subplots(nrows=2, ncols=3, sharex=True)
fig.subplots_adjust(
    top=0.981,
    bottom=0.078,
    left=0.065,
    right=0.99,
    hspace=0.178,
    wspace=0.211
)

# Create the axes for all plots
axPower, axPID, axMassflowOrPressure, axFrequency, axTemperature, axEnthalpy = axes.flatten()

# Initialize all variable plots
pltPowerCore, = axPower.plot([0], [0], '-', label='Core PK')
# pltPowerCoreMeasure, = axPower.plot([0], [0], '-', label='Core fluid')
pltPowerInje, = axPower.plot([0], [0], '-', label='Steam Generator')
pltPowerMech, = axPower.plot([0], [0], '-', label='Mechanical')
pltPowerElec, = axPower.plot([0], [0], '-', label='Electrical')
pltPIDfreqency, = axPID.plot([0], [0], '-', label='Valve turbine admission [deg]')
# pltTotalReactivity, = axPID.plot([0], [0], '-', label='Total reactivity [pcm]')
pltPIDextReact, = axPID.plot([0], [0], '-', label='External reactivity [pcm]')
pltPIDextSourceMod, = axPID.plot([0], [0], '-', label='External source mod [x10]')
pltMassFlowRate, = axMassflowOrPressure.plot([0], [0], '-', label=r"$\dot{m}$ SG")
# pltPressures = [
#    axMassflowOrPressure.plot([0], [0], '--', label=f"P SG{i}")[0] for i in range(1, 8+1)
# ]
pltPressure, = axMassflowOrPressure.plot([0], [0], '-', label=r"P$_{out}$ SG")
pltFrequency, = axFrequency.plot([0], [0], '-', label="Frequency sensor")
pltTemperatures = [
    axTemperature.plot([0], [0], '-', label=f"SG{i}")[0] for i in range(1, 8+1)
]
# pltTemperatureCoreIn, = axTemperature.plot([0], [0], '--', label="Core in")
# pltTemperatureCoreOut, = axTemperature.plot([0], [0], '--', label="Core out")
pltEnthalpys = [
    axEnthalpy.plot([0], [0], '-', label=f"SG{i}")[0] for i in range(1, 8+1)
]

# Texts
pltPowerCoreText = axPower.text(0, 0, "", va='bottom', ha='center')
pltPowerInjeText = axPower.text(0, 0, "", va='top', ha='center')

# Set axes labels
axTemperature.set_ylabel(r'Temperature [$K$]')
axPower.set_ylabel(r'Power [$MW$]')
axEnthalpy.set_ylabel(r'Specific enthalpy [$MJ/kg$]')
# axSteamQuality.set_ylabel(r'Steam quality [$-$]')
# axPressure.set_ylabel(r'Pressure [$MPa$]')
axPID.set_ylabel(r'PID output [$-$]')
axMassflowOrPressure.set_ylabel(r'Mass flow rate [$kg/s$] or Pressure [$bar$]')
axFrequency.set_ylabel(r'Frequency [$Hz$]')

# Set axes common options
for ax in axes.flatten():
    ax.set_xlabel(r'Time [$s$]')
    ax.grid(True)
    ax.legend()
 
def update(frame):
    """
    Function to update the plot values.
    """
    filename = sys.argv[1]

    # Test if file exists
    if (not os.path.isfile(filename)):
        return

    # Plot from the log
    path = "/".join(filename.split("/")[:-1])
    logFilename = path+"/log.GeN-Foam"
    timeLog, powerPK = extractFromLogPK(logFilename, "totalPower", 2)
    # _, totalReactivity = extractFromLogPK(logFilename, "totalReactivity", 2)
    # Tin  = extractFromLog(logFilename, "faceZone Edge_14_rotated TBulk = (.*) K")
    # Tout = extractFromLog(logFilename, "faceZone Edge_15_rotated TBulk = (.*) K")
    if (len(timeLog) == len(powerPK) and len(timeLog) != 0):
        powerPK = listModifier(powerPK, scale=1e-6 * 360/5)
        pltPowerCore.set_data(timeLog, powerPK)
        pltPowerCoreText.set_position((timeLog[-1], powerPK[-1]))
        pltPowerCoreText.set_text(f"Core={powerPK[-1]:.1f}")

        # pltTotalReactivity.set_data(timeLog, totalReactivity)


    # Read the CSV output file
    data = pd.read_csv(filename)

    # Extract the time
    time = data['time']
    lastTime = list(time)[-1]
 
    # --- Update all the graphs with the new values
    # Power in MW
    powerSG = listModifier(data['injectedPower'], scale=1e-6)
    pltPowerInje.set_data(time, powerSG)
    pltPowerInjeText.set_position((lastTime, powerSG[-1]))
    pltPowerInjeText.set_text(f"SG={powerSG[-1]:.1f}")

    # powerCoreFluid = listModifier(data['corePowerMeasure.y'], scale=1e-6)
    # pltPowerCoreMeasure.set_data(time, powerCoreFluid)
    # pltPowerCoreText.set_text(f"Core fluid={powerCoreFluid[-1]:.1f}")
    # pltPowerCoreText.set_position((lastTime, powerCoreFluid[-1]))

    pltPowerMech.set_data(time, listModifier(data['powerSensorMech.power'], scale=1e-6))
    pltPowerElec.set_data(time, listModifier(data['powerSensorElectric.P'], scale=1e-6))

    # PID outputs in deg or pcm
    pltPIDfreqency.set_data(time, listModifier(data['valveVapAdmissionTurb.theta'], scale=90))
    # pltPIDextReact.set_data(time, listModifier(data['externalReactivity'], scale=1e5))
    # pltPIDextSourceMod.set_data(time, listModifier(data['externalSourceModulation'], scale=10))

    # Frequency in Hz
    pltFrequency.set_data(time, listModifier(data['frequencySensor.f']))

    # Mass flow rate in kg/s
    pltMassFlowRate.set_data(time, listModifier(data['sensW.w']))

    # Pressure in bar
    pltPressure.set_data(time, listModifier(data[f'SG8.p'], scale=1e-5))
    # for i, pltPressure in enumerate(pltPressures):
    #     pltPressure.set_data(time, listModifier(data[f'SG{i+1}.p'], scale=1e-5))

    # Temperature in K
    for i, pltTemperature in enumerate(pltTemperatures):
        pltTemperature.set_data(time, listModifier(data[f'steamGenTemperature{i+1}']))
    # pltTemperatureCoreIn.set_data(np.arange(timeLog[0], timeLog[-1], 1), Tin)
    # pltTemperatureCoreOut.set_data(np.arange(timeLog[0], timeLog[-1], 1), Tout)

    # Specific enthalpy in MJ/kg
    for i, pltEnthalpy in enumerate(pltEnthalpys):
        pltEnthalpy.set_data(time, listModifier(data[f'SG{i+1}.h[2]'], scale=1e-6))
    
    # Resize all the plots
    for ax in axes.flatten():
        ax.relim()
        ax.autoscale()

    # axFrequency.set_ylim((49.97, 50.03))
    # axPower.set_ylim(0)
 
# Main function for animation
animation = FuncAnimation(fig, update, interval=1000)

# Plot rendering
plt.show()

#=============================================================================*