"""
Script to analyse the load follow transient.

Author: Thomas Guilbaud, EPFL/Transmutex SA, 24/04/2023
"""

#=============================================================================*
# Imports
#=============================================================================*

import pandas as pd
import matplotlib.pyplot as plt
import sys
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

# Test give an output file to analyse
if (len(sys.argv) <= 1):
    print(f'\n    Usage: python3 {sys.argv[0]} path/to/fmuResults.csv\n')
    sys.exit(1)

# timeTransient=1000
t0 = 1000
wedge = 360/5

# Create figures
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(12, 12), dpi=400)


# Linestyle array -> 4 input files maximum
linestyles = ["-", "--", "-.", ":"]

for k, filename, ls in zip(range(len(sys.argv[1:])), sys.argv[1:], linestyles):
    print(f"Analyze {filename}")

    # Test output file is in csv format
    if (filename.split('.')[-1] != 'csv'):
        print(f'Incorrect file extension {filename} ! The file must be in csv format.')
        print(f'Usage: python3 {sys.argv[0]} path/to/fmuResults.csv')
        sys.exit(1)

    # Read output file
    data = pd.read_csv(filename)

    # Extract log file name from path
    path = "/".join(filename.split("/")[:-1])
    logFilename = path+"/log.GeN-Foam"
    
    # Extract values from the log
    massFlowSG   = extractFromLog(logFilename, "faceZone Edge_102_rotated massFlow = (.*) kg/s")
    inletSGTemp  = extractFromLog(logFilename, "faceZone Edge_94_rotated TBulk = (.*) K")
    outletSGTemp = extractFromLog(logFilename, "faceZone Edge_98_rotated TBulk = (.*) K")
    corePowerPK = extractFromLog(logFilename, "totalPower = (.*) W", scale=wedge/1e6)
    corePowerVolInt = extractFromLog(logFilename, "volIntegrate\(fluidRegion\) of powerDensityNeutronics = (.*)", scale=1e-6)

    innerCoreSurface = extractFromLog(logFilename, "faceZone Edge_14_rotated massFlow = .* kg/s over (.*) m2")
    outerCoreSurface = extractFromLog(logFilename, "faceZone Edge_35_rotated massFlow = .* kg/s over (.*) m2")

    innerCoreInletTemp  = extractFromLog(logFilename, "faceZone Edge_13_rotated TBulk = (.*) K")
    outerCoreInletTemp  = extractFromLog(logFilename, "faceZone Edge_34_rotated TBulk = (.*) K")
    innerCoreOutletTemp = extractFromLog(logFilename, "faceZone Edge_15_rotated TBulk = (.*) K")
    outerCoreOutletTemp = extractFromLog(logFilename, "faceZone Edge_36_rotated TBulk = (.*) K")

    coreInletTemp = [(Ti*Si + To*So)/(Si+So) for Ti, To, Si, So in zip(
        innerCoreInletTemp, outerCoreInletTemp, innerCoreSurface, outerCoreSurface
    )]
    coreOutletTemp = [(Ti*Si + To*So)/(Si+So) for Ti, To, Si, So in zip(
        innerCoreOutletTemp, outerCoreOutletTemp, innerCoreSurface, outerCoreSurface
    )]


    # Create figure for plot
    axTemperature, axPower, axReactivity, axBeamIntensity, axTurbineValve, axFrequency = axes.flatten()

    # Extract reactivity
    timeLog, totalReactivity = extractFromLogPK(logFilename, "totalReactivity", tmin=t0)
    axReactivity.plot(
        listModifier(timeLog, offset=-t0),
        totalReactivity,
        label='Total' if k == 0 else "",
        ls=ls
    )
        
    timeLog, dopplerReactivity = extractFromLogPK(logFilename, "Doppler", 4, tmin=t0)
    axReactivity.plot(
        listModifier(timeLog, offset=-t0),
        dopplerReactivity,
        label='Doppler' if k == 0 else "",
        ls=ls
    )
    
    timeLog, TStructReactivity = extractFromLogPK(logFilename, "-> TStruct ", 3, tmin=t0)
    axReactivity.plot(
        listModifier(timeLog, offset=-t0),
        TStructReactivity,
        label='Radial exp' if k == 0 else "",
        ls=ls
    )
    
    timeLog, TFuelReactivity = extractFromLogPK(logFilename, "-> TFuel", 3, tmin=t0)
    axReactivity.plot(
        listModifier(timeLog, offset=-t0),
        TFuelReactivity,
        label='Axial exp' if k == 0 else "",
        ls=ls
    )
    
    timeLog, TCoolReactivity = extractFromLogPK(logFilename, "-> TCool", 3, tmin=t0)
    axReactivity.plot(
        listModifier(timeLog, offset=-t0),
        TCoolReactivity,
        label='Coolant' if k == 0 else "",
        ls=ls
    )
    
    timeLog, TCladReactivity = extractFromLogPK(logFilename, "-> TClad", 3, tmin=t0)
    axReactivity.plot(
        listModifier(timeLog, offset=-t0),
        TCladReactivity,
        label='Cladding' if k == 0 else "",
        ls=ls
    )
    
    timeLog, drivelineReactivity = extractFromLogPK(logFilename, "-> driveline", 3, tmin=t0)
    axReactivity.plot(
        listModifier(timeLog, offset=-t0),
        drivelineReactivity,
        label='Driveline' if k == 0 else "",
        ls=ls
    )
    
    timeLog, externalReactivity = extractFromLogPK(logFilename, "-> extReactivity  =", 3, tmin=t0)
    if (len(externalReactivity) != 0):
        axReactivity.plot(
            listModifier(timeLog, offset=-t0),
            externalReactivity,
            label='Ext FMU' if k == 0 else "",
            ls=ls
        )


    # Extract time list
    time = listModifier(data['time'], offset=-t0)
    idxStart = time.index(0)
    time = time[idxStart:]

    # Extract powers
    powerInje = listModifier(data['injectedPower'], scale=1e-6)
    powerMech = listModifier(data['powerSensorMech.power'], scale=1e-6)
    powerElec = listModifier(data['powerSensorElectric.P'], scale=1e-6)

    # Plot other interesting parameters outside of the steam generator
    # axPower.plot(corePowerVolInt, ls=ls, label="Core" if k == 0 else "")
    axPower.plot(time, powerInje[idxStart:], ls=ls, label='SG' if k == 0 else "", marker='o', markevery=100, markersize=5)
    axPower.plot(time, powerMech[idxStart:], ls=ls, label="Mechanic" if k == 0 else "", marker='x', markevery=100, markersize=5)
    axPower.plot(time, powerElec[idxStart:], ls=ls, label="Electric" if k == 0 else "", marker='s', markevery=100, markersize=5)

    axTemperature.plot(coreInletTemp, label="Core inlet" if k == 0 else "", ls=ls)
    axTemperature.plot(coreOutletTemp, label="Core outlet" if k == 0 else "", ls=ls, marker='o', markevery=100, markersize=5)
    axTemperature.plot(time, list(data[f'steamGenTemperature1'])[idxStart:], label=f'SG water inlet' if k == 0 else "", ls=ls, marker='x', markevery=100, markersize=5)
    axTemperature.plot(time, list(data[f'steamGenTemperature8'])[idxStart:], label=f'SG steam outlet' if k == 0 else "", ls=ls, marker='s', markevery=100, markersize=5)

    axFrequency.plot(time, list(data['frequencySensor.f'])[idxStart:], ls=ls, label="Frequency sensor" if k == 0 else "")
    axTurbineValve.plot(time, list(data['valveVapAdmissionTurb.theta'])[idxStart:], ls=ls, label="Turbine valve admission" if k == 0 else "")
    axBeamIntensity.plot(time, list(data['externalSourceModulation'])[idxStart:], ls=ls, label="Source modulation" if k == 0 else "")

    # Plot vertical labels
    axTemperature.set_ylabel(r'Temperature [$K$]')
    axPower.set_ylabel(r'Power [$MW$]')
    axReactivity.set_ylabel(r'Reactivity [$pcm$]')
    axFrequency.set_ylabel(r'Frequency [$Hz$]')
    axTurbineValve.set_ylabel(r'Valve orientation [$deg$]')
    axBeamIntensity.set_ylabel(r'Source modulation')

    # Plotting options on all the subplots
    for ax in list(axes.flatten()):
        ax.legend(loc='center right')
        ax.grid(True)
        ax.set_xlabel(r"Time [$s$]")
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)

        # Reset color wheel for the next plots
        ax.set_prop_cycle(None)

fig.tight_layout()

fig.savefig("results_allplot_loadFollowSource.png")

# plt.show()

#=============================================================================*
