"""
Script to analyse the csv output from the FMU coupled simulation with GeN-Foam.

Author: Thomas Guilbaud, EPFL/Transmutex SA, 12/2022
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
    with open(filename, 'r') as f:
        for line in f:
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
    print('Usage: python Allplot_FMU.py path/to/fmuResults.csv')
    sys.exit(1)


# Create figures
fig1, axes1 = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(9, 6), dpi=300)
fig2, axes2 = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(9, 6), dpi=300)

# Linestyle array -> 3 input files maximum
linestyles = ["-", "--", "-."]

for k, filename, ls in zip(range(len(sys.argv[1:])), sys.argv[1:], linestyles):
    print(f"Analyze {filename}")

    # Test output file is in csv format
    if (filename.split('.')[-1] != 'csv'):
        print(f'Incorrect file extension {filename} ! The file must be in csv format.')
        print('Usage: python Allplot_FMU.py path/to/fmuResults.csv')
        sys.exit(1)

    # Read output file
    data = pd.read_csv(filename)

    # Extract time list
    time = data['time']

    print(f'Simulation time = {list(time)[-1]:.2f} s = {list(time)[-1]/60:.2f} min')

    # Create figure for plot
    axTemperature, axPower, axEnthalpy, axSteamQuality = axes1.flatten()
    axPressure, axPID, axMassflowrate, axFrequency = axes2.flatten()

    powerTotGFInput = []
    powerInje = listModifier(data['injectedPower'], scale=1e-6)
    powerMech = listModifier(data['powerSensorMech.power'], scale=1e-6)
    powerElec = listModifier(data['powerSensorElectric.P'], scale=1e-6)

    # Read parameters from the different sections of the steam generator
    for i in range(1, 8+1):
        # Extract values from SG sections
        T = list(data[f'steamGenTemperature{i}'])
        power = listModifier(data[f'heatSourceSG{i}.power'], scale=1e-6)
        enthalpy = listModifier(data[f'SG{i}.h[2]'], scale=1e-6)
        pressure = listModifier(data[f'SG{i}.p'], scale=1e-6)
        steamQuality = list(data[f'SG{i}.x[2]'])

        # Merge the individual sections of SG power to a have a total power
        if (i == 1):
            powerTotGFInput = power
        else:
            powerTotGFInput = list(map(lambda x, y: x+y, powerTotGFInput, power))

        # Test only for the first section of the SG
        if (i == 1):
            inletEnthalpy = listModifier(data[f'SG{i}.h[1]'], scale=1e-6)
            inletSteamQuality = list(data[f'SG{i}.x[1]'])
            inletTemperature = list(data[f'SG{i}.T[1]'])
            print(f'SG{i} {inletTemperature[-1]:10.2f} K {inletEnthalpy[-1]:24.3f} MJ/kg {inletSteamQuality[-1]:10.4f}')

        # Print info per SG sections
        print(
            f'SG{i}',
            f'{T[-1]:10.2f} K',
            f'{power[-1]/56.70*5/360:10.5f} MW',
            f'{enthalpy[-1]:10.3f} MJ/kg',
            f'{steamQuality[-1]:10.4f}',
            f'{pressure[-1]:10.4f} MPa'
        )

        steanGenLabel = f'SG{i}' if k==0 else ""

        # Plot
        axTemperature.plot(time, T, ls=ls, label=steanGenLabel)
        axPressure.plot(time, pressure, ls=ls, label=steanGenLabel)
        # axPower.plot(time, power, ls=ls, label=steanGenLabel)
        axEnthalpy.plot(time, enthalpy, ls=ls, label=steanGenLabel)
        axSteamQuality.plot(time, steamQuality, ls=ls, label=steanGenLabel)


    # Print global parameters at the end of the simulation
    print(f'Total GeN-Foam power = {powerTotGFInput[-1]:.5f} MW')
    print(f'Mechanical     power = {powerMech[-1]:.5f} MW')
    print(f'Electrical     power = {powerElec[-1]:.5f} MW')
    print(f"Target inlet  temperature {400+273.15} K")
    print(f"Target outlet temperature {480+273.15} K")


    path = "/".join(filename.split("/")[:-1])
    logFilename = path+"/log.GeN-Foam"
    powerSG = 0
    # for i in range(8):
    #     powerSG += extractFromLog(logFilename, f"volIntegrate\(SGenPrim{i+1}\) of multiply\(htc,subtract\(T,T.fixedTemperatureFMU\)\) = (.*)")[-1]
    print("From the log:")
    print("  Power core = {:.6f} MW".format(
        extractFromLog(logFilename, "volIntegrate\(fluidRegion\) of powerDensityNeutronics = (.*)")[-1]/1e6
    ))
    print(f"  Power SG   = {powerSG/1e6:.6f} MW")

    Si = extractFromLog(logFilename, "faceZone Edge_14_rotated massFlow = .* kg/s over (.*) m2")[-1]
    So = extractFromLog(logFilename, "faceZone Edge_35_rotated massFlow = .* kg/s over (.*) m2")[-1]

    innerCoreInletTemp  = extractFromLog(logFilename, "faceZone Edge_13_rotated TBulk = (.*) K")[-1]
    outerCoreInletTemp  = extractFromLog(logFilename, "faceZone Edge_34_rotated TBulk = (.*) K")[-1]
    innerCoreOutletTemp = extractFromLog(logFilename, "faceZone Edge_15_rotated TBulk = (.*) K")[-1]
    outerCoreOutletTemp = extractFromLog(logFilename, "faceZone Edge_36_rotated TBulk = (.*) K")[-1]

    coreInletTemp = (innerCoreInletTemp*Si + outerCoreInletTemp*So)/(Si+So)
    coreOutletTemp = (innerCoreOutletTemp*Si + outerCoreOutletTemp*So)/(Si+So)

    print(f"  Tin core    = {coreInletTemp:.3f} K")
    print(f"  Tout core   = {coreOutletTemp:.3f} K")

    timeLog, powerPK = extractFromLogPK(path+"/log.GeN-Foam", "totalPower", 2)
    powerPK = listModifier(powerPK, scale=72*1e-6)

    # Plot other intersting parameters outside of the steam genarator
    axPower.plot(timeLog, powerPK, ls=ls, label='Core PK' if k==0 else "")
    axPower.plot(time, powerTotGFInput, ls=ls, label='Thermal in SG' if k==0 else "")
    axPower.plot(time, powerMech, ls=ls, label="Mechanic" if k==0 else "")
    axPower.plot(time, powerElec, ls=ls, label="Electric" if k==0 else "")
    # axPower.plot(time, powerInje, ls=ls, label="Injected")

    axPID.plot(time, listModifier(data['valveVapAdmissionTurb.theta'], scale=90), ls=ls, label='Valve turbine admission [deg]' if k==0 else "")

    # axRate.plot(time, listModifier(data['derValveAdmission.y'], scale=90), ls=ls, label='Valve turbine admission [deg/s]' if k==0 else "")

    axMassflowrate.plot(time, data['sensW.w'], ls=ls, label="Outlet SG" if k==0 else "")

    axFrequency.plot(time, data['frequencySensor.f'], ls=ls, label="Frequency sensor" if k==0 else "")


    # Plot vertical labels
    axTemperature.set_ylabel(r'Temperature [$K$]')
    axPower.set_ylabel(r'Power [$MW$]')
    axEnthalpy.set_ylabel(r'Specific enthalpy [$MJ/kg$]')
    axSteamQuality.set_ylabel(r'Steam quality [$-$]')
    axPressure.set_ylabel(r'Pressure [$MPa$]')
    axPID.set_ylabel(r'PID output [$-$]')
    # axRate.set_ylabel(r'Rate [$-$]')
    axMassflowrate.set_ylabel(r'Mass flow rate [$kg/s$]')
    axFrequency.set_ylabel(r'Frequency [$Hz$]')

    # Plot axis ranges
    # axPressure.set_xlim((0, 300))
    # axPower.set_ylim(0)

    # Plotting options on all the subplots
    for ax in list(axes1.flatten())+list(axes2.flatten()):
        ax.legend()
        ax.grid(True)
        ax.set_xlabel(r"Time [$s$]")
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)

        # Reset color wheel for the next plots
        ax.set_prop_cycle(None)

fig1.tight_layout()
fig2.tight_layout()

fig1.savefig("fig_results_allplot_fmu1.png")
fig2.savefig("fig_results_allplot_fmu2.png")

# plt.show()

#=============================================================================*