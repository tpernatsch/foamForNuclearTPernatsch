"""
Post-processing script to plot power, mass flow rate, temperatures and
reactivity evolution over time.
"""

#==============================================================================*
### IMPORTS

import sys
import matplotlib.pyplot as plt


#==============================================================================*
### USAGE

if (len(sys.argv) == 1):
    print(f"\n    Usage: python3 {sys.argv[0]} path1/to/log.GeN-Foam path2/to/log.GeN-Foam ...\n")
    sys.exit(0)


#==============================================================================*
### SETTINGS

beta = 0.00313126 * 1e5

SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE-1)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


#==============================================================================*
### FUNCTIONS

def readTimes(loglines) -> list:
    times = []
    totalRunTime = 0
    for line in loglines:
        if "Time =" in line[:6]:
            times.append(float(line.split()[2]))
        if "ExecutionTime = " in line:
            totalRunTime = float(line.split()[2])
    print(f"    Total run time = {totalRunTime} s")
    return(times)


def readValues(times, loglines, keyword: str, pos: int, scale: float=1) -> list:
    values = []
    for line in loglines:
        if (keyword in line):
            values.append(scale * float(line.split()[pos]))

            if (len(values) == len(times)):
                return(values)
    return(values)


def readExpTimesAndValues(loglines, timeOffset: float, valueOffset: float, valueScale: float=1):
    times, values = [], []
    for line in loglines:
        try:
            times.append(float(line.split()[0]) + timeOffset)
            values.append(valueScale * float(line.split()[1]) + valueOffset)
        except:
            pass
    return(times, values)


def matchSize(list1: list, list2: list) -> None:
    while len(list1) > len(list2):
        del list1[-1]
    while len(list2) > len(list1):
        del list2[-1]


#==============================================================================*
### EXTRACT AND PLOT DATA

fig, axes = plt.subplots(2, 2, sharex=True, figsize=(10, 7))

axP, axF, axT, axR = axes.flatten()

for ax in axes.flatten():
    ax.set_xlabel(r'$t\;(s)$')
    ax.grid(True)

axP.set_ylabel(r'$P\;(W)$')
axF.set_ylabel(r'$\dot{m}\;(kg/s)$')
axT.set_ylabel(r'$T\;(K)$')
axR.set_ylabel(r'$\rho\;(\$)$')
axP.set_yscale('log')
axF.set_yscale('log')

for logname, ls in zip(sys.argv[1:], ['-', '-.', ':']):
    # Reset color wheel
    for ax in axes.flatten():
        ax.set_prop_cycle(None)  # Reset the property cycle

    print(f"Processing {logname} ...")

    with open(logname, "r") as log:
        loglines = log.readlines()

        times = readTimes(loglines)

        for i in range(10):
            del times[-1]

        flowPrimary = readValues(times, loglines, "faceZone pumpMiddleCutPrimary massFlow", 4)
        flowCore = readValues(times, loglines, "faceZone coreMiddleCut massFlow", 4)
        flowInnerCore = readValues(times, loglines, "faceZone innerCoreMiddleCut massFlow", 4)
        flowOuterCore = readValues(times, loglines, "faceZone outerCoreMiddleCut massFlow", 4)
        flowColdLegSecondaryInlet = readValues(times, loglines, "patch coldLegSecondaryInlet massFlow", 4)


        axF.plot(times, flowPrimary, label="primary", ls=ls)
        axF.plot(times, flowCore, label="core", ls=ls)
        axF.plot(times, flowInnerCore, label="innerCore", ls=ls)
        axF.plot(times, flowOuterCore, label="outerCore", ls=ls)

        TCoreIn = readValues(times, loglines, "faceZone coreInlet TBulk", 4)
        TCoreOut = readValues(times, loglines, "faceZone coreOutlet TBulk", 4)
        TInnerCoreOut = readValues(times, loglines, "faceZone innerCoreOutlet TBulk", 4)
        TOuterCoreOut = readValues(times, loglines, "faceZone outerCoreOutlet TBulk", 4)
        TIHXInPrimary = readValues(times, loglines, "faceZone primaryIHXInlet TBulk", 4)
        TIHXOutPrimary = readValues(times, loglines, "faceZone primaryIHXOutlet TBulk", 4)
        TIHXInSecondary = readValues(times, loglines, "faceZone secondaryIHXInlet TBulk", 4)
        TIHXOutSecondary = readValues(times, loglines, "faceZone secondaryIHXOutlet TBulk", 4)
        THotLegPrimary = readValues(times, loglines, "faceZone hotLegMiddleCut TBulk", 4)
        TColdLegPrimary = readValues(times, loglines, "faceZone coldLegMiddleCut TBulk", 4)


        axT.plot(times, TCoreIn, label="coreInlet", ls=ls)
        axT.plot(times, TCoreOut, label="coreOutlet", ls=ls)
        axT.plot(times, TInnerCoreOut, label="innerCoreOutlet", ls=ls)
        axT.plot(times, TOuterCoreOut, label="outerCoreOutlet", ls=ls)
        axT.plot(times, THotLegPrimary, label="THotLegPrimary", ls=ls)
        axT.plot(times, TColdLegPrimary, label="TColdLegPrimary", ls=ls)
        #axT.plot(times, TIHXInPrimary, label="IHXPrimaryInlet", ls=ls)
        #axT.plot(times, TIHXOutPrimary, label="IHXPrimaryOutlet", ls=ls)
        #axT.plot(times, TIHXInSecondary, label="IHXSecondaryInlet", ls=ls)
        #axT.plot(times, TIHXOutSecondary, label="IHXSecondaryOutlet", ls=ls)

        try:
            totalPower = readValues(times, loglines, "totalPower =", 2, 180)
            fissionPower = readValues(times, loglines, "-> fission", 3, 180)
            decayPower = readValues(times, loglines, "-> decay", 3, 180)
            axP.plot(times, totalPower, label="Total Power", ls=ls)
            axP.plot(times, fissionPower, label="Fission Power", ls=ls)
            axP.plot(times, decayPower, label="Decay Power", ls=ls)
        except:
            pass

        try:
            RTot = readValues(times, loglines, "totalReactivity", 2, 1.0/beta)
            RDoppler = readValues(times, loglines, "Doppler (fast)", 4, 1.0/beta)
            RTFuel = readValues(times, loglines, "-> TFuel", 3, 1.0/beta)
            RTClad = readValues(times, loglines, "-> TClad", 3, 1.0/beta)
            RRhoCool = readValues(times, loglines, "-> rhoCool", 3, 1.0/beta)
            RTStruct = readValues(times, loglines, "-> TStruct ", 3, 1.0/beta)
            RDriveline = readValues(times, loglines, "-> driveline", 3, 1.0/beta)
            RGEM = readValues(times, loglines, "-> GEM", 3, 1.0/beta)

            axR.plot(times, RTot, label="Total", ls=ls)
            axR.plot(times, RDoppler, label="Doppler", ls=ls)
            axR.plot(times, RTFuel, label="Fuel axial expansion", ls=ls)
            axR.plot(times, RTClad, label="Cladding temperature", ls=ls)
            axR.plot(times, RRhoCool, label="Coolant density", ls=ls)
            axR.plot(times, RTStruct, label="Core radial expansion", ls=ls)
            axR.plot(times, RDriveline, label="Driveline", ls=ls)
            axR.plot(times, RGEM, label="GEM", ls=ls)
        except:
            pass


#==============================================================================*

timeOffset = 910

with open("expFlow", "r") as log:
    loglines = log.readlines()
    expFlowTimes, expFlowValues = readExpTimesAndValues(loglines, timeOffset, 0, 1)
    axF.plot(expFlowTimes, expFlowValues, label="Exp. Total", ls='--', marker='x', markevery=5)

with open("expPIOTA2", "r") as log:
    loglines = log.readlines()
    expPIOTA2Times, expPIOTA2Values = readExpTimesAndValues(loglines, timeOffset, 273.15, 1)
    axT.plot(expPIOTA2Times, expPIOTA2Values, label="Exp. PIOTA2", ls='--', marker='x', markevery=5)

'''
with open("expPIOTA6", "r") as log:
    loglines = log.readlines()
    expPIOTA6Times, expPIOTA6Values = readExpTimesAndValues(loglines, timeOffset, 273.15, 1)
    axT.plot(expPIOTA6Times, expPIOTA6Values, label="Exp. PIOTA6", ls='--', marker='x', markevery=5)
'''

with open("expPow", "r") as log:
    loglines = log.readlines()
    expPowTimes, expPowValues = readExpTimesAndValues(loglines, timeOffset, 0, 1000000)
    axP.plot(expPowTimes, expPowValues, label="Exp. Primary", ls='--', marker='x', markevery=5)

with open("expReactivity", "r") as log:
    loglines = log.readlines()
    expRTimes, expRValues = readExpTimesAndValues(loglines, timeOffset, 0, 1)
    axR.plot(expRTimes, expRValues, label="Exp. Total", ls='--', marker='x', markevery=5)


#==============================================================================*

for ax in axes.flatten():
    ax.legend(loc='upper right')

if (times[0] >= 900):
    axR.set_xlim(900, 1200)

fig.tight_layout()
fig.savefig("FFTF_transient.png")
plt.show()


#==============================================================================*
