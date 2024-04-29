import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

# Usage
if (len(sys.argv) < 2):
    print(f"Usage: python3 {sys.argv[0]} log.GeN-Foam")
    sys.exit(1)

filnames = sys.argv[1:]


def normalizedList(l: list) -> list:
    first = l[0]
    return(
        [(e-first)/first for e in l]
    )

data = pd.read_csv("momentumSourceTest.csv")

timeModelica = data["time"]
momentumSourceModelica = data["model.root.system1.momentumModelica"]
powerModelica = data["model.root.system1.power"]

# desiredMomentumSource = data['model.root.system1.ramp.y']
# desiredPower = data["model.root.system1.pi.SP"]
# desiredPowerScaled = data["model.root.system1.pi.SPs"]
# powerScaled = data["model.root.system1.pi.PVs"]
# momentumSourceScaled = data["model.root.system1.pi.CSs"]



# plt.plot(time, powerScaled, label='normalized power variation')
# plt.plot(time, desiredPowerScaled, '--', label='normalized desired power')
fig, (axMomentum, axPower, axTemperature) = plt.subplots(nrows=3, figsize=(8, 6))

linestyles = ['-', '--', '-.', ':']

for filename, ls in zip(filnames, linestyles):
    timesLog, mFlowRates, TBulkCoreInlets, TBulkCoreOutlets, totalPowers = [], [], [], [], []

    with open(filename, "r") as file:
        print(f"Analyze {filename} ...")

        for line in file.readlines():
            if "Time =" in line and not "ExecutionTime" in line:
                time = float(line.split()[2])
            if "totalPower = " in line:
                totalPower = float(line.split()[2])
            if "faceZone TcoreInlet TBulk" in line:
                TBulkCoreInlet = float(line.split()[4])
            if "faceZone TcoreOutlet TBulk" in line:
                TBulkCoreOutlet = float(line.split()[4])
            if "faceZone massFlowSurface_z1 massFlow" in line:
                mFlowRate = float(line.split()[4])
                
                mFlowRates.append(mFlowRate)
                TBulkCoreInlets.append(TBulkCoreInlet)
                TBulkCoreOutlets.append(TBulkCoreOutlet)
                totalPowers.append(totalPower)
                timesLog.append(time)
            
        for i in range(len(timesLog)):
            timesLog[i] = timesLog[i] - timesLog[0]

        del totalPowers[0]
        del mFlowRates[0]
        del TBulkCoreInlets[0]
        del TBulkCoreOutlets[0]
        del timesLog[0]


    mFlowRatesScaled = normalizedList(mFlowRates)
    totalPowersScaled = normalizedList(totalPowers)


    # axMomentum.plot(timeModelica, momentumSourceModelica, label='Modelica')
    axMomentum.plot(timesLog, mFlowRatesScaled,  ls=ls, label='GeN-Foam')
    # plt.plot(time, desiredMomentumSource, label='normalized momentumSource variation')

    axMomentum.set_xlabel("Time [s]")
    axMomentum.set_ylabel("Normalized variation\nmomentum source")
    axMomentum.legend()

    # axPower.plot(timeModelica, powerModelica, label="Modelica")
    axPower.plot(timesLog, totalPowersScaled, ls=ls, label="GeN-Foam")
    
    axPower.set_xlabel("Time [s]")
    axPower.set_ylabel("Normalized variation\nmomentum source")
    axPower.legend()

    axTemperature.plot(timesLog, TBulkCoreInlets, ls=ls, label="GeN-Foam inlet")
    axTemperature.plot(timesLog, TBulkCoreOutlets, ls=ls, label="GeN-Foam outlet")
    
    axTemperature.set_xlabel("Time [s]")
    axTemperature.set_ylabel("Temperature [K]")
    axTemperature.legend()

fig.tight_layout()

plt.show()
