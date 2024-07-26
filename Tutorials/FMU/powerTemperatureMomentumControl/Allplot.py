import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

# Usage
if (len(sys.argv) != 2):
    print(f"Usage: python3 {sys.argv[0]} log.GeN-Foam")
    sys.exit(1)

filename = sys.argv[1]


def normalizedList(l: list) -> list:
    first = l[100]
    return(
        [(e-first)/first for e in l]
    )


volume = 4 # m3
fractionStructure = 0.6
volumetricArea = 1000 # m2/m3

data = pd.read_csv("momentumSourceTest.csv")

timeModelica = data["time"]
modMomentumSourceOut1 = data["model.root.system1.modMomentumSourceOut1"]
modMomentumSourceOut2 = data["model.root.system1.modMomentumSourceOut2"]
modPowerOut = data["model.root.system1.modPowerOut"]
modTemperatureOut = data["model.root.system1.modTemperatureOut"]
modTin = data["model.root.system1.modTin"]
modPowerIn = data["model.root.system1.modPowerIn"]


fig, (axMomentum, axPower, axTemperature) = plt.subplots(nrows=3, figsize=(8, 6))


timesLog, mFlowRateOutlet1s, mFlowRateOutlet2s, TBulkOutlet1s, TBulkOutlet2s = [], [], [], [], []

with open(filename, "r") as file:
    print(f"Analyze {filename} ...")

    for line in file.readlines():
        if "Time =" in line and not "ExecutionTime" in line:
            time = float(line.split()[2])
        if "patch outlet1 TBulk" in line:
            TBuldOutlet1 = float(line.split()[4])
        if "patch outlet2 TBulk" in line:
            TBulkOutlet2 = float(line.split()[4])
        if "patch outlet1 massFlow" in line:
            mFlowRateOutlet1 = float(line.split()[4])
        if "patch outlet2 massFlow" in line:
            mFlowRateOutlet2 = float(line.split()[4])
            
            mFlowRateOutlet1s.append(mFlowRateOutlet1)
            mFlowRateOutlet2s.append(mFlowRateOutlet2)
            TBulkOutlet1s.append(TBuldOutlet1)
            TBulkOutlet2s.append(TBulkOutlet2)
            timesLog.append(time)
        
    for i in range(len(timesLog)):
        timesLog[i] = timesLog[i] - timesLog[0]

    del mFlowRateOutlet1s[0]
    del mFlowRateOutlet2s[0]
    del TBulkOutlet1s[0]
    del TBulkOutlet2s[0]
    del timesLog[0]



axMomentum.plot(timeModelica, normalizedList(modMomentumSourceOut1), label='Modelica 1')
axMomentum.plot(timeModelica, normalizedList(modMomentumSourceOut2), label='Modelica 2')
axMomentum.plot(timesLog, normalizedList(mFlowRateOutlet1s),  ls='--', label='GeN-Foam 1')
axMomentum.plot(timesLog, normalizedList(mFlowRateOutlet2s),  ls='--', label='GeN-Foam 2')

axMomentum.set_xlabel("Time [s]")
axMomentum.set_ylabel("Normalized variation\nmomentum source")
axMomentum.legend()


axPower.plot(timeModelica, [p * volume for p in modPowerOut], label="Modelica Out")
axPower.plot(timeModelica, [p * fractionStructure/volumetricArea for p in modPowerIn], label="Modelica In (integrated heat flux)")

axPower.set_xlabel("Time [s]")
axPower.set_ylabel("Power [W]")
axPower.legend()

axTemperature.plot(timeModelica, modTemperatureOut, label="Modelica Tout")
axTemperature.plot(timeModelica, modTin, label="Modelica Tin")
axTemperature.plot(timesLog, TBulkOutlet1s, ls='--', label="GeN-Foam 1")
axTemperature.plot(timesLog, TBulkOutlet2s, ls='--', label="GeN-Foam 2")

axTemperature.set_xlabel("Time [s]")
axTemperature.set_ylabel("Temperature [K]")
axTemperature.set_ylim(290)
axTemperature.legend()

fig.tight_layout()

plt.show()
