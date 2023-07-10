import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys


if (len(sys.argv) < 2):
    print("Usage: python plot.py log.GeN-Foam")
    sys.exit(1)


for filename in sys.argv[1:]:
    with open(filename, "r") as file:
        timesLog, mFlowRates = [], []
        time = 0
        
        for line in file.readlines():
            if "Time =" in line and not "ExecutionTime" in line:
                time = float(line.split()[2])
            if "faceZone massFlowSurface_z1 massFlow" in line:
                mFlowRate = float(line.split()[4])
                
                mFlowRates.append(mFlowRate)
                timesLog.append(time)
            
        for i in range(len(timesLog)):
            timesLog[i] = timesLog[i] - timesLog[0]

        del mFlowRates[0]
        del timesLog[0]


mFlowRates = np.array(mFlowRates)
mFlowRatesScaled = (mFlowRates - mFlowRates[0])/mFlowRates[0]


data = pd.read_csv("momentumSourceTest.csv")

time = data["time"]
momentumSource = data["model.root.system1.momentumModelica"]
# desiredMomentumSource = data['model.root.system1.ramp.y']
power = data["model.root.system1.power"]
# desiredPower = data["model.root.system1.pi.SP"]

# desiredPowerScaled = data["model.root.system1.pi.SPs"]
# powerScaled = data["model.root.system1.pi.PVs"]

# momentumSourceScaled = data["model.root.system1.pi.CSs"]


# plt.plot(time, powerScaled, label='normalized power variation')
# plt.plot(time, desiredPowerScaled, '--', label='normalized desired power')
#plt.plot(time, momentumSourceScaled, label='normalized momentumSource variation')
plt.plot(timesLog, mFlowRatesScaled, label='normalized momentumSource variation')
# plt.plot(time, desiredMomentumSource, label='normalized momentumSource variation')

plt.xlabel("Time [s]")
plt.ylabel("Normalized variation")
plt.legend()

plt.show()
