### IMPORTS

import sys
import matplotlib.pyplot as plt

### MAIN

if (len(sys.argv) < 2):
    print("Usage: python plotMassFLowRate.py log.GeN-Foam")
    sys.exit(1)


for filename in sys.argv[1:]:
    with open(filename, "r") as file:
        times, mFlowRates = [], []
        time = 0
        
        for line in file.readlines():
            if "Time =" in line and not "ExecutionTime" in line :
                time = float(line.split()[2])
            if "faceZone massFlowSurface_z1 massFlow" in line :
                mFlowRate = float(line.split()[4])
                
                mFlowRates.append(mFlowRate)
                times.append(time)
            
        t0 = times[0]
        for i in range(len(times)):
            times[i] = times[i] - t0

        del times[0]
        del mFlowRates[0]
 
 

plt.figure()
plt.plot(times, mFlowRates)
plt.xlabel("iteration")
plt.ylabel("Mass flow rate")
plt.grid()
plt.show()
