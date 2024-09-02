"""
Author: Thomas Guilbaud, EPFL/Transmutex, 12/07/2023
"""

import matplotlib.pyplot as plt


def readFromFile(filename: str) -> dict: 
    """
    Read the result file from an FMU.
    filename: Name of the result file (i.e. fmu1_1_result.txt)
    return a 
    """
    dataDict = {}

    with open(filename, "r") as file:
        isReadHeader, isReadData = False, False
        for line in file.readlines():
            if ("char name(" in line):
                isReadHeader = True
            elif ("char description(" in line):
                isReadHeader = False
            elif ("float data_2(" in line):
                isReadData = True

            elif (isReadHeader and line != "\n"):
                colName = line.split("\n")[0]
                dataDict[colName] = []

            elif (isReadData):
                data = [float(d) for d in line.split()]
                for key, d in zip(dataDict, data):
                    dataDict[key].append(d)

    return(dataDict)


dataReactor = readFromFile("Reactor_0_result.txt")

plt.figure("Steam generator power")
steamGenPower = [0]*len(dataReactor["time"])
for i in range(8):
    plt.plot(
        dataReactor["time"], 
        [power/1e3 for power in dataReactor[f"gfSteamGenPowerOut{i+1}"]], 
        label=f"SG{i+1}"
    )
    steamGenPower = [pi+ps for pi, ps in zip(steamGenPower, dataReactor[f"gfSteamGenPowerOut{i+1}"])]
steamGenPower = [abs(p * 360/5 * 56.69862) / 1e6 for p in steamGenPower]
plt.xlabel("Time [s]")
plt.ylabel("Power extracted from the steam generator [MW]")
plt.legend()
plt.grid(True)
plt.savefig("results_powerSteamGen.png")


plt.figure("Steam generator temperature")
for i in range(8):
    plt.plot(
        dataReactor["time"], 
        dataReactor[f"gfSteamGenTemperature{i+1}"], 
        label=f"SG{i+1}"
    )
plt.xlabel("Time [s]")
plt.ylabel("Steam/water temperature in steam generator [K]")
plt.legend()
plt.grid(True)
plt.savefig("results_temperatureSteamGen.png")


plt.figure("Power")
corePower = [(pi+po)*360/5 * 198.28145 * 1e-6 for pi, po in zip(
    dataReactor["gfInnerCorePowerOut"], dataReactor["gfOuterCorePowerOut"]
)]
plt.plot(dataReactor["time"], corePower, label=f"Core")
plt.plot(dataReactor["time"], steamGenPower, label=f"SG")
plt.xlabel("Time [s]")
plt.ylabel("Power [MW]")
plt.legend()
plt.grid(True)
plt.savefig("results_power.png")


plt.figure("Reactivity")
plt.plot(
    dataReactor["time"], 
    [rho*1e5 for rho in dataReactor["gfExternalReactivity"]], 
    label=f"Ext FMU"
)
plt.xlabel("Time [s]")
plt.ylabel("Reactivity [pcm]")
plt.legend()
plt.grid(True)
plt.savefig("results_reactivity.png")
