"""
Script to analyse the system response to a new power request.

Author: Thomas Guilbaud, EPFL/Transmutex SA, 04/2023
"""

#=============================================================================*

import sys
import matplotlib.pyplot as plt
import pandas as pd

#=============================================================================*

if (len(sys.argv) < 2):
    print("Usage: python plot.py folder/path ...")
    sys.exit(1)

#=============================================================================*

def relError(l1: list, l2: list) -> list:
    return(
        list(map(lambda e1, e2: abs(e1-e2)/e1, l1, l2))
    )

def listModifier(liste: list, scale: float=1, offset: float=0) -> list:
    return(
        [scale*value+offset for value in liste]
    )

#=============================================================================*

class SList:
    """ Super list to manipulate the values from the log. """

    def __init__(
        self,
        keyWord: str,
        ax=None,
        pos: int=2,
        exception: str="",
        linestyle: str="-",
        label: str=""
    ) -> None:
        self.list = []
        self.currentValue = 0
        self.keyWord = keyWord
        self.exception = exception
        self.pos = pos

        # Plotting
        self.ax = ax
        self.linestyle = linestyle
        self.label = label

    def getFromLine(self, line: str) -> None:
        self.currentValue = float(line.split()[self.pos])

    def testAndGetFromLine(self, line: str) -> None:
        if (
            self.keyWord in line and (
                self.exception not in line or self.exception == ""
            )
        ):
            self.getFromLine(line)

    def addValue(self) -> None:
        self.list.append(self.currentValue)

    def modifyList(self, scale: float=1, offset: float=0) -> None:
        self.list = listModifier(self.list, scale=scale, offset=offset)

    def plot(self, time: list) -> None:
        self.ax.plot(
            time[1:],
            self.list[:-1],
            linestyle=self.linestyle,
            label=self.label
        )


#=============================================================================*

linestyles = ['-', '--', '-.']

fig, axes = plt.subplots(3, sharex=True)

(axPower, axReactivity, axTemperature) = axes.flatten()

folders = sys.argv[1:]

for folder, linestyle in zip(folders, linestyles):
    records = {
        "time": SList("Time = ", exception="ExecutionTime"),
        "powers": SList("totalPower", ax=axPower, linestyle=linestyle, label="Point-Kinetics"),
        "TFuels": SList("TFuel = ", ax=axTemperature, linestyle=linestyle, label="Fuel"),
        "totRhos": SList("totalReactivity", ax=axReactivity, linestyle=linestyle, label="Total"),
        "fuelRhos": SList("-> TFuel", pos=3, ax=axReactivity, linestyle=linestyle, label="TFuel"),
        "extRhos": SList("-> extReactivity", pos=3, ax=axReactivity, linestyle=linestyle, label="Ext react")
    }

    # Read file
    with open(f"{folder}/log.GeN-Foam", "r") as file:
        # Read lines
        for line in file.readlines():
            for record in records.values():
                record.testAndGetFromLine(line)

            # When end of time step, add values to the list
            if ("ClockTime" in line):
                for record in records.values():
                    record.addValue()

    # Offset the time
    times = records["time"]
    times.modifyList(offset=-100)
    times = times.list

    # Plot
    records["powers"].plot(times)
    records["totRhos"].plot(times)
    records["extRhos"].plot(times)
    records["fuelRhos"].plot(times)
    records["TFuels"].plot(times)


    # Read from FMU output
    data = pd.read_csv(f"{folder}/ExternalReactivityController.csv")

    time = listModifier(data['time'], offset=-100)
    idxTime = time.index(0)

    powerGF = list(data['power'])
    powerRamp = list(data['ramp.y'])
    diffPower = relError(powerRamp, powerGF)
    extRhoFMUs = listModifier(data['pid.y'], scale=1e5)

    # print(
    #     f"At time t=0 (idx {idxTime}): {time[idxTime:idxTime+4]}\n",
    #     f"rhoTot={records['totRhos'].list[0:4]}\n",
    #     f"PK            = {records['powers'].list[0:4]}\n",
    #     f"fieldIntegral = {powerGF[idxTime:idxTime+4]}\n",
    #     f"PID output    = {list(data['pid.y'])[idxTime:idxTime+4]}\n",
    #     f"ramp          = {powerRamp[idxTime:idxTime+4]}\n",
    #     f"Rel error power = {diffPower[idxTime:idxTime+4]}\n"
    # )

    # axPower.plot(time, powerGF, label="fieldIntegralToFMU")
    axPower.plot(time, powerRamp, label='Command')


axPower.set_xlim(0)
# axPower.set_ylim((9.8e6, 11.3e6))

# axReactivity.plot(time, extRhoFMUs, label="PID output")

#=============================================================================*

# Plot options
for ax in axes.flatten():
    ax.set_xlabel("Time [s]")
    ax.grid(True)
    ax.legend()

axPower.set_ylabel("Power [W]")
axReactivity.set_ylabel("Reactivity [pcm]")
axTemperature.set_ylabel("Fuel temperature [K]")

fig.tight_layout()
fig.savefig(f"fig_results_{'_'.join([folder.replace('/', '') for folder in folders])}.png")

# plt.show()

#=============================================================================*