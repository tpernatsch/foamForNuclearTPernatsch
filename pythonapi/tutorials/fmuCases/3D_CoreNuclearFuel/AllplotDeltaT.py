"""
>>> python3 AllplotDeltaT.py steadyState/Core/log.GeN-Foam steadyState/FuelPin0/log.offbeat
"""

#=============================================================================*

import matplotlib.pyplot as plt
import sys
import numpy as np

#=============================================================================*

if (len(sys.argv) == 1):
    print(f"\n    Usage: python3 {sys.argv[0]} <path1/to/log.GeN-Foam> ...\n")
    sys.exit(0)

#=============================================================================*

def relErr(th: list, sim: list) -> list:
    return([abs(thi-simi)/thi * 100 for thi, simi in zip(th, sim)])


#=============================================================================*

filenames = sys.argv[1:]
linestyles = ["-", "--", "-.", ":", (5, (10, 3)), (0, (5, 10))]
markers = ["o", "", "*", "x"]

minDeltaT: float = 0.05

fig, ax = plt.subplots(dpi=200)

# Read each file
for i, filename in enumerate(filenames):
    with open(filename, 'r') as file:
        filenameNoPath: str = filename.split('/')[-1]
        print(f"Analyzing {filename} ... ")
        tRange = []
        deltaTtemp, evaluatedDeltaTtemp = 0, 0
        deltaT, evaluatedDeltaT = [], []

        # Extract data line by line
        for line in file:
            if ("Time = " == line[:7]):
                tRange.append(float(line.split()[2]))
                if (len(tRange) > 1):
                    deltaT.append(deltaTtemp)
                    evaluatedDeltaT.append(evaluatedDeltaTtemp)
            elif ("adjustTimeStep" in line):
                evaluatedDeltaTtemp = float(line.split()[21])
            elif ("deltaT = " in line):
                deltaTtemp = float(line.split()[2])

        # Plot data
        if (len(evaluatedDeltaT) > 0):
            if (len(tRange) != len(evaluatedDeltaT)):
                tRange = tRange[:-1]

            if (i == 0):
                ax.plot(tRange, deltaT, label=r"Used $\Delta t$", marker='o')

            ax.plot(tRange, evaluatedDeltaT, label=fr"Evaluated $\Delta t$ {filenameNoPath}", marker='*')

        else:
            print("No value found in", filename)


# Add options to the plot
if (max(tRange) > 100):
    ax.fill_between(x=[100, 110], y1=minDeltaT/5, y2=1e2, color='tab:grey', alpha=0.3)
ax.hlines(y=minDeltaT, xmin=0, xmax=max(tRange), color='tab:grey', alpha=0.5, ls='--')
ax.text(s=fr'$\Delta t_{{min}}$ = {minDeltaT}', x=max(tRange), y=minDeltaT, ha='right', va='top')
ax.legend()
ax.set_xlabel("Time [s]")
ax.set_ylabel(r"Step size $\Delta t$ [s]")
ax.set_yscale('log')
# ax.set_ylim(minDeltaT/5, 1e2)

# Save
fig.tight_layout()
fig.savefig("fig_results_deltaT.png")

# plt.show()

#=============================================================================*
