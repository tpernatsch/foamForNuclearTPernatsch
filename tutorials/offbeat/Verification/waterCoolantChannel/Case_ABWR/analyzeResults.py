# ----------------------------------------------------------------------------------------------- #
# Python script for the comparison of OFFBEAT coolant channel verification results
# ----------------------------------------------------------------------------------------------- #

import numpy as np
import re
import matplotlib.pyplot as plt
from scipy.interpolate import RBFInterpolator, InterpolatedUnivariateSpline

# ----------------------------------------------------------------------------------------------- #

def fileReader(filePath, firstLine, numberOfColumns):
    myData = []

    with open(filePath, 'r') as myFile:
        myLines = myFile.readlines()
        
        pattern = re.compile(r'\s+')
        
        for myLine in myLines[firstLine - 1:]:
            splitLine = pattern.split(myLine.strip())
            if len(splitLine) == numberOfColumns:
                floatLine = list(map(float, splitLine))
                myData.append(floatLine)
                
    dataArray = np.array(myData)
    return dataArray

# Read TRACE file
filePathTRACE = 'TRACE_ABWR_Twall.txt'
firstLineTRACE = 4
numberOfColumnsTRACE = 12
dataTRACE = fileReader(filePathTRACE, firstLineTRACE, numberOfColumnsTRACE)

# Read OFFBEAT file
filePathOFFBEAT = 'postProcessing/patchProbes_cladOuter/0/T'
firstLineOFFBEAT = 13
numberOfColumnsOFFBEAT = 12
dataOFFBEAT = fileReader(filePathOFFBEAT, firstLineOFFBEAT, numberOfColumnsOFFBEAT)

# Resample both to 1 Hz
newTimeValues = np.arange(0.0, 28800, 1.0)
newDataTRACE = newTimeValues
newDataOFFBEAT = newTimeValues
for myColumns in range(1, numberOfColumnsTRACE):
    ius = InterpolatedUnivariateSpline(dataTRACE[:, 0], dataTRACE[:, myColumns])
    newDataTRACE = np.column_stack((newDataTRACE, ius(newTimeValues)))
    ius = InterpolatedUnivariateSpline(dataOFFBEAT[:, 0], dataOFFBEAT[:, myColumns])
    newDataOFFBEAT = np.column_stack((newDataOFFBEAT, ius(newTimeValues)))
    
# ----------------------------------------------------------------------------------------------- #
# Plots
# ----------------------------------------------------------------------------------------------- #

# Get a colormap and extract colors
cmap = plt.get_cmap("tab10")  
colors = [cmap(i / 10) for i in range(10)]

cm=1/2.54
fig, axs = plt.subplots(2, 1, figsize=(15*cm, 15*cm), sharex=True)

nNodes = np.array([5])
# nNodes = np.arange(1,10)

maxTempDifference = -10

for i, plottedNode in enumerate(nNodes):
    plottedNode = int(plottedNode)

    # Plot 1: TRACE vs. OFFBEAT wall temperature for "plottedNode"
    axs[0].plot(newDataOFFBEAT[:, 0], newDataOFFBEAT[:, plottedNode], ls='--', color=colors[i], label=f'node {plottedNode}')
    axs[0].plot(newDataTRACE[:, 0], newDataTRACE[:, plottedNode], ls='-', color=colors[i])

    # Plot 2: TRACE - OFFBEAT wall temperature difference for "plottedNode"
    myDifference = newDataTRACE[:, plottedNode] - newDataOFFBEAT[:, plottedNode]
    axs[1].plot(newDataOFFBEAT[:, 0], myDifference, color='b', label='Difference')

    maxTempDifference = max(max(myDifference), maxTempDifference)

axs[0].set_ylabel('Wall temperature (K)')
axs[0].legend(frameon=False, bbox_to_anchor=(1.1,0.5), loc='center left')
axs[0].grid(ls=":", alpha=0.3)
axs[0].set_title("OFFBEAT: dashed lines, TRACE: solid lines")

axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel('Wall temperature difference (K)')
axs[1].grid(ls=":", alpha=0.3)

plt.tight_layout()
plt.savefig('results')

# ----------------------------------------------------------------------------------------------- #
# Analysis
# ----------------------------------------------------------------------------------------------- #

with open("results", "w") as f:
    f.write(f"Maximum difference OFFBEAT vs TRACE: {maxTempDifference:.2f} K")