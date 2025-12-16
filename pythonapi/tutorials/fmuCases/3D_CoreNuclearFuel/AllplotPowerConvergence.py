"""
Script to plot the power evolution in GeN-Foam coupled with fuel pin FMU.

Author: Thomas Guilbaud, EPFL/Transmutex SA, 2023
"""

#=============================================================================*
# Imports

import matplotlib.pyplot as plt
import sys

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels

#=============================================================================*
# Usage

if (len(sys.argv) == 1):
    print(f"\n    Usage: python3 {sys.argv[0]} <ext> <path1/to/log.GeN-Foam> ...\n")
    sys.exit(0)


#=============================================================================*
# Functions

def relErr(th: list, sim: list) -> list:
    return([abs(thi-simi)/thi * 100 for thi, simi in zip(th, sim)])


def extractDataFromLogfile(filename: str):
    try:
        with open(filename, 'r') as file:
            text = f"| {filename.split('/')[0]:40} "
            tRange = []
            dtRange = []
            totalPowerHeatFluxTemp, totalPowerRhoCpdtTemp, totalPowerNeutronicsTemp = 0, 0, 0
            totalPowerHeatFlux, totalPowerRhoCpdt, totalPowerNeutronics = [], [], []
            timeTemp, deltaTtemp, clocktime = 0, 0, 0
            isAdd = False

            # Extract data line by line
            for line in file:
                if ("Time = " == line[:7]):
                    if (timeTemp > 0):
                        tRange.append(timeTemp)
                        dtRange.append(deltaTtemp)
                        totalPowerHeatFlux.append(totalPowerHeatFluxTemp)
                        totalPowerRhoCpdt.append(totalPowerRhoCpdtTemp)
                        totalPowerNeutronics.append(totalPowerNeutronicsTemp)
                    timeTemp = float(line.split()[-1])
                elif ("step_size = " in line):
                    deltaTtemp = float(line.split(";")[2].split()[2])
                elif (line[:3] == "End" and len(line) <= 4):
                    tRange.append(timeTemp)
                    totalPowerHeatFlux.append(totalPowerHeatFluxTemp)
                    totalPowerRhoCpdt.append(totalPowerRhoCpdtTemp)
                    totalPowerNeutronics.append(totalPowerNeutronicsTemp)
                elif ("Integrated power in cellZone innerCore:" in line):
                    isAdd = False
                elif ("Integrated power in cellZone outerCore:" in line):
                    isAdd = True
                elif ("heat flux  = " in line):
                    totalPowerHeatFluxTemp = float(line.split()[3]) / 1e6 + (totalPowerHeatFluxTemp if isAdd else 0)
                elif ("rhoCpdTdt  = " in line):
                    totalPowerRhoCpdtTemp = float(line.split()[2]) / 1e6 + (totalPowerRhoCpdtTemp if isAdd else 0)
                elif ("neutronics = " in line):
                    totalPowerNeutronicsTemp = float(line.split()[2]) / 1e6 + (totalPowerNeutronicsTemp if isAdd else 0)
                elif ("ClockTime" in line):
                    clocktime = float(line.split()[6])

            # Plot data
            if (len(totalPowerHeatFlux) > 0):
                text += f"| {tRange[-1]:8.1f} | {clocktime:8.1f} | {clocktime/3600:8.2f} | {clocktime/3600 * 200/tRange[-1]:8.2f} | {totalPowerHeatFlux[-1]:.3f}, {totalPowerRhoCpdt[-1]:.3f}, {totalPowerNeutronics[-1]:.3f}"

                return(text, tRange, totalPowerHeatFlux, totalPowerRhoCpdt, totalPowerNeutronics)

            else:
                return(f"No power found in {filename}", [], [], [], [])

    except:
        return(f"No file named {filename} or issues with this file", [], [], [], [])

    return("", [], [], [], [])


def plot(
    ax,
    tRange: list[float],
    values: list[float],
    axError=None,
    funcError=lambda p: p,
    ls: str="",
    label: str=""
) -> None:

    ax.plot(tRange, values, ls=ls, label=label, markersize=4)

    if (axError != None):
        axError.plot(tRange, [funcError(p) for p in values], ls=ls, label=label, markersize=4)


def addPlotOptions(
    ax,
    tBT: float=0, dtBT: float=0,
    xmin: float=0, xmax: float=0, ymin: float=0, ymax: float=0,
    targetPower: float=0,
    label: str=""
) -> None:
    if (tBT > 0 and dtBT > 0 and ymax != ymin):
        ax.fill_between(x=[tBT, tBT+dtBT], y1=ymin, y2=ymax, color='tab:grey', alpha=0.3)

    if (targetPower != None and xmin != xmax):
        ax.hlines(xmin=xmin, xmax=xmax, y=[targetPower], color="tab:grey", alpha=0.5, label=label)


def createZoom(ax, zoomLoc, xlim, ylim):
    axZoom = ax.inset_axes(
        zoomLoc,
        xlim=xlim, ylim=ylim,
        # xticklabels=[], yticklabels=[]
    )
    _, c1 = ax.indicate_inset_zoom(axZoom, edgecolor="black")

    c1[0].set_visible(False)
    c1[1].set_visible(False)
    c1[2].set_visible(True)
    c1[3].set_visible(True)

    return(axZoom)


def extractBT(string: str) -> float:
    try:
        return(float(string.split("BT")[0].split("_")[-1]))
    except:
        return(0)


#=============================================================================*

ext = ""
filenames = sys.argv[1:]
if ("log" not in sys.argv[1]):
    filenames = sys.argv[2:]
    ext = f"_{sys.argv[1]}"


# fuelFraction: float = 0.363516 # 0.21733
targetPower: float = 300 # MW
tBT: float = 100 # s, put large number to avoid plotting the beam trip band
dtBT: float = extractBT(filenames[0]) # s

for filename in filenames:
    if (dtBT != extractBT(filename)):
        dtBT = 0
        break

linestyles = ["-", "--", "-.", ":", (5, (10, 3)), (0, (5, 3))]
markers = ["", "o", "*", "x"]
labels = [" ".join(name.split("/")[0].split("_")[-2:]) for name in filenames]
if (dtBT == 0):
    labels = [f"BT {extractBT(name)} s" for name in filenames]


fig, axes = plt.subplots(nrows=3, figsize=(6, 8))
(axHeatFlux, axRhoCpdTdt, axNeutronics) = axes.flatten()

figError, axesError = plt.subplots(nrows=3, figsize=(6, 8))
(axHeatFluxError, axRhoCpdTdtError, axNeutronicsError) = axesError.flatten()

# Add zoom window
zoomLoc = [0.7, 0.15, 0.28, 0.7]
zoomXlim=(tBT - 2, tBT + (dtBT+5 if dtBT > 0 else 22))

maxHeatFlux = 330

axZoomHeatFlux = createZoom(axHeatFlux, zoomLoc, xlim=zoomXlim, ylim=(0, maxHeatFlux))
axZoomRhoCpdTdt = createZoom(axRhoCpdTdt, zoomLoc, xlim=zoomXlim, ylim=(-230, 230))
axZoomNeutronics = createZoom(axNeutronics, zoomLoc, xlim=zoomXlim, ylim=(0, 310))


# Read each file
isAboveBeamTrip = False
print()
print(f"| {'Filename':40} | {'Last':^8} | {'Clock':^8} | {'Clock':^8} | {'Total':^8} |")
print(f"| {'':40} | {'time [s]':^8} | {'Time [s]':^8} | {'Time [h]':^8} | {'time [h]':^8} |")
print(f"|:{'':-^40}-|:{'':-^8}:|:{'':-^8}:|:{'':-^8}:|:{'':-^8}:|")
for i, label, filename in zip(range(len(filenames)), labels, filenames):
    ls = linestyles[i % len(linestyles)]

    text, tRange, totalPowerHeatFlux, totalPowerRhoCpdt, totalPowerNeutronics = extractDataFromLogfile(filename)

    # if (i == 1):
    #     totalPowerHeatFlux = [val * 3.635160e-01 for val in totalPowerHeatFlux]
    #     totalPowerRhoCpdt = [val * 3.635160e-01 for val in totalPowerRhoCpdt]
    #     totalPowerNeutronics = [val * 3.635160e-01 for val in totalPowerNeutronics]

    if (len(tRange) > 0 and tRange[-1] > tBT):
        isAboveBeamTrip = True

    print(text)

    maxHeatFlux = max([maxHeatFlux] + totalPowerHeatFlux)

    plot(axHeatFlux, tRange, totalPowerHeatFlux, ls=ls, label=label, axError=axHeatFluxError, funcError=lambda p: abs(1-p/targetPower)*100)
    plot(axRhoCpdTdt, tRange, totalPowerRhoCpdt, ls=ls, label=label, axError=axRhoCpdTdtError, funcError=lambda p: abs(p/targetPower)*100)
    plot(axNeutronics, tRange, totalPowerNeutronics, ls=ls, label=label, axError=axNeutronicsError, funcError=lambda p: abs(1-p/targetPower)*100)

    plot(axZoomHeatFlux, tRange, totalPowerHeatFlux, ls=ls, label=label)
    plot(axZoomRhoCpdTdt, tRange, totalPowerRhoCpdt, ls=ls, label=label)
    plot(axZoomNeutronics, tRange, totalPowerNeutronics, ls=ls, label=label)



# Add options to the plot
for ax in [axHeatFlux, axZoomHeatFlux]:
    addPlotOptions(ax, tBT=tBT, dtBT=dtBT, xmin=0, xmax=200, ymin=0, ymax=maxHeatFlux, targetPower=targetPower, label="Equilibrium power")
for ax in [axRhoCpdTdt, axZoomRhoCpdTdt]:
    addPlotOptions(ax, tBT=tBT, dtBT=dtBT, xmin=0, xmax=200, ymin=-250, ymax=250, targetPower=0)
for ax in [axNeutronics, axZoomNeutronics]:
    addPlotOptions(ax, tBT=tBT, dtBT=dtBT, xmin=0, xmax=200, ymin=0, ymax=maxHeatFlux, targetPower=targetPower)

for ax in [axHeatFluxError, axRhoCpdTdtError, axNeutronicsError]:
    addPlotOptions(ax, tBT=tBT, dtBT=dtBT, xmin=0, xmax=200, ymin=5e-3, ymax=100)


axHeatFlux.set_ylim((0, maxHeatFlux if maxHeatFlux < 1000 else 330))
axRhoCpdTdt.set_ylim((-300, 300))
axNeutronics.set_ylim((0, 310))

# Vertical axis label
axHeatFlux.set_ylabel("Heat flux [MW]")
axRhoCpdTdt.set_ylabel(r"Enthalpy $\left( \rho c_p \frac{dT}{dt} \right)$ [MW]")
axNeutronics.set_ylabel("Neutronics [MW]")

axHeatFluxError.set_ylabel("Relative difference heat flux [%]")
axRhoCpdTdtError.set_ylabel(r"Relative difference $\rho c_p \frac{dT}{dt}$ [%]")
axNeutronicsError.set_ylabel("Relative difference neutronics [%]")

axHeatFlux.legend(bbox_to_anchor=(0, 1.08), loc="lower left", ncols=3)
axHeatFluxError.legend(bbox_to_anchor=(0, 1.08), loc="lower left", ncols=3)

# Common options
for ax in axes.flatten():
    ax.set_xlabel("Time [s]")
    ax.set_xlim((0, 250))
    # ax.set_ylabel("Intergrated power [MW]")
    # ax.xaxis.set_ticks(np.arange(0, max(tRange)+0.1, 0.1))
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    # ax.legend(bbox_to_anchor=(0, 1.08), loc="lower left", ncols=3)

for ax in axesError.flatten():
    ax.set_xlabel("Time [s]")
    ax.set_yscale('log')
    ax.set_xlim((0, 200))
    ax.set_ylim((5e-3, 100))

# Save
fig.tight_layout()
fig.savefig(f"fig_results_totalPower{ext}.png")
# fig.savefig(f"fig_results_totalPower{ext}.pdf")
figError.tight_layout()
figError.savefig(f"fig_results_totalPowerError{ext}.png")
# figError.savefig(f"fig_results_totalPowerError{ext}.pdf")

print()

# plt.show()

#=============================================================================*
