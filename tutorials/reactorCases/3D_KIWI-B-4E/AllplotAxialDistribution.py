"""
Script to plot axial distributions from multiple files. Plot
- Coolant and fuel temperatures
- Coolant density
- Neutron fluxes (group 0 and 1)
- Power density
- Pressure
- Velocity

Thomas Guilbaud, EPFL, 13/08/2023
"""

#==============================================================================*
# Imports
#==============================================================================*

import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


#==============================================================================*
# Constants
#==============================================================================*

SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 15

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)   # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)   # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


#==============================================================================*
# Reactor Parameters
#==============================================================================*

dzCore: float = 1.32 # m
fuelVolumeFraction: float = 0.672179


#==============================================================================*
# Usage
#==============================================================================*

if (len(sys.argv) <= 1):
    print(f"\n  Usage: python3 {sys.argv[0]} path/to/data_<>.xy path2/to/data_<>.xy ...\n")
    sys.exit(1)


#==============================================================================*
# Usefull Functions and Classes
#==============================================================================*

def listModifier(l: list, scale: float=1, offset: float=0) -> list:
    return([scale*e+offset for e in l])


def doubleArrow(
        ax, x, y, dx, head_width: float=50, head_length: float=0.05, color: str='black'
    ):
    ax.arrow(
        x=x, y=y, dx=dx, dy=0, length_includes_head=True,
        head_width=head_width, head_length=head_length, color=color
    )
    ax.arrow(
        x=dx, y=y, dx=-dx, dy=0, length_includes_head=True,
        head_width=head_width, head_length=head_length, color=color
    )

def plotRefPressure(ax, color: str='tab:purple'):
    """
    Place experimental measures of pressure in MPa from LA-3185-MS page 66
    at t = 20800 s
    """
    return(
        ax.scatter(
            [-0.1, 1.4],
            [3.951, 3.372], # [3.861, 3.427],
            color=color,
            marker='x',
            linewidths=2
        )
    )

def plotRefTemperatureCoolant(ax, color: str='tab:orange'):
    """
    Place experimental measures of coolant temperature in K from LA-3185-MS
    page 66 at t = 20800 s
    """
    return(
        ax.scatter(
            [-0.1, 1.4],
            [93.9, 1955.6], # [83.3, 1972],
            color=color,
            marker='x',
            linewidths=2
        )
    )

def plotRefTemperatureFuel20400s(ax, color: str='tab:red'):
    """ Place experimental measures of fuel temperature in K from LA-3185-MS
    page  at t = 20400 s"""
    return(
        ax.scatter(
            [
                0.007184827, 0.051138366, 0.203163399, 0.353516277,
                0.509448103, 0.659833269, 0.815814375, 0.966354182,
                1.120617251
            ],
            [
                197.3703022, 359.9080283, 898.9766668, 1238.636455,
                1496.085303, 1796.686747, 1994.520225, 2108.052751,
                2161.981344
            ],
            color=color,
            marker='x',
            label="Experimental measures\n"+r"of T$_{fuel}$ at t=20400s",
            linewidths=2
        )
    )


def plotRefTemperatureFuel20800s77inch(ax, color: str='tab:red'):
    """ Place experimental measures of fuel temperature in K from LA-3185-MS
    page 174 at t = 20800 s and 7.7 inches from center"""
    return(
        ax.scatter(
            [0.0508, 0.2032, 0.3556, 0.508, 0.6604, 0.762, 0.9652, 1.1176],
            [394.8191313268603, 931.4015583216104, 1232.9704966133552,
             1507.672676042831, 1910.3547118491426, 2002.7632591135157,
             2166.9073451177037, 1956.6004291150812
            ],
            color=color,
            marker='x',
            label="Experimental measures\n"+r"of T$_{fuel}$ at t=20800s",
            linewidths=2
        )
    )

def plotRefTemperatureFuel20800s97inch(ax, color: str='tab:red'):
    """ Place experimental measures of fuel temperature in K from LA-3185-MS
    page 174 at t = 20800 s and 9.7 inches from center"""
    return(
        ax.scatter(
            [
                0.011305961, 0.049234335, 0.20089857, 0.354596042, 0.506047548,
                0.659738302, 0.813321572, 0.960283388, 1.120349279
            ],
            [
                211.5019132, 355.0787416, 873.998247, 1178.923642, 1458.668523,
                1756.041035, 1932.567421, 2164.468899, 2129.527279
            ],
            color=color,
            marker='x',
            label="Experimental measures\n"+r"of T$_{fuel}$ at t=20800s",
            linewidths=2
        )
    )

def plotRefPowerDensity(ax, color: str='tab:grey'):
    """ From 1965 LANL - LA-3323-MS, Table I Thermal and thermoelatic input
    parameters """
    inch = 2.54e-2 # m
    Btuinch3s = 1055.056 / (inch**3) / 1e9
    return(
        ax.scatter(
            [val * inch for val in [4, 12, 20, 24, 28, 36, 48]], # in m
            [6.721789e-01 * val * Btuinch3s for val in [22.0, 37.1, 43.7, 44.8, 43.6, 37.0, 15.4]], # in kW/cm3
            color=color,
            marker='x',
            label="Experimental measures\n"+r"of Power density at t=20800s",
            linewidths=2
        )
    )


class Figure:
    """
    Figure object to simplify the creation of graphs.
    """
    def __init__(self, title: str, yLabel: str, isYlimZero: bool=False, isYlog: bool=False) -> None:
        self.fig, self.ax = plt.subplots(figsize=(6, 4.5), dpi=300)
        # self.fig.subplots_adjust(left=0.15, right=0.85)
        # self.ax.set_title(title)
        self.ax.set_xlim((0, dzCore))
        self.ax.set_xlabel("Reactor axial position [m]")
        self.ax.set_ylabel(yLabel)
        # self.ax.grid()
        self.isYlimZero = isYlimZero

        if (isYlog):
            self.ax.set_yscale('log')

        self.allYList = []

    def plot(self, x: list, y: list, label: str, ls: str, color: str=None) -> None:
        self.allYList += y
        self.ax.plot(x, y, label=label, ls=ls, color=color)

    def finalAndSave(self, filename: str) -> None:
        # Core limits
        if (len(self.allYList) > 0):
            minValue, maxValue = min(self.allYList), max(self.allYList)
            avgValue = (maxValue+minValue)/2
            # self.ax.vlines(
            #     x=[0, dzCore],
            #     ymin=0 if self.isYlimZero else minValue,
            #     ymax=maxValue,
            #     ls='--', color='grey'
            # )
            self.ax.vlines(
                x=[dzCore/2], ymin=0, ymax=maxValue,
                color='grey', alpha=0.5, zorder=-3
            )

            # Hydrogen flow arrow
            self.ax.text(
                x=dzCore/2,
                y=minValue+0.25*(maxValue-minValue)/2,
                s=r'H$_2$ Flow',
                ha='center'
            )
            self.ax.arrow(
                x=dzCore/3,
                y=0.995*(minValue+0.23*(maxValue-minValue)/2),
                dx=dzCore/3,
                dy=0,
                length_includes_head=True,
                width=0.003*avgValue,
                head_length=dzCore/30,
                head_width=0.020*avgValue,
                color='black'
            )

            self.ax.legend()

            if (self.isYlimZero):
                self.ax.set_ylim(0)
            else:
                self.ax.set_ylim((0.9*minValue, 1.1*maxValue))

        self.fig.tight_layout()
        self.fig.savefig(filename)


#==============================================================================*
# Extract Results and Plot
#==============================================================================*

# Create figures and axis
figPower = Figure(
    "Core axial power density distribution",
    r"Power density [kW/cm$^3$]",
    isYlimZero=True
)
figFlux = Figure(
    "Core axial flux distribution",
    r"Neutron flux [n/cm$^2$/s]",
    isYlog=True,
    isYlimZero=False
)
figTemperature = Figure(
    "Core axial temperature distribution",
    r"Temperature [K]",
    isYlimZero=True
)
figDensity = Figure(
    "Core axial hydrogen density distribution",
    r"Hydrogen density [kg/m$^3$]",
    isYlog=True
)
figVelocity = Figure(
    "Core axial hydrogen velocity magnitude distribution",
    r"Velocity magnitude [m/s]",
    isYlimZero=True
)
figPressure = Figure(
    "Core axial pressure distribution",
    r"Pressure [MPa]"
)

# --- Create figure summary setup
figSummary, axPressure = plt.subplots(1, figsize=(8, 5), dpi=200)
figSummary.subplots_adjust(right=0.6, top=0.78) # right=0.75, top=0.8)
axPressure.xaxis.set_ticks(np.arange(-10, 10, 0.2))
axPressure.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
# Create new axis for power, temperature and velocity quantities
axPowerDensity = axPressure.twinx()
axTemperature = axPressure.twinx()
axVelocity = axPressure.twinx()
axCoolantDensity = axPressure.twinx()
# Move new axis to avoid overlap
axVelocity.spines.right.set_position(("axes", 1.25))
axPowerDensity.spines.right.set_position(("axes", 1.5))
axCoolantDensity.spines.right.set_position(("axes", 1.75))

for label in (
    axPressure.get_xticklabels() +
    axPowerDensity.get_yticklabels() +
    axTemperature.get_yticklabels() +
    axVelocity.get_yticklabels() +
    axPressure.get_yticklabels() +
    axCoolantDensity.get_yticklabels()
):
    label.set_fontsize(MEDIUM_SIZE)


linestyles = ['-', '--', '-.', ':']

for filename, ls in zip(sys.argv[1:], linestyles):
    # File name
    folder = filename.split('/')[0]
    timeStep = filename.split('/')[-2]
    label = filename.split('/')[-1].split('_')[0]
    print(f"Process {folder} / {label} at time {timeStep}...")

    # Update label if multiple files provided
    if (len(sys.argv[1:]) >= 2):
        label += " t="+timeStep

    # Extract all the parameter names from the filename (suppose OpenFOAM
    # autogeneration)
    names = ['z']+filename.split('/')[-1].split('.xy')[0].split('_')[1:]

    # Create empty list for each parameters
    data = {name: [] for name in names}

    # Open the file
    with open(filename, 'r') as file:
        # Extract all the data line by line
        for line in file.readlines():
            # Extract and format
            line = [float(e) for e in line.split()]

            # Fill in lists
            for i, name in enumerate(names):
                data[name].append(line[i])

        # Inverse Z-axis to have the values increasing with Z
        zRange = [abs(dzCore/2-z) for z in data['z']]

        # Rescale variables
        for d in data:
            if ("flux" in d):
                # Convert from n/m2/s to n/cm2/s
                data[d] = listModifier(data[d], scale=1e-4)
                figFlux.plot(zRange, data[d], label=label+" "+d, ls=ls)

            elif ("power" in d):
                # Convert from W to MW
                data[d] = listModifier(data[d], scale=fuelVolumeFraction*1e-9)
                figPower.plot(zRange, data[d], label=label, ls=ls)

            elif ("T" == d[0]):
                figTemperature.plot(zRange, data[d], label=label+" "+d.split(".")[0], ls=ls)

            elif ("rho" == d[:3]):
                figDensity.plot(zRange, data[d], label=label, ls=ls)

            elif ("magU" == d):
                figVelocity.plot(zRange, data[d], label=label, ls=ls)

            elif ("p" == d):
                # Convert from Pa to MPa
                data[d] = listModifier(data[d], scale=1e-6)
                figPressure.plot(zRange, data[d], label=label, ls=ls)

    if ("neutroRegion" in filename):
        pRho, = axCoolantDensity.plot(zRange, data['rhoCool'], '-.', color='tab:cyan', label='H$_2$ density')
        axCoolantDensity.set_ylim(0)

    # Plot summary figure only if fluidRegion
    if ("fluidRegion" in filename):
        minP, maxP = min(data['p']), max(data['p'])
        minT, maxT = min(data['T']), max(data['T'])

        # Reference values
        # plotRefPressure(axPressure, color="tab:purple")
        # plotRefTemperatureCoolant(axTemperature, color="tab:blue")
        if ("R4" in filename):
            sTmatrix = plotRefTemperatureFuel20800s77inch(axTemperature, color="tab:red")

        elif ("R5" in filename):
            sTmatrix = plotRefTemperatureFuel20800s97inch(axTemperature, color="tab:red")

        elif ("R1" in filename):
            sPowerDensity = plotRefPowerDensity(axPowerDensity, color="tab:grey")

            axTemperature.scatter(
                [val * 2.54e-2 for val in [4, 12, 20, 24, 28, 36, 48]],
                [val * 0.55555 for val in [415, 1074, 1986, 2459, 2916, 3714, 4437]],
                color="tab:blue",
                marker='x',
                label=r"$T_{bulk}$"
            )
            sTmatrix = axTemperature.scatter(
                [val * 2.54e-2 for val in [4, 12, 20, 24, 28, 36, 48]],
                [val * 0.55555 for val in [1000, 1880, 2780, 3240, 3630, 4270, 4630]],
                color="tab:red",
                marker='x',
                label=r"$T_{fuel}$"
            )

        # Plot
        pPower, = axPowerDensity.plot(zRange, data['powerDensityStructure'], color="tab:grey", label="Power Density", ls=(0, (3, 1, 1, 1, 1, 1)))
        # axPowerDensity.fill_between(
        #     x= zRange,
        #     y1= data['powerDensityStructure'],
        #     color= "tab:gray",
        #     alpha= 0.1
        # )
        pTbulk, = axTemperature.plot(zRange, data['T'], '-', color='tab:blue', label=r'T$_{bulk}$')
        pTsurf, = axTemperature.plot(zRange, data['Tsurface.lumpedNuclearStructure'], '--', color='tab:orange', label=r'T$_{surface}$')
        pTmatrix, = axTemperature.plot(zRange, data['Tmatrix.lumpedNuclearStructure'], ':', color='tab:red', label=r'T$_{matrix}$')
        pP, = axPressure.plot(zRange, data['p'], '-.', color='tab:purple', label='Pressure')
        pU, = axVelocity.plot(zRange, data['magU'], '-.', color='tab:green', label='Velocity')

        ## Plot details
        # axTemperature.vlines(
        #     x=[0, dzCore/2, dzCore], ymin=0, ymax=2*maxT,
        #     color='grey', linestyles='dashed'
        # )
        axTemperature.vlines(
            x=[dzCore/2], ymin=0, ymax=2*maxT,
            color='grey', alpha=0.5, zorder=-3
        )
        # axTemperature.text(x=dzCore/2, y=minT, s='Core', ha='center')
        # doubleArrow(axTemperature, 0, 0.9*minT, dzCore)
        axTemperature.text(x=dzCore/2, y=minT+0.25*(maxT-minT)/2, s=r'H$_2$ Flow', ha='center')
        axTemperature.arrow(
            x=dzCore/3, y=0.9*(minT+0.25*(maxT-minT)/2), dx=dzCore/3, dy=0,
            length_includes_head=True,
            head_width=50,
            head_length=0.05,
            color='black'
        )

        # Axis titles
        axPressure.set_xlabel('Core axial position [m]', fontsize=MEDIUM_SIZE)
        axPressure.set_ylabel('Pressure [MPa]', fontsize=MEDIUM_SIZE)
        axPressure.set_xlim(0.000001, dzCore)
        axPressure.set_ylim(0.99*minP, 1.01*maxP)
        axTemperature.set_ylabel('Temperature [K]', fontsize=MEDIUM_SIZE)
        axTemperature.set_ylim((0, 2300))
        # axTemperature.set_ylim((0, 2700))
        axVelocity.set_ylabel('Velocity [m/s]', fontsize=MEDIUM_SIZE)
        axVelocity.set_ylim(0, 700)
        axPowerDensity.set_ylabel('Power density [kW/cm$^3$]', fontsize=MEDIUM_SIZE)
        axPowerDensity.set_ylim(0)
        axCoolantDensity.set_ylabel('Coolant density [kg/m$^3$]', fontsize=MEDIUM_SIZE)


# Plot reference values
# plotRefPressure(figPressure.ax, color="tab:blue")
# plotRefTemperatureCoolant(figTemperature.ax, color="tab:blue")
suffix = "R1"
if ("R4" in filename):
    plotRefTemperatureFuel20800s77inch(figTemperature.ax, color="black")
    suffix = "R4"
elif ("R5" in filename):
    plotRefTemperatureFuel20800s97inch(figTemperature.ax, color="black")
    suffix = "R5"


# Add vertical lines, legend and save
figPower.finalAndSave(f"fig_results_axialPowerDistribution_{suffix}.png")
figFlux.finalAndSave(f"fig_results_axialFluxDistribution_{suffix}.png")
figTemperature.finalAndSave(f"fig_results_axialTemperatureDistribution_{suffix}.png")
figDensity.finalAndSave(f"fig_results_axialDensityDistribution_{suffix}.png")
figVelocity.finalAndSave(f"fig_results_axialVelocityDistribution_{suffix}.png")
figPressure.finalAndSave(f"fig_results_axialPressureDistribution_{suffix}.png")

# Legend and grid
axPressure.legend(
    handles=[pTbulk, pTsurf, pTmatrix, pP, pU, pPower, pRho, sTmatrix],
    bbox_to_anchor=(0, 1.01, 1, 0.2),
    loc="lower left",
    mode='expand',
    ncol=2
)

figSummary.tight_layout()
# Save figure
figSummary.savefig(f'fig_results_axialDistributionSummary_{suffix}.png')


#==============================================================================*
