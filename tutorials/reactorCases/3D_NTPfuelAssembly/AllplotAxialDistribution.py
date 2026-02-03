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
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
        

#==============================================================================*
# Reactor Parameters
#==============================================================================*

dzCore = 1.32 # m


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


def doubleArrow(ax, x, y, dx, head_width, color):
    ax.arrow(
        x=x, y=y, dx=dx, dy=0, length_includes_head=True, 
        head_width=head_width, color=color
    )
    ax.arrow(
        x=dx, y=y, dx=-dx, dy=0, length_includes_head=True, 
        head_width=head_width, color=color
    )

def plotRefPressure(ax, color: str='tab:purple'):
    """ Place experimental measures of pressure in MPa """
    return(
        ax.scatter(
            [-0.1, 1.4], 
            [3.861, 3.427],
            color=color,
            marker='x',
            linewidths=2
        )
    )

def plotRefTemperatureCoolant(ax, color: str='tab:orange'):
    """ Place experimental measures of coolant temperature in K """
    return(
        ax.scatter(
            [-0.1, 1.4], 
            [83.3, 1972],
            color=color,
            marker='x',
            linewidths=2
        )
    )

def plotRefTemperatureFuel(ax, color: str='tab:red'):
    """ Place experimental measures of fuel temperature in K """
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


class Figure:
    """
    Figure object to simplify the creation of graphs.
    """
    def __init__(self, title: str, yLabel: str, isYlimZero: bool=False, isYlog: bool=False) -> None:
        self.fig, self.ax = plt.subplots(figsize=(7, 5))
        self.fig.subplots_adjust(left=0.15, right=0.85)
        self.ax.set_title(title)
        self.ax.set_xlabel("Reactor axial position [m]")
        self.ax.set_ylabel(yLabel)
        self.ax.grid()
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
            self.ax.vlines(
                x=[0, dzCore], 
                ymin=minValue, 
                ymax=maxValue, 
                ls='--', color='grey'
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

        self.fig.savefig(filename)


#==============================================================================*
# Extract Results and Plot
#==============================================================================*

# Create figures and axis
figPower = Figure(
    "Core axial power density distribution",
    r"Power density [MW/m$^3$]",
    isYlimZero=True
)
figFlux = Figure(
    "Core axial flux distribution",
    r"Neutron flux [n/cm$^2$/s]",
    isYlog=True
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
figSummary, axPressure = plt.subplots(1, figsize=(8, 5))
figSummary.subplots_adjust(right=0.65, top=0.78) # right=0.75, top=0.8)
axPressure.xaxis.set_ticks(np.arange(-10, 10, 0.2))
axPressure.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
# Create new axis for power, temperature and velocity quantities
axPowerDensity = axPressure.twinx()
axTemperature = axPressure.twinx()
axVelocity = axPressure.twinx()
# Move new axis to avoid overlap
axVelocity.spines.right.set_position(("axes", 1.2))
axPowerDensity.spines.right.set_position(("axes", 1.4))


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
        zRange = [dzCore-z for z in data['z']]

        # Rescale variables
        for d in data:
            if ("flux" in d):
                # Convert from n/m2/s to n/cm2/s
                data[d] = listModifier(data[d], scale=1e-4)
                figFlux.plot(zRange, data[d], label=label+" "+d, ls=ls)
                
            elif ("power" in d):
                # Convert from W to MW
                data[d] = listModifier(data[d], scale=1e-6)
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

    # Plot summary figure only if fluidRegion
    if ("fluidRegion" in filename):
        minP, maxP = min(data['p']), max(data['p'])

        for label in (
            axPressure.get_xticklabels() +
            axPowerDensity.get_yticklabels() +
            axTemperature.get_yticklabels() +
            axVelocity.get_yticklabels() +
            axPressure.get_yticklabels()
        ):
	        label.set_fontsize(MEDIUM_SIZE)

        # Reference values
        plotRefPressure(axPressure)
        plotRefTemperatureCoolant(axTemperature)
        sTmatrix = plotRefTemperatureFuel(axTemperature)

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
        axPressure.vlines(
            x=[0, dzCore], ymin=0.9*minP, ymax=1.1*maxP, 
            color='grey', linestyles='dashed'
        )
        axPressure.text(x=dzCore/2, y=minP, s='Core', ha='center')
        doubleArrow(axPressure, 0, 0.995*minP, dzCore, 0.025, 'black')
        axPressure.text(x=dzCore/2, y=minP+0.25*(maxP-minP)/2, s=r'H$_2$ Flow', ha='center')
        axPressure.arrow(
            x=dzCore/3, y=0.995*(minP+0.25*(maxP-minP)/2), dx=dzCore/3, dy=0, 
            length_includes_head=True, head_width=0.025, color='black'
        )

        # Axis titles
        axPressure.set_xlabel('Reactor axial position [m]', fontsize=MEDIUM_SIZE)
        axPressure.set_ylabel('Pressure [MPa]', fontsize=MEDIUM_SIZE)
        axPressure.set_ylim(0.99*minP, 1.01*maxP)
        axTemperature.set_ylabel('Temperature [K]', fontsize=MEDIUM_SIZE)
        axTemperature.set_ylim(0)
        axVelocity.set_ylabel('Velocity [m/s]', fontsize=MEDIUM_SIZE)
        axVelocity.set_ylim(0, 700)
        axPowerDensity.set_ylabel('Power density [MW/m$^3$]', fontsize=MEDIUM_SIZE)
        axPowerDensity.set_ylim(0)

        # Legend and grid
        axPressure.legend(
            handles=[pTbulk, pTsurf, pTmatrix, pP, pU, pPower, sTmatrix], 
            bbox_to_anchor=(0, 1.01, 1, 0.2), 
            loc="lower left", 
            mode='expand', 
            ncol=2
        )
        # axTemperature.grid(True)

        # Save figure
        figSummary.savefig('fuelElementAxialDistribution.png')


# Plot reference values
plotRefPressure(figPressure.ax, color="tab:blue")
plotRefTemperatureCoolant(figTemperature.ax, color="tab:blue")
plotRefTemperatureFuel(figTemperature.ax, color="black")

# Add vertical lines, legend and save
figPower.finalAndSave("axialPowerDistribution.png")
figFlux.finalAndSave("axialFluxDistribution.png")
figTemperature.finalAndSave("axialTemperatureDistribution.png")
figDensity.finalAndSave("axialDensityDistribution.png")
figVelocity.finalAndSave("axialVelocityDistribution.png")
figPressure.finalAndSave("axialPressureDistribution.png")

plt.show()

#==============================================================================*
