"""
Script to plot the residuals during an OpenFOAM simulation.

Author: Thomas Guilbaud, EPFL/Transmutex SA, 12/09/2023
"""
#==============================================================================*

# Import
import matplotlib.pyplot as plt
import os
import sys
import numpy as np
from scipy.optimize import curve_fit

#==============================================================================*

if (len(sys.argv) != 2):
    print(f"\n    Usage: python3 {sys.argv[0]} path/to/log.GeN-Foam\n")
    sys.exit(0)

filename = sys.argv[1]

if (not os.path.isfile(filename)):
    print(f"\nFile path '{filename}' does not exist !\n")
    sys.exit(1)

#==============================================================================*

class Residuals:
    def __init__(self, parameters: list, regionName: str="") -> None:
        self.parameters: dict = {}
        for param in parameters:
            self.parameters[param] = []

        self.filename: str = ""
        self.regionName: str = regionName


    def extractFromFile(self, filename: str) -> bool:
        """
        Extraction of the residuals in the log file.

        return True if at least one data has been extracted.
        """

        self.filename = filename

        with open(filename, 'r') as file:
            print(f"Process {self.filename} in region {self.regionName}")

            for line in file:
                for param in self.parameters:
                    if (f"Solving for {param}" in line):
                        self.parameters[param].append(
                            float(line.split()[11].split(",")[0])
                        )
                    elif (f"Residual {param}" in line):
                        self.parameters[param].append(
                            float(line.split()[3].split(",")[0])
                        )
                    elif ("Final neutroResidual" in line and param == "neutroResidual"):
                        self.parameters[param].append(
                            float(line.split()[6])
                        )
                    elif (param == "keff" and "keff_ = " in line):
                        self.parameters[param].append(float(line.split()[2]))
                    elif (param == "ExecutionTime" and "ExecutionTime = " in line):
                        self.parameters[param].append(float(line.split()[2]))
                    elif (param == "T (avg min max)" and "T (avg min max)" in line):
                        self.parameters[param].append(float(line.split()[7]))

        return(len(self.parameters[param]) > 0)


    def plot(self, ax, paramToFit: str=None, axZoom=None, lastPoint: float=-1) -> None:
        """ Plot all the residuals """
        print("Plot", list(self.parameters.keys()))

        plotTitle = self.filename.split(".")[-1] + " " + self.regionName
        ax.set_title(plotTitle)
        for param in self.parameters:
            lenParam = len(self.parameters[param])
            if (lenParam > 0 and param != "ExecutionTime"):
                ax.plot(self.parameters[param], label=param)
                if (axZoom != None and lastPoint > 0):
                    axZoom.set_title(plotTitle+f"\nlast {lastPoint} points")
                    axZoom.plot(self.parameters[param][-lastPoint:], label=param)
                    axZoom.legend()
                ax.legend()

            if (param == paramToFit and lenParam > 0):
                print()
                lastKeffs = self.parameters[param][-10:]
                if (param == "keff"):
                    lastDiffRho = [(1/ki-1/kf) * 1e5 for ki, kf in zip(lastKeffs[:-1], lastKeffs[1:])]
                else:
                    lastDiffRho = [kf-ki for ki, kf in zip(lastKeffs[:-1], lastKeffs[1:])]
                avgReactivityVar = np.mean(lastDiffRho)
                print(f"d{param}/diter = {avgReactivityVar:.4f} pcm/iteration")

                nFit: int = lastPoint # 100
                if (lenParam > nFit):
                    def func(x, a, b, d):
                        return d + a * np.exp(-b * x)

                    xdata = np.array(list(range(len(self.parameters[param]))))
                    currentIteration = max(xdata)
                    lastKeff = self.parameters[param][-1]

                    try:
                        popt, pcov = curve_fit(
                            func, xdata[-nFit:], self.parameters[param][-nFit:],
                            p0=[0, 0, lastKeff],
                            bounds=([-1000, -1000, 0.9], [1000, 1000, 4000])
                        )
                        finalKeff = func(currentIteration + 1e9, *popt)
                        diffKeff = (1/self.parameters[param][-1] - 1/finalKeff) * 1e5

                        a, b, d = popt
                        minThreshold: float = 1e-5
                        iterationToThreshold =  1/b * np.log(abs(a) / minThreshold)
                        remainingIteration = iterationToThreshold - currentIteration
                    except:
                        print("Fiting an exponential didn't succeed.")

                    try:
                        ax.plot(xdata[-nFit:], func(xdata[-nFit:], *popt), '--', label='fit: $%5.3f e^{-%5.3f x} + %5.5f$' % tuple(popt))
                        axZoom.plot(func(xdata[-lastPoint:], *popt), '--', label='fit: $%5.3f e^{-%5.3f x} + %5.5f$' % tuple(popt))
                        ax.hlines(xmin=0, xmax=max(xdata), y=finalKeff, ls='--', color="tab:grey", alpha=0.5)
                        if (abs(finalKeff - lastKeff) < abs(lastKeff - self.parameters[param][-nFit])):
                            axZoom.hlines(xmin=0, xmax=lastPoint, y=finalKeff, ls='--', color="tab:grey", alpha=0.5)
                        ax.legend()

                        currentTime = max(self.parameters['ExecutionTime'])
                        timeToThreshold = iterationToThreshold / currentIteration * currentTime
                        remainingTime = timeToThreshold - currentTime

                        print(f"Min iteration to be < {minThreshold*1e5:.2f} pcm")

                        print(f"{'Param':10} | {'Iteration':10} | {'Time [s]':15} | {'keff':12}")
                        print(f"{'':->10}-|-{'':->10}-|-{'':->15}-|-{'':->12}")
                        print(f"{'Estimated':10} | {iterationToThreshold:10.0f} | {f'{timeToThreshold:.1f} = {timeToThreshold/3600:.1f}h':>15} | {finalKeff:12.5f}")
                        print(f"{'Current':10} | {currentIteration:10.0f} | {f'{currentTime:.1f} = {currentTime/3600:.1f}h':>15} | {lastKeff:12.5f}")
                        print(f"{'Remain':10} | {remainingIteration:10.0f} | {f'{remainingTime:.1f} = {remainingTime/3600:.1f}h':>15} | {f'{diffKeff:.2f} pcm':>12}")
                    except:
                        pass

#==============================================================================*

# Residuals to plot
ResidualsNeutronics = Residuals(
    parameters=["fluxStar0", "fluxStar1", "fluxStar20", "fluxStar21", "angularFlux_0_0", "angularFlux_1_0", "neutroResidual"],
    regionName="Neutronics"
)
ResidualsFluid = Residuals(
    parameters=["p_rgh", "h", "k", "epsilon"],
    regionName="Thermal-hydraulics"
)
ResidualsKeff = Residuals(
    parameters=["keff", "ExecutionTime"], #, "T (avg min max)"],
    regionName=r"$k_{eff}$"
)

# Extract from file
isResNeutro = ResidualsNeutronics.extractFromFile(filename=filename)
isResFluid  = ResidualsFluid.extractFromFile(filename=filename)
isResKeff  = ResidualsKeff.extractFromFile(filename=filename)

#==============================================================================*

# Create figure
fig, axes = plt.subplots(ncols=3, nrows=2, figsize=(14, 8))
axNeutronics, axFluid, axKeff, axNeutronicsZoom, axFluidZoom, axKeffZoom = axes.flatten()

# Plot the resisuals
ResidualsNeutronics.plot(ax=axNeutronics, axZoom=axNeutronicsZoom, lastPoint=200)
ResidualsFluid.plot(ax=axFluid, axZoom=axFluidZoom, lastPoint=200)
# ResidualsKeff.plot(ax=axKeff, axZoom=axKeffZoom, lastPoint=2000, paramToFit="T (avg min max)")
ResidualsKeff.plot(ax=axKeff, axZoom=axKeffZoom, lastPoint=100, paramToFit="keff")

# Options
for ax in axes.flatten():
    ax.set_xlabel("Iterations")
    ax.set_ylabel("Residual")
    ax.set_yscale('log')
    ax.grid()

axKeff.set_ylabel(r"$k_{eff}$")
axKeffZoom.set_ylabel(r"$k_{eff}$")
axKeff.set_yscale('linear')
axKeffZoom.set_yscale('linear')

fig.tight_layout()

# Save the figure
fig.savefig(f"fig_results_residuals_{filename.split('/')[0]}.png")

print("End")


#==============================================================================*
