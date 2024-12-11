"""
Extract and plot the neutron flux from GeN-Foam and compare to analytic
solution. Extract the keff from uniform/reactorState

Author: Thomas Guilbaud
"""

#==============================================================================*
# Imports

import sys
import numpy as np
from fluidfoam import readmesh, readscalar
import matplotlib.pyplot as plt
from scipy import optimize

#==============================================================================*
# Functions

def extractFieldOverLine(region: str, timeStep: str, fieldName: str):
    # Read the mesh
    x, y, z = readmesh(".", region=region, structured=False)

    # Extract the field
    field = readscalar(".", timeStep, fieldName, region=region, structured=False)

    return(z, field)


def getKeff(filename: str) -> float:
    with open(filename, 'r') as f:
        for line in f.readlines():
            if ("keff" in line):
                return(float(line.split()[-1].split(';')[0]))

    return(0)


def analyticSolution(fuelLength, reflWidth):
    DR = 0.2
    DC = 0.1275
    SigmaAR = 3
    LR = np.sqrt(DR/SigmaAR)

    func = lambda B: np.tan(B*fuelLength/2) - DR/(DC*LR*B) / np.tanh(reflWidth/LR)

    sol = optimize.root(func, [1])
    B = sol.x[0]

    def phi(z):
        if (np.abs(z) <= fuelLength/2):
            return(np.cos(B*z))
        return(np.cos(B*fuelLength/2) * np.sinh((fuelLength/2 + reflWidth - np.abs(z))/LR)) / np.sinh(reflWidth/LR)

    return(phi)


def plotFlux(
        timeStep: str,
        fieldName: str
    ) -> None:
    """
    isAddExtrapolationDistance:
        Turn to True to see the extrapolated flux (see "Possible exercises"
        section)
    """

    # Extract the flux
    z, flux0 = extractFieldOverLine(
        region="neutroRegion",
        timeStep=timeStep,
        fieldName=fieldName
    )

    fuelLength = 1.197
    reflWidth = 0.5
    fuelLengthNoRefl = 1.5

    # Normalize the flux
    normFlux = [flux0_ / max(flux0) for flux0_ in flux0]

    # Analytic solution
    fluxTh = analyticSolution(fuelLength, reflWidth)
    fluxThNoRefl = lambda z_: np.cos(z_ * np.pi/fuelLengthNoRefl) if np.abs(z_) < fuelLengthNoRefl/2 else 0

    # Compute the relative error between GeN-Foam and the analytic solution
    relError = [(normFlux_ / fluxTh(z_) - 1) * 100 for z_, normFlux_ in zip(z, normFlux)]

    # Plot
    fig, axes = plt.subplots(2, 1, figsize=(6, 6))
    axFlux, axError = axes.flatten()

    ymax = max(normFlux) * 1.05

    axFlux.plot(z, normFlux, label="GeN-Foam")
    axFlux.plot(z, [fluxTh(z_) for z_ in z], label="Analytic", ls="--")
    axFlux.plot(z, [fluxThNoRefl(z_) for z_ in z], label="Analytic no refl", ls="-.")

    # Core boundaries
    axFlux.vlines(
        x=[-fuelLength/2, fuelLength/2], ymin=0, ymax=ymax,
        color="tab:grey", alpha=0.5
    )
    axFlux.vlines(
        x=[-fuelLengthNoRefl/2, fuelLengthNoRefl/2], ymin=0, ymax=ymax,
        color="tab:grey", alpha=0.5, ls='-.'
    )
    axFlux.set_ylabel("Neutron flux [a.u]")
    axFlux.set_ylim((0, ymax))
    axFlux.legend()

    axError.plot(z, relError)
    axError.set_ylabel("Relative error [%]")

    for ax in axes:
        ax.set_xlabel("Axial position [m]")
        ax.set_xlim(min(z), max(z))

    axError.hlines(xmin=min(z), xmax=max(z), y=[0], alpha=0.5, color="tab:grey")

    fig.tight_layout()

    fig.savefig("results_neutronFlux.png")

    plt.show()


#==============================================================================*
# Main

if __name__ == "__main__":
    # Usage
    if (len(sys.argv) == 1):
        print(f"\n    Usage: python3 {sys.argv[0]} timeStep fieldName\n")
        sys.exit(0)

    # Extract user parameters
    timeStep = sys.argv[1]
    fieldName = sys.argv[2]
    isTest = False
    if (len(sys.argv) >= 4):
        isTest = bool(sys.argv[3])

    # Plot and extract keff
    if (not isTest):
        plotFlux(timeStep, fieldName)

        # Extract keff
        keff = getKeff(f"{timeStep}/uniform/reactorState")

        print(f"\nAt t={timeStep}, keff = {keff}, expected keff = 1.0\n")


    # Test only keff
    if (isTest):
        expectedKeff = 1.0
        error = 1e-5

        keff = getKeff(f"{timeStep}/uniform/reactorState")

        print(f"keff with {error} relative error")

        print(f"|      | Simulated | Expected |")
        print(f"|:-----|:---------:|:--------:|")
        print(f"| keff | {keff:9} | {expectedKeff:8} |\n")

        if (abs(expectedKeff-keff)/expectedKeff <= error):
            sys.exit(0)
        else:
            sys.exit(1)



#==============================================================================*
