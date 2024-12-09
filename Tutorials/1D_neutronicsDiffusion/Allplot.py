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


def plotFlux(
        timeStep: str,
        fieldName: str,
        isAddExtrapolationDistance: bool=False
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

    fuelLength = 1.5

    # Normalize the flux
    normFlux = [flux0_ / max(flux0) for flux0_ in flux0]

    # Analytic solution
    fluxTh = lambda z_: np.sin(z_ * np.pi/fuelLength)

    # Compute the relative error between GeN-Foam and the analytic solution
    relError = [(normFlux_ / fluxTh(z_) - 1) * 100 for z_, normFlux_ in zip(z, normFlux)]

    # Plot
    fig, axes = plt.subplots(2, 1, figsize=(6, 6))
    axFlux, axError = axes.flatten()

    ymax = max(normFlux) * 1.05

    axFlux.plot(z, normFlux, label="GeN-Foam")
    axFlux.plot(z, fluxTh(z), label="Analytic", ls="--")
    axFlux.set_ylabel("Neutron flux [a.u]")
    axFlux.set_ylim((0, ymax))
    axFlux.legend()

    axError.plot(z, relError)
    axError.set_ylabel("Relative error [%]")

    for ax in axes:
        ax.set_xlabel("Axial position [m]")
        ax.set_xlim(min(z), max(z))


    if (isAddExtrapolationDistance):
        # Estimated extrapolated limit
        d = 0.1275*3*0.7104

        print(f"\nEstimated extrapolation length = {d} m\n")

        # Right extrapolation
        slope = (normFlux[-1]-normFlux[-2])/(z[-1]-z[-2])
        zExtra = np.linspace(z[-1], z[-1]+d, 100)
        extFlux = lambda z_: slope*(z_ - z[-1]) + normFlux[-1]
        axFlux.plot(
            zExtra, extFlux(zExtra),
            ls=':', color='tab:grey',
            label="Extrapolated limit"
        )

        # Left extrapolation
        slope = (normFlux[0]-normFlux[1])/(z[0]-z[1])
        zExtra = np.linspace(z[0]-d, z[0], 100)
        extFlux = lambda z_: slope*(z_ - z[0]) + normFlux[0]
        axFlux.plot(zExtra, extFlux(zExtra), ls=':', color='tab:grey')

        # Core boundaries
        axFlux.vlines(
            x=[0, fuelLength], ymin=0, ymax=ymax,
            color="tab:grey", alpha=0.5
        )

        # Reset axis limits
        axError.hlines(xmin=min(z)-d, xmax=max(z)+d, y=[0], alpha=0.5, color="tab:grey")
        for ax in axes:
            ax.set_xlim(min(z) - d, max(z) + d)

        axFlux.legend()

    else:
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
    isTest = bool(sys.argv[3])

    # Plot and extract keff
    if (not isTest):
        plotFlux(timeStep, fieldName, isAddExtrapolationDistance=False)

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
