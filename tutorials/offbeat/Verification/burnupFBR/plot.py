import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import numpy as np
import pandas as pd

# import scienceplots 
# plt.style.use("ieee")


pathBurnup = "postProcessing/burnup/0/volFieldValue.dat"
pathNuclides = "postProcessing/nuclides/0/volFieldValue.dat"
# serpentPath = "../_serpent/nuclides.csv"
serpentPath = "./resultsSERPENT.csv"

serpentData = pd.read_csv(serpentPath)
dataBurnup = pd.read_csv(pathBurnup, skiprows=3,  sep="\t")
dataNuclides = pd.read_csv(pathNuclides, skiprows=3,  sep="\t")

# Get list of nuclides and burnup OFFBEAT
nuclides = [w.replace("volAverage(","").replace(")","").replace("N_","") for w in dataNuclides.columns[1:]]
burnup = dataBurnup["volAverage(Bu_kgU)"].to_numpy()/1000

# Get list of elements 
elements = np.unique([''.join([char for char in item if not char.isdigit()]) for item in nuclides])
# Remove metastable americium
elements = [elem for elem in elements if elem != 'Amm']

colors = ['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e']

totPuSERPENT = 0
totPuOFFBEAT = 0

# --- Plot
for element in elements:

    # Mask nuclides list for this element
    theseNuclides = [n for n in nuclides if element in n]

    minY = 1e300
    maxY = 0

    cm = 1/2.54
    fig, ax = plt.subplots(1, figsize=(8*cm,8*cm))


    for i, nuclide in enumerate(theseNuclides):

        color = colors[i]

        # ---  Serpent
        MM = "".join([c for c in nuclide if c.isdigit()])
        ncl = nuclide.replace(MM,"")
        if (ncl=="Amm"): 
            ncl="Am"
            MM="242,m"
        string = fr"${{^{{{MM}}}}}{ncl}$"
        ax.plot(serpentData["# Bu_MWdkg"], serpentData[f"{nuclide}_dens_1/cm3"], lw=1, alpha=0.4, 
                ls='', marker='o', markersize=5, label=string, color = color)

        minY = min(minY, min(serpentData[f"{nuclide}_dens_1/cm3"][1:]))
        maxY = max(maxY, max(serpentData[f"{nuclide}_dens_1/cm3"]))

        # ---  OFFBEAT
        concentration = dataNuclides[f"volAverage(N_{nuclide})"]

        ax.plot(burnup, concentration, ls="-", color = color)

        # --- compute tot Pu
        if ("Pu" in element):
            totPuOFFBEAT += concentration.iloc[-1]
            totPuSERPENT += serpentData[f"{nuclide}_dens_1/cm3"].iloc[-1]

    ax.set_xlabel("Burnup (MWd/kgHM)")
    ax.set_ylabel(r"Nuclide Density (cm$^{-3}$)")
    # ax.legend(frameon=False, fontsize = 8)
    # ax.legend(frameon=False, fontsize = 8, ncol=len(theseNuclides), bbox_to_anchor=(0.5, 1), loc="center")
    # ax.legend(frameon=False, fontsize = 8, ncol=1, bbox_to_anchor=(1, 0.5), loc="center right")
    legend = ax.legend(frameon=False, ncol=2, fontsize = 9, loc="lower right")
    ax.grid(ls=":")

    ax.set_yscale("log")
    # ax.set_xlim(-10, 300)
    if element=="Pu":
        ax.set_ylim(minY/150, 10 * maxY)
    else:
        ax.set_ylim(minY/10, 10 * maxY)
            

    # Set inward ticks
    ax.tick_params(direction='in', which='both')

    fig.tight_layout()
    fig.savefig(f"plt_ch4_SERPENTverif_{element}.png")

# Print in file diff in Pu between SERPENT and OFFBEAT
with open("results","w") as f:
    error = (totPuOFFBEAT-totPuSERPENT)/totPuSERPENT
    f.write(f"Error Tot. [Pu] OFFBEAT/SERPENT: {error:.2e} %\n")
    if error*100 < 0.5:
        f.write(f"Test passed")
    else:
        f.write(f"Test failed")

