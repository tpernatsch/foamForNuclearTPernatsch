# MultiGroupXS Python Package

Author: Thomas Guilbaud, 2022/07/17

This package aims to generate multi-group cross-sections (MGXS) for the deterministic multi-physics code [GeN-Foam](https://gitlab.com/foam-for-nuclear/GeN-Foam) using the Monte Carlo code [OpenMC](https://docs.openmc.org/en/stable/).

Package updated to OpenMC version 0.13.4.


## How to install

Requires the installation of OpenMC with the `pip3 install .` command, and the numpy package.

```bash
python3 setup.py install --user
# or
pip3 install .
```

Execute the following command to check the installation.

```bash
pip show MultiGroupXS
```


## How to use

In any Python script, import the MultiGroupXS package as below (see [KiwiB4E/core/OpenMC/simple/main.py](KiwiB4E/core/OpenMC/simple/main.py)):
```Python
from MultiGroupXS import *
```

It is now possible to use `MultiGroupXS` and `MultiGroupXSManager`. Here is an example of use:
```Python
# Basic import
import openmc
from MultiGroupXS import *

# Instantiate an empty Tallies object
tallies = openmc.Tallies()

# Instantiate a 2-group EnergyGroups object
energyGroups = openmc.mgxs.EnergyGroups([0., 0.625, 20.0e6]) # eV

# Instantiate a 7-group of delayed neutrons
delayedGroups = list(range(1,7))

# Create the MultiGroupXSManager on the master universe of the geometry (uReactor)
mgxsManager = MultiGroupXSManager(
    uReactor, energyGroups, delayedGroups
)

# Add a new MGXS specific to cCore
mgxsManager.AddDomain(
    MultiGroupXS(
        cCore,              # Domain: openmc.Cell, openmc.Universe, openmc.Material
        energyGroups,       # Energy group (openmc.mgxs.EnergyGroups)
        delayedGroups,      # Delayed neutrons
        fuelFraction=1,     # Fuel fraction inside the domain
        genfoamName='core'  # Name of the cellZone in GeN-Foam
    )
)

# Add a new MGXS specific to cCore
mgxsManager.AddDomain(
    MultiGroupXS(
        mReflector,
        energyGroups,
        delayedGroups,
        fuelFraction=0,
        genfoamName='reflector'
    )
)

# Add all the MGXS tallies to the main tally object of OpenMC
mgxsManager.AddToTallies(tallies)

# Export tallies in XML for simulation
tallies.export_to_xml()


# Run the simulation ...


# Load the StatePoint result file
batches = nInactiveCycles + nActiveCycles
sp = openmc.StatePoint(f"statepoint.{batches}.h5")

# Read and generate the MGXS for GeN-Foam
mgxsManager.LoadFromStatepoint(sp)

# Generate the nuclearData file for GeN-Foam
mgxsManager.PrintForGeNFoam('nuclearData')
```


## How to uninstall

Execute the following command to uninstall the MultiGroupXS package.

```bash
pip3 uninstall MultiGroupXS
```
