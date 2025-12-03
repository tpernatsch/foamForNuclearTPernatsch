
#==============================================================================*
# Imports

import openmc
from foamForNuclear.openmcTools import MultiGroupXS, MultiGroupXSManager


#==============================================================================*
# Model

model = openmc.Model()


# Define materials
#-----------------

fuel = openmc.Material(name='UO2 (2.4%)')
fuel.temperature = 574
fuel.set_density('g/cm3', 10.29769)
fuel.add_nuclide('U234', 4.4843e-6)
fuel.add_nuclide('U235', 5.5815e-4)
fuel.add_nuclide('U238', 2.2408e-2)
fuel.add_nuclide('O16', 4.5829e-2)

clad = openmc.Material(name='Zircaloy')
clad.temperature = 574
clad.set_density('g/cm3', 6.55)
clad.add_nuclide('Zr90', 2.1827e-2)
clad.add_nuclide('Zr91', 4.7600e-3)
clad.add_nuclide('Zr92', 7.2758e-3)
clad.add_nuclide('Zr94', 7.3734e-3)
clad.add_nuclide('Zr96', 1.1879e-3)

water = openmc.Material(name='Borated water')
water.temperature = 574
water.set_density('g/cm3', 0.740582)
water.add_nuclide('H1', 4.9457e-2)
water.add_nuclide('O16', 2.4672e-2)
water.add_nuclide('B10', 8.0042e-6)
water.add_nuclide('B11', 3.2218e-5)
water.add_s_alpha_beta('c_H_in_H2O')

# Define the materials file
model.materials = (fuel, clad, water)


# Define the geometry
#--------------------

# Instantiate ZCylinder surfaces
pitch = 2 * 0.63
fuel_or = openmc.ZCylinder(x0=0, y0=0, r=0.39218, name='Fuel OR')
clad_or = openmc.ZCylinder(x0=0, y0=0, r=0.45720, name='Clad OR')
left = openmc.XPlane(x0=-pitch/2, name='left', boundary_type='reflective')
right = openmc.XPlane(x0=pitch/2, name='right', boundary_type='reflective')
front = openmc.YPlane(y0=-pitch/2, name='front', boundary_type='reflective')
back = openmc.YPlane(y0=pitch/2, name='back', boundary_type='reflective')

# Instantiate Cells
fuel_pin = openmc.Cell(name='Fuel', fill=fuel)
cladding = openmc.Cell(name='Cladding', fill=clad)
water = openmc.Cell(name='Water', fill=water)

# Use surface half-spaces to define regions
fuel_pin.region = -fuel_or
cladding.region = +fuel_or & -clad_or
water.region = +clad_or & +left & -right & +front & -back

# Create root universe
model.geometry.root_universe = openmc.Universe(0, name='root universe')
model.geometry.root_universe.add_cells([fuel_pin, cladding, water])


# Settings
#---------

model.settings.batches = 1000
model.settings.inactive = 100
model.settings.particles = 10000
model.settings.source = openmc.IndependentSource(
    space=openmc.stats.Box([-pitch/2, -pitch/2, -1],
                            [pitch/2, pitch/2, 1]),
    constraints={'fissionable': True}
)
model.settings.temperature = {
    'method': 'interpolation',
    'multipole': True
}

plot = openmc.Plot.from_geometry(model.geometry)
plot.pixels = (300, 300)
plot.color_by = 'material'
model.plots.append(plot)


# Tallies
#--------

# Instantiate an empty Tallies object
tallies = openmc.Tallies()

# Instantiate a 2-group EnergyGroups object
energyGroups = openmc.mgxs.EnergyGroups([0., 0.625, 20.0e6]) # eV

# Instantiate a 7-group of delayed neutrons
delayedGroups = list(range(1,7))

# Create the MultiGroupXSManager on the master universe of the geometry (uReactor)
mgxsManager = MultiGroupXSManager(
    model.geometry.root_universe, energyGroups, delayedGroups
)

# Add a new MGXS specific
mgxsManager.AddDomain(
    MultiGroupXS(
        fuel_pin,              # Domain: openmc.Cell, openmc.Universe, openmc.Material
        energyGroups,       # Energy group (openmc.mgxs.EnergyGroups)
        delayedGroups,      # Delayed neutrons
        fuelFraction=1,     # Fuel fraction inside the domain
        genfoamName='fuel'  # Name of the cellZone in GeN-Foam
    )
)
mgxsManager.AddDomain(
    MultiGroupXS(
        cladding,              # Domain: openmc.Cell, openmc.Universe, openmc.Material
        energyGroups,       # Energy group (openmc.mgxs.EnergyGroups)
        delayedGroups,      # Delayed neutrons
        fuelFraction=0,     # Fuel fraction inside the domain
        genfoamName='cladding'  # Name of the cellZone in GeN-Foam
    )
)
mgxsManager.AddDomain(
    MultiGroupXS(
        water,              # Domain: openmc.Cell, openmc.Universe, openmc.Material
        energyGroups,       # Energy group (openmc.mgxs.EnergyGroups)
        delayedGroups,      # Delayed neutrons
        fuelFraction=0,     # Fuel fraction inside the domain
        genfoamName='water'  # Name of the cellZone in GeN-Foam
    )
)

# Add all the MGXS tallies to the main tally object of OpenMC
mgxsManager.AddToTallies(tallies)


model.tallies = tallies

# Export model
model.export_to_xml()


#==============================================================================*
# Run

openmc.run(threads=10)


#==============================================================================*