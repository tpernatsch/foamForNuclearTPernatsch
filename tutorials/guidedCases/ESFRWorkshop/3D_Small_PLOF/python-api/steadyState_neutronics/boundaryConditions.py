"""
boundaryConditions.py — Initial conditions and boundary conditions (neutronics)
===============================================================================

This file defines the neutron flux field and its boundary conditions, as well
as the temperature and density fields that the neutronics solver uses to look
up cross-sections from the pre-computed tables.

NEUTRON FLUX (defaultFlux)
--------------------------
The neutron flux represents the number of neutrons passing through a unit area
per unit time. Solving the neutron diffusion equation for flux gives us the
spatial power distribution inside the reactor.

The initial value (internalField) is just a starting guess — the eigenvalue
solver will iterate until it finds the true flux shape (the fundamental mode).

BOUNDARY CONDITIONS
-------------------
At every external surface of the reactor model we apply a "black" boundary
condition (zero flux). This means we assume all neutrons that reach the
boundary are absorbed and never return — a conservative approximation that
avoids having to model the full surrounding reflector in detail.

TEMPERATURE AND DENSITY FIELDS (TFuel, TClad, rhoCool)
-------------------------------------------------------
The neutron cross-sections depend on the local temperature and coolant density.
These fields must be initialised so the solver knows which cross-section table
entries to use at the start of the calculation.

In this stand-alone neutronics case the temperatures are frozen at hot-zero-
power values (inlet temperature throughout).  In the coupled TH+neutronics
case these fields are updated every iteration by the thermal-hydraulics solver.
"""

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc

from mesh import nMesh   # neutronics mesh object defined in mesh.py

# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------
INLET_TEMPERATURE   = 668   # Sodium inlet temperature [K]
COOLANT_DENSITY     = 860   # Sodium density at inlet conditions [kg/m³]

# ---------------------------------------------------------------------------
# Initial time folder  (t = 0, i.e. starting conditions)
# ---------------------------------------------------------------------------
timeFolder0 = ffn.TimeFolder(time=0)

# ---------------------------------------------------------------------------
# Neutron flux field
# ---------------------------------------------------------------------------
# GeN-Foam uses a multi-group diffusion approach.  "defaultFlux" is the
# template applied to every energy group (flux0, flux1, …).
defaultFlux = ffn.Field("defaultFlux", region=nMesh.region)
defaultFlux.dimensions   = ffn.Dimension(default='flux')   # [neutrons / m² / s]
defaultFlux.internalField = 1   # initial guess; normalised automatically during eigenvalue iteration

# Black boundary condition: zero neutron flux at every outer surface
# (any neutron reaching the boundary is considered lost to the system)
defaultFlux.set_boundary_condition('wall',   bc.FixedValue(0))
defaultFlux.set_boundary_condition('top',    bc.FixedValue(0))
defaultFlux.set_boundary_condition('bottom', bc.FixedValue(0))

# ---------------------------------------------------------------------------
# Fuel temperature  (TFuel)
# ---------------------------------------------------------------------------
# Volume-averaged fuel pellet temperature used for Doppler cross-section
# interpolation.  Held at inlet temperature in this stand-alone case.
TFuel = ffn.Field("TFuel", region=nMesh.region)
TFuel.dimensions   = ffn.Dimension(default='T')   # [K]
TFuel.internalField = INLET_TEMPERATURE
TFuel.set_boundary_condition('wall',   bc.ZeroGradient())
TFuel.set_boundary_condition('top',    bc.ZeroGradient())
TFuel.set_boundary_condition('bottom', bc.ZeroGradient())

# ---------------------------------------------------------------------------
# Cladding temperature  (TClad)
# ---------------------------------------------------------------------------
# Volume-averaged cladding temperature used for cross-section interpolation.
TClad = ffn.Field("TClad", region=nMesh.region)
TClad.dimensions   = ffn.Dimension(default='T')   # [K]
TClad.internalField = INLET_TEMPERATURE
TClad.set_boundary_condition('wall',   bc.ZeroGradient())
TClad.set_boundary_condition('top',    bc.ZeroGradient())
TClad.set_boundary_condition('bottom', bc.ZeroGradient())

# ---------------------------------------------------------------------------
# Coolant density  (rhoCool)
# ---------------------------------------------------------------------------
# Sodium density used for coolant-void cross-section interpolation.
rhoCool = ffn.Field("rhoCool", region=nMesh.region)
rhoCool.dimensions   = ffn.Dimension(default='rho')   # [kg/m³]
rhoCool.internalField = COOLANT_DENSITY
rhoCool.set_boundary_condition('wall',   bc.ZeroGradient())
rhoCool.set_boundary_condition('top',    bc.ZeroGradient())
rhoCool.set_boundary_condition('bottom', bc.ZeroGradient())

# Collect all fields into the initial time folder
timeFolder0.append(defaultFlux)
timeFolder0.append(TFuel)
timeFolder0.append(TClad)
timeFolder0.append(rhoCool)
