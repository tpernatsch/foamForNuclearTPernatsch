"""
boundaryConditions.py — Initial conditions and boundary conditions (neutronics)
===============================================================================

TASK 1 — Set the neutron flux boundary conditions
--------------------------------------------------
The neutron flux (defaultFlux) is the key field in the neutronics calculation.
It represents the spatial distribution of neutrons in the reactor.

You need to set boundary conditions on all outer surfaces of the mesh.

Physical reasoning:
  - At every external boundary of our model, neutrons that reach the edge
    are assumed to be lost (absorbed by the surrounding shielding).
  - This is called a "black" or "vacuum" boundary condition.
  - Mathematically: zero neutron flux at the boundary → bc.FixedValue(0)

Fill in the three TODO lines below.
"""

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc

from mesh import nMesh

# ---------------------------------------------------------------------------
INLET_TEMPERATURE = 668   # Sodium inlet temperature [K]
COOLANT_DENSITY   = 860   # Sodium density at inlet conditions [kg/m³]
# ---------------------------------------------------------------------------

timeFolder0 = ffn.TimeFolder(time=0)

# Neutron flux field
defaultFlux = ffn.Field("defaultFlux", region=nMesh.region)
defaultFlux.dimensions    = ffn.Dimension(default='flux')
defaultFlux.internalField = 1   # initial guess — normalised automatically by the eigenvalue solver

# TODO: Set the black boundary condition (zero flux) on the three outer surfaces.
#       Use bc.FixedValue(0) for all of them.
#       The boundary names are: 'wall', 'top', 'bottom'
defaultFlux.set_boundary_condition('wall',   ???)   # outer radial boundary
defaultFlux.set_boundary_condition('top',    ???)   # top axial boundary
defaultFlux.set_boundary_condition('bottom', ???)   # bottom axial boundary

# ---------------------------------------------------------------------------
# Temperature and density fields for cross-section interpolation
# (held fixed at inlet conditions in this stand-alone neutronics case)
# ---------------------------------------------------------------------------
TFuel = ffn.Field("TFuel", region=nMesh.region)
TFuel.dimensions   = ffn.Dimension(default='T')
TFuel.internalField = INLET_TEMPERATURE
TFuel.set_boundary_condition('wall',   bc.ZeroGradient())
TFuel.set_boundary_condition('top',    bc.ZeroGradient())
TFuel.set_boundary_condition('bottom', bc.ZeroGradient())

TClad = ffn.Field("TClad", region=nMesh.region)
TClad.dimensions   = ffn.Dimension(default='T')
TClad.internalField = INLET_TEMPERATURE
TClad.set_boundary_condition('wall',   bc.ZeroGradient())
TClad.set_boundary_condition('top',    bc.ZeroGradient())
TClad.set_boundary_condition('bottom', bc.ZeroGradient())

rhoCool = ffn.Field("rhoCool", region=nMesh.region)
rhoCool.dimensions   = ffn.Dimension(default='rho')
rhoCool.internalField = COOLANT_DENSITY
rhoCool.set_boundary_condition('wall',   bc.ZeroGradient())
rhoCool.set_boundary_condition('top',    bc.ZeroGradient())
rhoCool.set_boundary_condition('bottom', bc.ZeroGradient())

timeFolder0.append(defaultFlux)
timeFolder0.append(TFuel)
timeFolder0.append(TClad)
timeFolder0.append(rhoCool)
