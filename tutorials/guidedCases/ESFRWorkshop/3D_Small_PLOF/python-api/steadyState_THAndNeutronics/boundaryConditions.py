"""
boundaryConditions.py — Initial and boundary conditions (coupled TH + neutronics)
==================================================================================

This file defines fields for both physics domains:
  - neutroRegion : neutron flux (defaultFlux)
  - fluidRegion  : coolant temperature (T), velocity (U), pressure (p, p_rgh)

All boundary conditions are the same as in the stand-alone cases.
The novelty here is that the two solvers exchange information:
  neutronics → TH : power density (where the heat is generated)
  TH → neutronics : coolant and fuel temperatures, coolant density (feedback)

This feedback loop is what makes the calculation "coupled" and physically
self-consistent: the power distribution depends on temperatures, which depend
on the power distribution.
"""

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc

from mesh import nMesh, thMesh

# ---------------------------------------------------------------------------
INLET_TEMPERATURE  = 668
INLET_VELOCITY     = 4.559114016
REFERENCE_PRESSURE = 100_000
COOLANT_DENSITY    = 860   # Sodium density at inlet conditions [kg/m³]
K_INIT             = 0.005
EPS_INIT           = 0.005
# ---------------------------------------------------------------------------

timeFolder0 = ffn.TimeFolder(time=0)

# ---------------------------------------------------------------------------
# Neutron flux
# ---------------------------------------------------------------------------
defaultFlux = ffn.Field("defaultFlux", region=nMesh.region)
defaultFlux.dimensions    = ffn.Dimension(default='flux')
defaultFlux.internalField = 1   # initial guess — normalised automatically by the eigenvalue solver
defaultFlux.set_boundary_condition('wall',   bc.FixedValue(0))
defaultFlux.set_boundary_condition('top',    bc.FixedValue(0))
defaultFlux.set_boundary_condition('bottom', bc.FixedValue(0))

# ---------------------------------------------------------------------------
# Temperature and density fields for cross-section interpolation
# (updated each iteration by the TH solver through the coupling)
# ---------------------------------------------------------------------------
TCool = ffn.Field("TCool", region=nMesh.region)
TCool.dimensions   = ffn.Dimension(default='T')
TCool.internalField = INLET_TEMPERATURE
TCool.set_boundary_condition('wall',   bc.ZeroGradient())
TCool.set_boundary_condition('top',    bc.ZeroGradient())
TCool.set_boundary_condition('bottom', bc.ZeroGradient())

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

# ---------------------------------------------------------------------------
# Coolant temperature
# ---------------------------------------------------------------------------
T = ffn.Field("T", region=thMesh.region)
T.dimensions    = ffn.Dimension(default='T')
T.internalField = INLET_TEMPERATURE
T.set_boundary_condition('wall',    bc.ZeroGradient())
T.set_boundary_condition('baffle0', bc.ZeroGradient())
T.set_boundary_condition('baffle1', bc.ZeroGradient())
T.set_boundary_condition('top',     bc.ZeroGradient())
T.set_boundary_condition('bottom',  bc.FixedValue(INLET_TEMPERATURE))

# ---------------------------------------------------------------------------
# Coolant velocity
# ---------------------------------------------------------------------------
U = ffn.Field("U", region=thMesh.region)
U.dimensions    = ffn.Dimension(default='U')
U.internalField = ffn.Vector(0, 0, INLET_VELOCITY)
U.set_boundary_condition('wall',    bc.Slip())
U.set_boundary_condition('baffle0', bc.Slip())
U.set_boundary_condition('baffle1', bc.Slip())
U.set_boundary_condition('top',     bc.ZeroGradient())
U.set_boundary_condition('bottom',  bc.FixedValue(U.internalField))

# ---------------------------------------------------------------------------
# Pressure
# ---------------------------------------------------------------------------
p = ffn.Field("p", region=thMesh.region)
p.dimensions    = ffn.Dimension(default='p')
p.internalField = REFERENCE_PRESSURE
p.set_boundary_condition('wall',    bc.ZeroGradient())
p.set_boundary_condition('baffle0', bc.Calculated(0))
p.set_boundary_condition('baffle1', bc.Calculated(0))
p.set_boundary_condition('top',     bc.FixedValue(REFERENCE_PRESSURE))
p.set_boundary_condition('bottom',  bc.ZeroGradient())

p_rgh = ffn.Field("p_rgh", region=thMesh.region)
p_rgh.dimensions    = ffn.Dimension(default='p')
p_rgh.internalField = REFERENCE_PRESSURE
p_rgh.set_boundary_condition('wall',    bc.ZeroGradient())
p_rgh.set_boundary_condition('baffle0', bc.ZeroGradient())
p_rgh.set_boundary_condition('baffle1', bc.ZeroGradient())
p_rgh.set_boundary_condition('top',     bc.FixedValue(REFERENCE_PRESSURE))
p_rgh.set_boundary_condition(
    'bottom', bc.FixedFluxExtrapolatedPressure(value=REFERENCE_PRESSURE, gradient=0)
)

# ---------------------------------------------------------------------------
# Buoyancy density ratio and turbulence fields (fluid region)
# ---------------------------------------------------------------------------
rhok = ffn.Field("rhok", region=thMesh.region)
rhok.dimensions   = ffn.Dimension(0, 0, 0, 0, 0, 0, 0)
rhok.internalField = 1
rhok.set_boundary_condition('wall',    bc.ZeroGradient())
rhok.set_boundary_condition('baffle0', bc.Calculated(0))
rhok.set_boundary_condition('baffle1', bc.Calculated(0))
rhok.set_boundary_condition('top',     bc.ZeroGradient())
rhok.set_boundary_condition('bottom',  bc.ZeroGradient())

k = ffn.Field("k", region=thMesh.region)
k.dimensions   = ffn.Dimension(0, 2, -2, 0, 0, 0, 0)
k.internalField = K_INIT
k.set_boundary_condition('wall',    bc.ZeroGradient())
k.set_boundary_condition('baffle0', bc.ZeroGradient())
k.set_boundary_condition('baffle1', bc.ZeroGradient())
k.set_boundary_condition('top',     bc.InletOutlet(inletValue=K_INIT, value=K_INIT))
k.set_boundary_condition('bottom',  bc.CustomPatch(parameters={
    "type": "turbulentIntensityKineticEnergyInlet", "intensity": 0.05,
    "value": f"uniform {K_INIT}",
}))

epsilon = ffn.Field("epsilon", region=thMesh.region)
epsilon.dimensions   = ffn.Dimension(0, 2, -3, 0, 0, 0, 0)
epsilon.internalField = EPS_INIT
epsilon.set_boundary_condition('wall',    bc.ZeroGradient())
epsilon.set_boundary_condition('baffle0', bc.ZeroGradient())
epsilon.set_boundary_condition('baffle1', bc.ZeroGradient())
epsilon.set_boundary_condition('top',     bc.InletOutlet(inletValue=EPS_INIT, value=EPS_INIT))
epsilon.set_boundary_condition('bottom',  bc.CustomPatch(parameters={
    "type": "turbulentMixingLengthDissipationRateInlet", "mixingLength": 0.001,
    "value": f"uniform {EPS_INIT}",
}))

nut = ffn.Field("nut", region=thMesh.region)
nut.dimensions   = ffn.Dimension(0, 2, -1, 0, 0, 0, 0)
nut.internalField = 0
nut.set_boundary_condition('wall',    bc.ZeroGradient())
nut.set_boundary_condition('baffle0', bc.ZeroGradient())
nut.set_boundary_condition('baffle1', bc.ZeroGradient())
nut.set_boundary_condition('top',     bc.Calculated(0))
nut.set_boundary_condition('bottom',  bc.Calculated(0))

alphat = ffn.Field("alphat", region=thMesh.region)
alphat.dimensions   = ffn.Dimension(1, -1, -1, 0, 0, 0, 0)
alphat.internalField = 0
alphat.set_boundary_condition('wall',    bc.ZeroGradient())
alphat.set_boundary_condition('baffle0', bc.Calculated(0))
alphat.set_boundary_condition('baffle1', bc.Calculated(0))
alphat.set_boundary_condition('top',     bc.Calculated(0))
alphat.set_boundary_condition('bottom',  bc.Calculated(0))

# ---------------------------------------------------------------------------
timeFolder0.append(defaultFlux)
timeFolder0.append(TCool)
timeFolder0.append(TFuel)
timeFolder0.append(TClad)
timeFolder0.append(rhoCool)
timeFolder0.append(T)
timeFolder0.append(U)
timeFolder0.append(p)
timeFolder0.append(p_rgh)
timeFolder0.append(rhok)
timeFolder0.append(k)
timeFolder0.append(epsilon)
timeFolder0.append(nut)
timeFolder0.append(alphat)
