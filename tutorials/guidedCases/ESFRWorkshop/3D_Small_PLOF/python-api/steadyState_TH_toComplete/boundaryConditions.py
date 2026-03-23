"""
boundaryConditions.py — Initial conditions and boundary conditions  [TO COMPLETE]
==================================================================================

TASK 1 — Set the temperature boundary conditions
-------------------------------------------------
The temperature field T represents the sodium coolant temperature.

Physical reasoning:
  - At the inlet (bottom), sodium enters at a fixed temperature → FixedValue
  - At all other boundaries, temperature is not imposed — it results from
    the flow and heat transfer → ZeroGradient (no constraint)

Boundaries: 'wall', 'baffle0', 'baffle1', 'top', 'bottom'
Inlet temperature: INLET_TEMPERATURE (defined below)

Fill in the five bc.??? calls for T.

TASK 2 — Set the velocity boundary conditions
----------------------------------------------
The velocity field U is a 3D vector representing the coolant flow direction
and speed.  Sodium flows upward (positive z-direction).

Physical reasoning:
  - At the inlet (bottom): fixed upward velocity → bc.FixedValue(U.internalField)
  - At the outlet (top): the velocity adjusts freely → bc.ZeroGradient()
  - At walls and baffles: the coolant cannot flow through the wall,
    but can slide along it (porous-medium assumption) → bc.Slip()

Fill in the five bc.??? calls for U.
"""

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc

from mesh import thMesh

# ---------------------------------------------------------------------------
INLET_TEMPERATURE  = 668
INLET_VELOCITY     = 4.559114016   # nominal upward velocity [m/s]
REFERENCE_PRESSURE = 100_000       # [Pa]
K_INIT             = 0.005         # turbulent kinetic energy initial value [m²/s²]
EPS_INIT           = 0.005         # turbulent dissipation rate initial value [m²/s³]
# ---------------------------------------------------------------------------

timeFolder0 = ffn.TimeFolder(time=0)

# ---------------------------------------------------------------------------
# Temperature
# ---------------------------------------------------------------------------
T = ffn.Field("T", region=thMesh.region)
T.dimensions    = ffn.Dimension(default='T')
T.internalField = INLET_TEMPERATURE

# TODO: set boundary conditions for T
#   Hint: use bc.FixedValue(INLET_TEMPERATURE) at the inlet,
#         and bc.ZeroGradient() everywhere else
T.set_boundary_condition('wall',    ???)
T.set_boundary_condition('baffle0', ???)
T.set_boundary_condition('baffle1', ???)
T.set_boundary_condition('top',     ???)
T.set_boundary_condition('bottom',  ???)   # inlet: fixed temperature

# ---------------------------------------------------------------------------
# Velocity
# ---------------------------------------------------------------------------
U = ffn.Field("U", region=thMesh.region)
U.dimensions    = ffn.Dimension(default='U')
U.internalField = ffn.Vector(0, 0, INLET_VELOCITY)

# TODO: set boundary conditions for U
#   Hint: use bc.Slip() at walls/baffles, bc.ZeroGradient() at outlet (top),
#         bc.FixedValue(U.internalField) at the inlet (bottom)
U.set_boundary_condition('wall',    ???)
U.set_boundary_condition('baffle0', ???)
U.set_boundary_condition('baffle1', ???)
U.set_boundary_condition('top',     ???)   # outlet: free
U.set_boundary_condition('bottom',  ???)   # inlet: fixed velocity

# ---------------------------------------------------------------------------
# Pressure fields  (provided — no task here)
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
# Additional fields  (provided — no task here)
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

timeFolder0.append(T)
timeFolder0.append(U)
timeFolder0.append(p)
timeFolder0.append(p_rgh)
timeFolder0.append(rhok)
timeFolder0.append(k)
timeFolder0.append(epsilon)
timeFolder0.append(nut)
timeFolder0.append(alphat)
