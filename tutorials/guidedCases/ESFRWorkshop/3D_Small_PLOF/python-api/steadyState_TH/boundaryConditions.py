"""
boundaryConditions.py — Initial conditions and boundary conditions (thermal-hydraulics)
========================================================================================

This file defines the initial state and boundary conditions for all fluid fields.

FIELDS
------
T                    — coolant temperature [K]
U                    — coolant velocity vector [m/s]
p                    — thermodynamic pressure [Pa]
p_rgh                — dynamic pressure (p minus the hydrostatic head ρgh) [Pa]
                       GeN-Foam uses p_rgh internally for numerical stability in
                       buoyancy-driven flows; p is kept for post-processing.
powerDensityStructure — nuclear power deposited in the solid structures [W/m³]
                       (zero in this stand-alone TH case; updated by neutronics
                       in the coupled case)
rhok                 — normalised buoyancy density field (dimensionless) used
                       internally by the solver to handle buoyancy
k                    — turbulent kinetic energy [m²/s²]   ─┐ only active if
epsilon              — turbulent dissipation rate [m²/s³]  │ turbulenceModel
nut                  — turbulent viscosity [m²/s]          │ is switched from
alphat               — turbulent thermal diffusivity [kg/m/s] ─┘ laminar to k-ε

BOUNDARY CONDITIONS OVERVIEW
-----------------------------
  'bottom'  — sodium inlet (bottom of the diagrid plenum)
  'top'     — sodium outlet (top of upper reflector / gas plenums)
  'wall'    — outer radial boundary of the core (no-flow)
  'baffle0', 'baffle1' — internal walls between adjacent assemblies

Temperature (T):
  - Inlet (bottom): fixed at the inlet temperature (Dirichlet)
  - All other boundaries: zero-gradient (Neumann) — sodium temperature
    adjusts freely, no constraint imposed at walls or outlet

Velocity (U):
  - Inlet (bottom): fixed upward velocity (Dirichlet)
  - Outlet (top): zero-gradient — velocity adjusts at the outlet
  - Walls and baffles: slip condition — tangential velocity allowed,
    no normal flow through the wall (appropriate for a porous-medium model)

Pressure (p, p_rgh):
  - Outlet (top): fixed reference pressure (Dirichlet)
  - Inlet (bottom): zero-gradient or extrapolated flux (Neumann)
  - Walls and baffles: zero-gradient
"""

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc

from mesh import thMesh

# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------
INLET_TEMPERATURE  = 668           # sodium inlet temperature [K]
INLET_VELOCITY     = 4.559114016   # mean axial inlet velocity [m/s]
                                   # (corresponds to nominal mass flow rate)
REFERENCE_PRESSURE = 100_000       # reference pressure at the outlet [Pa]

# Turbulence initial values (only used if switching to k-ε model)
K_INIT       = 0.005   # turbulent kinetic energy initial guess [m²/s²]
EPS_INIT     = 0.005   # dissipation rate initial guess [m²/s³]

# ---------------------------------------------------------------------------
timeFolder0 = ffn.TimeFolder(time=0)

# ---------------------------------------------------------------------------
# Temperature field
# ---------------------------------------------------------------------------
T = ffn.Field("T", region=thMesh.region)
T.dimensions   = ffn.Dimension(default='T')
T.internalField = INLET_TEMPERATURE   # core starts cold (hot-zero-power)

T.set_boundary_condition('wall',    bc.ZeroGradient())
T.set_boundary_condition('baffle0', bc.ZeroGradient())
T.set_boundary_condition('baffle1', bc.ZeroGradient())
T.set_boundary_condition('top',     bc.ZeroGradient())
T.set_boundary_condition('bottom',  bc.FixedValue(INLET_TEMPERATURE))

# ---------------------------------------------------------------------------
# Velocity field
# ---------------------------------------------------------------------------
U = ffn.Field("U", region=thMesh.region)
U.dimensions   = ffn.Dimension(default='U')
U.internalField = ffn.Vector(0, 0, INLET_VELOCITY)   # upward axial flow

U.set_boundary_condition('wall',    bc.Slip())
U.set_boundary_condition('baffle0', bc.Slip())
U.set_boundary_condition('baffle1', bc.Slip())
U.set_boundary_condition('top',     bc.ZeroGradient())
U.set_boundary_condition('bottom',  bc.FixedValue(U.internalField))

# ---------------------------------------------------------------------------
# Thermodynamic pressure
# ---------------------------------------------------------------------------
p = ffn.Field("p", region=thMesh.region)
p.dimensions   = ffn.Dimension(default='p')
p.internalField = REFERENCE_PRESSURE

p.set_boundary_condition('wall',    bc.ZeroGradient())
p.set_boundary_condition('baffle0', bc.Calculated(0))
p.set_boundary_condition('baffle1', bc.Calculated(0))
p.set_boundary_condition('top',     bc.FixedValue(REFERENCE_PRESSURE))
p.set_boundary_condition('bottom',  bc.ZeroGradient())

# ---------------------------------------------------------------------------
# Modified pressure (p_rgh = p − ρgh)
# ---------------------------------------------------------------------------
p_rgh = ffn.Field("p_rgh", region=thMesh.region)
p_rgh.dimensions   = ffn.Dimension(default='p')
p_rgh.internalField = REFERENCE_PRESSURE

p_rgh.set_boundary_condition('wall',    bc.ZeroGradient())
p_rgh.set_boundary_condition('baffle0', bc.ZeroGradient())
p_rgh.set_boundary_condition('baffle1', bc.ZeroGradient())
p_rgh.set_boundary_condition('top',     bc.FixedValue(REFERENCE_PRESSURE))
p_rgh.set_boundary_condition(
    'bottom', bc.FixedFluxExtrapolatedPressure(value=REFERENCE_PRESSURE, gradient=0)
)

# ---------------------------------------------------------------------------
# Buoyancy density ratio (rhok)
# ---------------------------------------------------------------------------
# Internal field for GeN-Foam's buoyancy algorithm; initialised to 1 (neutral).
rhok = ffn.Field("rhok", region=thMesh.region)
rhok.dimensions   = ffn.Dimension(0, 0, 0, 0, 0, 0, 0)  # dimensionless
rhok.internalField = 1
rhok.set_boundary_condition('wall',    bc.ZeroGradient())
rhok.set_boundary_condition('baffle0', bc.Calculated(0))
rhok.set_boundary_condition('baffle1', bc.Calculated(0))
rhok.set_boundary_condition('top',     bc.ZeroGradient())
rhok.set_boundary_condition('bottom',  bc.ZeroGradient())

# ---------------------------------------------------------------------------
# Turbulence fields (k-ε model)
# ---------------------------------------------------------------------------
# The default turbulence model is LAMINAR — these fields are not solved in that
# mode.  They are included here so the case runs without errors if students
# switch turbulenceProperties.simulationType to 'RAS' (k-ε or similar).

# Turbulent kinetic energy k [m²/s²]
k = ffn.Field("k", region=thMesh.region)
k.dimensions   = ffn.Dimension(0, 2, -2, 0, 0, 0, 0)
k.internalField = K_INIT
k.set_boundary_condition('wall',    bc.ZeroGradient())
k.set_boundary_condition('baffle0', bc.ZeroGradient())
k.set_boundary_condition('baffle1', bc.ZeroGradient())
k.set_boundary_condition('top',     bc.InletOutlet(inletValue=K_INIT, value=K_INIT))
k.set_boundary_condition('bottom',  bc.CustomPatch(parameters={
    "type":      "turbulentIntensityKineticEnergyInlet",
    "intensity": 0.05,   # 5% turbulence intensity at inlet
    "value":     f"uniform {K_INIT}",
}))

# Turbulent dissipation rate ε [m²/s³]
epsilon = ffn.Field("epsilon", region=thMesh.region)
epsilon.dimensions   = ffn.Dimension(0, 2, -3, 0, 0, 0, 0)
epsilon.internalField = EPS_INIT
epsilon.set_boundary_condition('wall',    bc.ZeroGradient())
epsilon.set_boundary_condition('baffle0', bc.ZeroGradient())
epsilon.set_boundary_condition('baffle1', bc.ZeroGradient())
epsilon.set_boundary_condition('top',     bc.InletOutlet(inletValue=EPS_INIT, value=EPS_INIT))
epsilon.set_boundary_condition('bottom',  bc.CustomPatch(parameters={
    "type":         "turbulentMixingLengthDissipationRateInlet",
    "mixingLength": 0.001,   # characteristic mixing length [m]
    "value":        f"uniform {EPS_INIT}",
}))

# Turbulent kinematic viscosity nut [m²/s]  — computed by the turbulence model
nut = ffn.Field("nut", region=thMesh.region)
nut.dimensions   = ffn.Dimension(0, 2, -1, 0, 0, 0, 0)
nut.internalField = 0
nut.set_boundary_condition('wall',    bc.ZeroGradient())
nut.set_boundary_condition('baffle0', bc.ZeroGradient())
nut.set_boundary_condition('baffle1', bc.ZeroGradient())
nut.set_boundary_condition('top',     bc.Calculated(0))
nut.set_boundary_condition('bottom',  bc.Calculated(0))

# Turbulent thermal diffusivity alphat [kg/m/s]  — computed by the turbulence model
alphat = ffn.Field("alphat", region=thMesh.region)
alphat.dimensions   = ffn.Dimension(1, -1, -1, 0, 0, 0, 0)
alphat.internalField = 0
alphat.set_boundary_condition('wall',    bc.ZeroGradient())
alphat.set_boundary_condition('baffle0', bc.Calculated(0))
alphat.set_boundary_condition('baffle1', bc.Calculated(0))
alphat.set_boundary_condition('top',     bc.Calculated(0))
alphat.set_boundary_condition('bottom',  bc.Calculated(0))

# ---------------------------------------------------------------------------
timeFolder0.append(T)
timeFolder0.append(U)
timeFolder0.append(p)
timeFolder0.append(p_rgh)
timeFolder0.append(rhok)
timeFolder0.append(k)
timeFolder0.append(epsilon)
timeFolder0.append(nut)
timeFolder0.append(alphat)
