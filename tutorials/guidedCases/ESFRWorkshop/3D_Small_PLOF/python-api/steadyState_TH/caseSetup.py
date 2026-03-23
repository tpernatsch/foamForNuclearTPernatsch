"""
caseSetup.py — Thermal-hydraulics solver configuration
======================================================

This file configures the thermal-hydraulics solver for the ESFR coolant.

SODIUM COOLANT PROPERTIES
--------------------------
Liquid sodium is used as coolant in the ESFR.  Its thermophysical properties
(density, heat capacity, viscosity, thermal conductivity) depend on temperature
and are given as polynomial fits valid in the reactor operating range.

POROUS MEDIUM APPROACH
-----------------------
The reactor core is not meshed pin-by-pin.  Instead, each hexagonal assembly
is treated as a porous medium: the solid structures (fuel pins, cladding,
wrapper tube) are "homogenised" into effective friction and heat-transfer
coefficients.  This makes the simulation tractable on a workstation.

Key parameters per assembly:
  volumeFraction  — fraction of the assembly cross-section occupied by coolant
  Dh              — hydraulic diameter of the coolant channel [m]

FUEL PIN MODEL
--------------
For fuel assemblies (inner/outer core), a 1D radial pin model is solved inside
each assembly to compute the fuel and cladding temperatures.  These are then
passed to the neutronics solver as feedback parameters.
"""

import foamForNuclear as ffn

from mesh               import thMesh
from boundaryConditions import timeFolder0, INLET_TEMPERATURE

# ---------------------------------------------------------------------------
# Thermal-hydraulics solver  (single-phase sodium)
# ---------------------------------------------------------------------------
thSolver = ffn.ThermalHydraulicsSolver(
    region            = "fluidRegion",
    solver            = "onePhase",
    removeBaffles     = True,   # merge internal assembly-wall baffles after decomposition
    mesh              = thMesh,
    isSetFvSolutionToDefault = False
)

# ---------------------------------------------------------------------------
# Sodium thermophysical properties  (temperature-dependent polynomials)
# ---------------------------------------------------------------------------
# Sodium properties evaluated at the inlet temperature as a starting point;
# the solver updates them each iteration based on the local temperature field.
thSolver.thermophysicalProperties = ffn.thermophysicalProperty.SodiumPolynomial()
thSolver.thermophysicalProperties.Cp    = ffn.Polynome(1646.97, -0.831587,  4.31182e-04)
thSolver.thermophysicalProperties.rho   = ffn.Polynome(860)        # density [kg/m³]
thSolver.thermophysicalProperties.mu    = ffn.Polynome(2e-4)       # dynamic viscosity [Pa·s]
thSolver.thermophysicalProperties.kappa = ffn.Polynome(60)         # thermal conductivity [W/m/K]

# Turbulence model: laminar for this simplified workshop case
thSolver.turbulenceProperties.simulationType = 'laminar'

# ---------------------------------------------------------------------------
# Structure zones  (all non-fuel assemblies — treated as passive heat sinks)
# ---------------------------------------------------------------------------
# These assemblies conduct heat into the coolant but have no internal power source.
STRUCTURE_ZONES = [
    "diagrid", "upperGasPlenum", "upperReflector",
    "lowerGasPlenum", "lowerReflector",
    "radialReflector", "follower", "controlRod"
]

passiveStructures = ffn.StructureProperty(
    STRUCTURE_ZONES,
    volumeFraction = 0.718520968,   # coolant volume fraction [-]
    Dh             = 0.00365        # hydraulic diameter [m]
)
passiveStructures.add_passive_structure(
    volumetricArea = 5,         # heat-transfer surface area per unit volume [m²/m³]
    rhoCp          = 4.8e6,     # volumetric heat capacity of the structure [J/m³/K]
    T              = INLET_TEMPERATURE
)
thSolver.add_structure_property(passiveStructures)

# ---------------------------------------------------------------------------
# Fuel assembly zones  (inner and outer core — nuclear heat source)
# ---------------------------------------------------------------------------
# The pin model resolves the temperature profile across the fuel pellet,
# gap, and cladding for each assembly.
CORE_ZONES = ["innerCore", "outerCore"]

core = ffn.StructureProperty(
    zones          = CORE_ZONES,
    volumeFraction = 0.718520968,
    Dh             = 0.00365
)
core.add_power_model(
    ffn.NuclearFuelPin(
        fuelInnerRadius = 0.0012,    # inner fuel pellet radius [m] (central void)
        fuelOuterRadius = 0.004715,  # outer fuel pellet radius [m]
        cladInnerRadius = 0.004865,  # inner cladding radius [m]
        cladOuterRadius = 0.005365,  # outer cladding radius [m]
        fuelMeshSize    = 30,        # radial mesh points in fuel
        cladMeshSize    = 5,         # radial mesh points in cladding
        fuelRho  = 10480, fuelCp  = 250, fuelK  = 3,     # fuel thermo-physical properties
        cladRho  = 7500,  cladCp  = 500, cladK  = 20,    # cladding thermo-physical properties
        gapH     = 3000,             # fuel-clad gap conductance [W/m²/K]
        fuelT    = INLET_TEMPERATURE,
        cladT    = INLET_TEMPERATURE,
        powerDensity = 410031689.2970038 
    )
)
core.add_passive_structure(volumetricArea=5, rhoCp=4.8e6, T=INLET_TEMPERATURE)
thSolver.add_structure_property(core)

# ---------------------------------------------------------------------------
# Friction (drag) model — Reynolds-power law correlation
# ---------------------------------------------------------------------------
# Models the pressure drop of sodium flowing through the assembly bundle.
# dP/dz ∝ coeff × Re^exp × (ρ U²/2Dh)
ALL_ZONES = CORE_ZONES + STRUCTURE_ZONES
thSolver.add_drag_model(
    ffn.ReynoldsPower(coeff=0.687, exp=-0.25, zones=ALL_ZONES)
)

# ---------------------------------------------------------------------------
# Heat transfer model — by flow regime (laminar/turbulent)
# ---------------------------------------------------------------------------
# The Nusselt number (dimensionless heat transfer coefficient) depends on
# whether the flow is laminar or turbulent (characterised by the Reynolds number).
heatTransfer = ffn.HeatTransferByRegime("lamTurb", zones=ALL_ZONES)
heatTransfer.add_regime(
    ffn.NusseltReynoldsPrandtlPower(   # laminar: Nu = 4 (constant for slug flow in Na)
        const=4, coeff=0, expRe=0, expPr=0, zones=['laminar']
    )
)
heatTransfer.add_regime(
    ffn.NusseltReynoldsPrandtlPower(   # turbulent: Kazimi-Carelli correlation
        const=4.82, coeff=0.0185, expRe=0.827, expPr=0.827, zones=['turbulent']
    )
)
thSolver.add_heat_transfer_model(heatTransfer)

# Regime map: transition from laminar to turbulent flow
lamTurb = ffn.RegimeMapOneParameter(name="lamTurb", parameter="Re")
lamTurb.add_regime("laminar",   0,    1000)
lamTurb.add_regime("turbulent", 2300, 2301)   # sharp transition for simplicity
thSolver.add_regime_map_model(lamTurb)

# ---------------------------------------------------------------------------
# Numerical discretisation schemes
# ---------------------------------------------------------------------------
thSolver.fvSchemes.divSchemes['default']                          = 'none'
thSolver.fvSchemes.divSchemes['div(phi,alpha)']                   = 'Gauss vanLeer'
thSolver.fvSchemes.divSchemes['div(phir,alpha)']                  = 'Gauss vanLeer'
thSolver.fvSchemes.divSchemes[r'div\(phi.*,U.*\)']               = 'Gauss upwind'
thSolver.fvSchemes.divSchemes['div(alphaRhoPhi,U)']               = 'Gauss upwind'
thSolver.fvSchemes.divSchemes['div(alphaRhoPhiNu,U)']             = 'Gauss linear'
thSolver.fvSchemes.divSchemes['div(alphaRhoPhi,K)']               = 'Gauss upwind'
thSolver.fvSchemes.divSchemes[r'div\(alphaRhoPhi.*,k.*\)']       = 'Gauss upwind'
thSolver.fvSchemes.divSchemes[r'div\(alphaRhoPhi.*,epsilon.*\)'] = 'Gauss upwind'
thSolver.fvSchemes.divSchemes[r'div\(alphaRhoPhi.*,(h|e).*\)']  = 'Gauss upwind'

# ---------------------------------------------------------------------------
# Linear solvers  (GAMG for pressure, smoothSolver for enthalpy)
# ---------------------------------------------------------------------------
fvSolution = ffn.fvSolution()
fvSolution.append("p_rgh",      ffn.fvSolutionSolver(solver='GAMG',         smoother='DIC',          tolerance=1e-7, relTol=0))
fvSolution.append("p_rghFinal", ffn.fvSolutionSolver(solver='GAMG',         smoother='DIC',          tolerance=1e-7, relTol=0))
fvSolution.append("e",          ffn.fvSolutionSolver(solver='smoothSolver',  smoother='symGaussSeidel', tolerance=1e-7, relTol=0, minIter=0))
fvSolution.append("h",          ffn.fvSolutionSolver(solver='smoothSolver',  smoother='symGaussSeidel', tolerance=1e-7, relTol=0, minIter=0))
fvSolution.append('".*"',       ffn.fvSolutionSolver(solver='PBiCGStab',    preconditioner='diagonal', tolerance=1e-7, relTol=0.001))
thSolver.fvSolution = fvSolution

# PIMPLE algorithm settings  (pressure-velocity coupling)
thSolver.pimpleOptions = ffn.PimpleOptions(
    nCorrectors             = 2,     # pressure correctors per outer iteration
    nNonOrthogonalCorrectors= 0,
    nOuterCorrectors        = 6,     # outer PIMPLE iterations (pseudo-time stepping)
    solveEnergy             = True,
    solveFluidMechanics     = True,
    momentumMode            = "faceCentered",
    correctUntilConvergence = None,
    porousInterfaceSharpness= None,
    minMagU                 = None,
    minNOuterCorrectors     = None,
)

# Under-relaxation to stabilise convergence of enthalpy and temperature
thSolver.add_relaxation_on_equation('h', 0.7)
thSolver.add_relaxation_on_field('T', 0.7)

# ---------------------------------------------------------------------------
# Solver list and coupling (no neutronics in this stand-alone case)
# ---------------------------------------------------------------------------
solvers  = ffn.Solvers([thSolver])

for solver in solvers:
    solver.decomposeParDict.numberOfSubdomains = 4
    solver.decomposeParDict.method             = 'simple'
    solver.decomposeParDict.simpleCoeffs['n']  = ffn.Vector(1, 1, 4)

coupling = ffn.Coupling(solvers=solvers)

model = ffn.Model(
    solvers     = solvers,
    coupling    = coupling,
    timeFolders = [timeFolder0]
)

# ---------------------------------------------------------------------------
# Simulation settings
# ---------------------------------------------------------------------------
settings = model.settings
settings.application       = 'GeN-Foam'
settings.endTime           = 100
settings.deltaT            = 1
settings.writeControl      = 'runTime'
settings.writeInterval     = 100
settings.writePrecision    = 7
settings.runTimeModifiable = True
settings.adjustTimeStep    = False
