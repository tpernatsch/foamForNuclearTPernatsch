"""
caseSetup.py — Coupled TH + neutronics solver configuration
============================================================

This file configures both physics solvers and — crucially — defines how
they exchange information with each other.

MULTI-PHYSICS COUPLING
-----------------------
The neutronics and thermal-hydraulics solvers are coupled through field
transfers executed at the end of each global iteration:

  neutronics → TH  (power deposition)
  ┌─────────────────────────────────────────────────────────────┐
  │ powerDensity         → powerDensityStructure               │
  │   nuclear heat generation rate [W/m³] mapped onto the      │
  │   fluid mesh and used to heat the coolant                   │
  │ secondaryPowerDensity → powerDensityLiquid                  │
  │   gamma/neutron heating directly in the coolant            │
  └─────────────────────────────────────────────────────────────┘

  TH → neutronics  (temperature and density feedback)
  ┌─────────────────────────────────────────────────────────────┐
  │ T                    → TCool     coolant temperature        │
  │ thermo:rho           → rhoCool   coolant density           │
  │ T.fuelAvForNeutronics → TFuel    volume-averaged fuel temp  │
  │ T.cladAvForNeutronics → TClad    volume-averaged clad temp  │
  │ T.passiveStructure   → TStructMech  passive structure temp  │
  └─────────────────────────────────────────────────────────────┘

These temperature and density fields change the neutron cross sections
(Doppler broadening, coolant void effect, structural expansion effects),
which in turn changes the power distribution → fully self-consistent solution.
"""

import foamForNuclear as ffn

from mesh               import nMesh, thMesh
from boundaryConditions import timeFolder0, INLET_TEMPERATURE

# ===========================================================================
# Neutronics solver
# ===========================================================================
neutronicsSolver = ffn.NeutronicsSolver(
    solver            = "diffusionNeutronics",
    mesh              = nMesh,
    region            = nMesh.region,
    power             = 8e8,
    keff              = 0.9388902,
    isMeshDeformation = False,   # no structural deformation in this case (no TM solver)
)
neutronicsSolver.neutronTransportOptions.maxNeutronIterations = 50

# ===========================================================================
# Thermal-hydraulics solver
# ===========================================================================
thSolver = ffn.ThermalHydraulicsSolver(
    region            = "fluidRegion",
    solver            = "onePhase",
    removeBaffles     = True,
    mesh              = thMesh,
    isSetFvSolutionToDefault = False
)

thSolver.thermophysicalProperties = ffn.thermophysicalProperty.SodiumPolynomial()
thSolver.thermophysicalProperties.Cp    = ffn.Polynome(1646.97, -0.831587, 4.31182e-04)
thSolver.thermophysicalProperties.rho   = ffn.Polynome(860)
thSolver.thermophysicalProperties.mu    = ffn.Polynome(2e-4)
thSolver.thermophysicalProperties.kappa = ffn.Polynome(60)

thSolver.turbulenceProperties.simulationType = 'laminar'

STRUCTURE_ZONES = [
    "diagrid", "upperGasPlenum", "upperReflector",
    "lowerGasPlenum", "lowerReflector",
    "radialReflector", "follower", "controlRod"
]
CORE_ZONES = ["innerCore", "outerCore"]
ALL_ZONES  = CORE_ZONES + STRUCTURE_ZONES

passiveStructures = ffn.StructureProperty(STRUCTURE_ZONES, volumeFraction=0.718520968, Dh=0.00365)
passiveStructures.add_passive_structure(volumetricArea=5, rhoCp=4.8e6, T=INLET_TEMPERATURE)
thSolver.add_structure_property(passiveStructures)

core = ffn.StructureProperty(zones=CORE_ZONES, volumeFraction=0.718520968, Dh=0.00365)
core.add_power_model(ffn.NuclearFuelPin(
    fuelInnerRadius=0.0012, fuelOuterRadius=0.004715,
    cladInnerRadius=0.004865, cladOuterRadius=0.005365,
    fuelMeshSize=30, cladMeshSize=5,
    fuelRho=10480, fuelCp=250, fuelK=3,
    cladRho=7500,  cladCp=500, cladK=20,
    gapH=3000, fuelT=INLET_TEMPERATURE, cladT=INLET_TEMPERATURE
))
core.add_passive_structure(volumetricArea=5, rhoCp=4.8e6, T=INLET_TEMPERATURE)
thSolver.add_structure_property(core)

thSolver.add_drag_model(ffn.ReynoldsPower(coeff=0.687, exp=-0.25, zones=ALL_ZONES))

heatTransfer = ffn.HeatTransferByRegime("lamTurb", zones=ALL_ZONES)
heatTransfer.add_regime(ffn.NusseltReynoldsPrandtlPower(const=4, coeff=0, expRe=0, expPr=0, zones=['laminar']))
heatTransfer.add_regime(ffn.NusseltReynoldsPrandtlPower(const=4.82, coeff=0.0185, expRe=0.827, expPr=0.827, zones=['turbulent']))
thSolver.add_heat_transfer_model(heatTransfer)

lamTurb = ffn.RegimeMapOneParameter(name="lamTurb", parameter="Re")
lamTurb.add_regime("laminar", 0, 1000)
lamTurb.add_regime("turbulent", 2300, 2301)
thSolver.add_regime_map_model(lamTurb)

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

fvSolution = ffn.fvSolution()
fvSolution.append("p_rgh",      ffn.fvSolutionSolver(solver='GAMG',        smoother='DIC',           tolerance=1e-7, relTol=0))
fvSolution.append("p_rghFinal", ffn.fvSolutionSolver(solver='GAMG',        smoother='DIC',           tolerance=1e-7, relTol=0))
fvSolution.append("e",          ffn.fvSolutionSolver(solver='smoothSolver', smoother='symGaussSeidel', tolerance=1e-7, relTol=0, minIter=0))
fvSolution.append("h",          ffn.fvSolutionSolver(solver='smoothSolver', smoother='symGaussSeidel', tolerance=1e-7, relTol=0, minIter=0))
fvSolution.append('".*"',       ffn.fvSolutionSolver(solver='PBiCGStab',   preconditioner='diagonal', tolerance=1e-7, relTol=0.001))
thSolver.fvSolution = fvSolution

thSolver.pimpleOptions = ffn.PimpleOptions(
    nCorrectors=2, nNonOrthogonalCorrectors=0, nOuterCorrectors=6,
    solveEnergy=True, solveFluidMechanics=True, momentumMode="faceCentered",
    correctUntilConvergence=None, porousInterfaceSharpness=None,
    minMagU=None, minNOuterCorrectors=None,
)
thSolver.add_relaxation_on_equation('h', 0.7)
thSolver.add_relaxation_on_field('T', 0.7)

# ===========================================================================
# Solver list
# ===========================================================================
solvers = ffn.Solvers([thSolver, neutronicsSolver])

for solver in solvers:
    solver.decomposeParDict.numberOfSubdomains = 4
    solver.decomposeParDict.method             = 'simple'
    solver.decomposeParDict.simpleCoeffs['n']  = ffn.Vector(1, 1, 4)

# ===========================================================================
# Multi-physics coupling — field transfers between solvers
# ===========================================================================
coupling = ffn.Coupling(solvers=solvers)

# Power from neutronics → deposited as heat in the TH coolant
coupling.add_field_transfer(neutronicsSolver, thSolver, "powerDensity",          "powerDensityStructure")
coupling.add_field_transfer(neutronicsSolver, thSolver, "secondaryPowerDensity", "powerDensityLiquid")

# Temperature and density feedback from TH → neutronics cross sections
coupling.add_field_transfer(thSolver, neutronicsSolver, "T",                     "TCool")
coupling.add_field_transfer(thSolver, neutronicsSolver, "thermo:rho",            "rhoCool")
coupling.add_field_transfer(thSolver, neutronicsSolver, "T.fuelAvForNeutronics", "TFuel")
coupling.add_field_transfer(thSolver, neutronicsSolver, "T.cladAvForNeutronics", "TClad")
coupling.add_field_transfer(thSolver, neutronicsSolver, "T.passiveStructure",    "TStructMech")

# ===========================================================================
# Model and settings
# ===========================================================================
model = ffn.Model(
    solvers     = solvers,
    coupling    = coupling,
    timeFolders = [timeFolder0]
)

settings = model.settings
settings.application       = 'GeN-Foam'
settings.endTime           = 100
settings.deltaT            = 1
settings.writeControl      = 'runTime'
settings.writeInterval     = 100
settings.writePrecision    = 7
settings.runTimeModifiable = True
settings.adjustTimeStep    = False
