"""

"""
#==============================================================================*
# Imports

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh
import foamForNuclear.porous_medium as porous

import matplotlib.pyplot as plt


#==============================================================================*
# Mesh

thMesh = mesh.BlockMesh(region="fluidRegion")

pPipe1 = thMesh.create_cube("pPipe", 0, 0, 0, 0.1, 0.1, 0.5, nz=10)
heater = thMesh.extrude_top([pPipe1], "heater", dz=1, nz=20)
pHX = thMesh.extrude_top([heater], "pHX", dz=1, nz=20)
pPipe2 = thMesh.extrude_top([pHX], "pPipe", dz=0.5, nz=10)

sPipe1 = thMesh.create_cube("sPipe", 0.2, 0, 0, 0.3, 0.1, 0.5, nz=10)
sPipe2 = thMesh.extrude_top([sPipe1], "sPipe", dz=1, nz=20)
sHX = thMesh.extrude_top([sPipe2], "sHX", dz=1, nz=20)
sPipe3 = thMesh.extrude_top([sHX], "sPipe", dz=0.5, nz=10)

outletP = ffn.mesh.Face("outletP")
outletP.add_sub_face(pPipe2.topFace())

outletS = ffn.mesh.Face("outletS")
outletS.add_sub_face(sPipe3.topFace())

inletP = ffn.mesh.Face("inletP")
inletP.add_sub_face(pPipe1.bottomFace())

inletS = ffn.mesh.Face("inletS")
inletS.add_sub_face(sPipe1.bottomFace())

walls = ffn.mesh.Face("walls", boundaryType="wall")
for face in [pPipe1, heater, pHX, pPipe2, sPipe1, sPipe2, sHX, sPipe3]:
    walls.add_sub_face(face.leftFace())
    walls.add_sub_face(face.rightFace())
    walls.add_sub_face(face.frontFace())
    walls.add_sub_face(face.backFace())

thMesh.add_boundary(outletP)
thMesh.add_boundary(outletS)
thMesh.add_boundary(inletP)
thMesh.add_boundary(inletS)
thMesh.add_boundary(walls)


#==============================================================================*
# Fields

timeFolder0 = ffn.timeFolder.TimeFolder(0)


Tliquid = ffn.fields.Field("T.liquid", region="fluidRegion")
Tliquid.dimensions = ffn.fields.Dimension(default='T')
Tliquid.internalField = 653.15
Tliquid.set_boundary_condition("inletP", bc.FixedValue(Tliquid.internalField))
Tliquid.set_boundary_condition("inletS", bc.FixedValue(Tliquid.internalField))
Tliquid.set_boundary_condition("outletP", bc.ZeroGradient())
Tliquid.set_boundary_condition("outletS", bc.ZeroGradient())
Tliquid.set_boundary_condition("walls", bc.ZeroGradient())

Tvapour = ffn.fields.Field("T.vapour", region="fluidRegion")
Tvapour.dimensions = ffn.fields.Dimension(default='T')
Tvapour.internalField = 653.15
for boundary in thMesh.faces:
    Tvapour.set_boundary_condition(boundary, bc.ZeroGradient())

Uliquid = ffn.fields.Field("U.liquid", region="fluidRegion")
Uliquid.dimensions = ffn.fields.Dimension(default='U')
Uliquid.internalField = ffn.common.Vector(0, 0, 2)
Uliquid.set_boundary_condition("inletP", bc.FixedValue(Uliquid.internalField))
Uliquid.set_boundary_condition("inletS", bc.FixedValue(Uliquid.internalField))
Uliquid.set_boundary_condition("outletP", bc.ZeroGradient())
Uliquid.set_boundary_condition("outletS", bc.ZeroGradient())
Uliquid.set_boundary_condition("walls", bc.Slip())

Uvapour = ffn.fields.Field("U.vapour", region="fluidRegion")
Uvapour.dimensions = ffn.fields.Dimension(default='U')
Uvapour.internalField = ffn.common.Vector(0, 0, 0.01)
Uvapour.set_boundary_condition("inletP", bc.ZeroGradient())
Uvapour.set_boundary_condition("inletS", bc.ZeroGradient())
Uvapour.set_boundary_condition("outletP", bc.ZeroGradient())
Uvapour.set_boundary_condition("outletS", bc.ZeroGradient())
Uvapour.set_boundary_condition("walls", bc.Slip())

alphaVapour = ffn.fields.Field("alpha.vapour", region="fluidRegion")
alphaVapour.dimensions = ffn.fields.Dimension()
alphaVapour.internalField = 0
for boundary in thMesh.faces:
    alphaVapour.set_boundary_condition(boundary, bc.ZeroGradient())

p = ffn.fields.Field("p", region="fluidRegion")
p.dimensions = ffn.fields.Dimension(default="p")
p.internalField = 1e5
for boundary in thMesh.faces:
    p.set_boundary_condition(boundary, bc.Calculated(p.internalField))

p_rgh = ffn.fields.Field("p_rgh", region="fluidRegion")
p_rgh.dimensions = ffn.fields.Dimension(default="p")
p_rgh.internalField = 1e5
p_rgh.set_boundary_condition("inletP", bc.FixedFluxPressure(p_rgh.internalField))
p_rgh.set_boundary_condition("inletS", bc.FixedFluxPressure(p_rgh.internalField))
p_rgh.set_boundary_condition("outletP", bc.FixedValue(p_rgh.internalField))
p_rgh.set_boundary_condition("outletS", bc.FixedValue(p_rgh.internalField))
p_rgh.set_boundary_condition("walls", bc.FixedFluxPressure(p_rgh.internalField))


timeFolder0.append(Tliquid)
timeFolder0.append(Tvapour)
timeFolder0.append(Uliquid)
timeFolder0.append(Uvapour)
timeFolder0.append(alphaVapour)
timeFolder0.append(p)
timeFolder0.append(p_rgh)


#==============================================================================*
# Solvers

thSolver = ffn.solvers.thermal_hydraulics.TwoPhaseThermalHydraulicsSolver(
    region=thMesh.region,
    removeBaffles=False,
    mesh=thMesh
)
# Fluid properties (liquid)
thSolver.fluid1.turbulenceProperties.simulationType = "laminar"
thSolver.fluid1 = ffn.solvers.thermal_hydraulics.Fluid(
    name="liquid",
    thermophysicalProperties=ffn.thermo.SodiumPolynomial(),
    stateOfMatter="liquid",
    thermoResidualAlpha=1e-5,
    dispersedDiameterModel=porous.dispersed_diameter.Constant(value=0.005),
)

# Fluid properties (vapour)
thSolver.fluid2.turbulenceProperties.simulationType = "laminar"
thSolver.fluid2 = ffn.solvers.thermal_hydraulics.Fluid(
    name="vapour",
    thermophysicalProperties=ffn.thermo.SodiumVapourPerfectGas(),
    stateOfMatter="gas",
    thermoResidualAlpha=0.01,
    dispersedDiameterModel=porous.dispersed_diameter.Constant(value=0.005),
)

# Structure properties
pipes = porous.Structure(
    zones=["pPipe", "sPipe", "pHX", "sHX"],
    volumeFraction=0.5,
    Dh=0.005
)

heater = porous.Structure(
    zones=['heater'],
    volumeFraction=0.5,
    Dh=0.005
)
heater.powerModel = porous.power_models.FixedPower(
    volumetricArea=300,
    T=653.15,
    Cp=500,
    rho=7700,
    powerDensity=1.1e9,
)

heatExchanger = porous.HeatExchangerModel(
    name="HeatExchanger1",
    primary="pHX",
    secondary="sHX",
    volumetricArea=250,
    wallConductance=2.5e4
)

thSolver.structures.append(pipes)
thSolver.structures.append(heater)
thSolver.structures.append(heatExchanger)

# Regime map
regimeMap = porous.regime_map.OneParameter(
    name="regimeMap01",
    interpolationMode="quadratic",
    parameter="normalized.alpha.vapour",
    regimeBounds={
        "regime0": [0, 0.79],
        "regime1": [0.99, 1],
    }
)
thSolver.regimeMapModels.append(regimeMap)

# Drag Models
thSolver.fluid_fluid.dragModels.append(
    porous.drag.SchillerNaumann()
)
thSolver.fluid1_structure.dragModels.append(
    porous.drag.ReynoldsPower(coeff=0.316, exp=-0.25, zones=["pPipe", "sPipe", "pHX", "sHX", "heater"])
)

# Heat Transfer Models
thSolver.fluid1_structure.heatTransferModels.append(
    porous.heat_transfer.NusseltReynoldsPrandtlPower(const=4, coeff=0.015, expRe=0.7, expPr=0.7, zones=["heater", "pHX", "sHX"])
)
thSolver.fluid_fluid.heatTransferModels.append(
    porous.heat_transfer.NusseltReynoldsPrandtlPower(const=10, coeff=0, expRe=0, expPr=0, phaseName="liquid")
)
thSolver.fluid_fluid.heatTransferModels.append(
    porous.heat_transfer.NusseltReynoldsPrandtlPower(const=10, coeff=0, expRe=0, expPr=0, phaseName="vapour")
)

# Two Phase Drag Multiplier
thSolver.fluid_fluid.twoPhaseDragMultiplierModel = porous.two_phase_drag_multiplier.LottesFlinn(multiplierFluid="liquid")

# Pair Geometry Models
thSolver.fluid_fluid.pairGeometryModel = porous.pair_geometry.PairGeometryModel(
    dispersionModel=porous.pair_geometry.dispersion.ByRegime(
        # name="RegimeMapFluidFluid",
        regimeMap="regimeMap01",
        regimes=[
            {"name": "regime0", "type": "constant", "dispersedPhase": "vapour"},
            {"name": "regime1", "type": "constant", "dispersedPhase": "liquid"},
        ]
    ),
    interfacialAreaDensityModel=porous.pair_geometry.interfacial_area_density.Spherical()
)

thSolver.fluid1_structure.pairGeometryModel = porous.pair_geometry.PairGeometryModel(
    contactPartitionModel=porous.pair_geometry.contact_partition.ByRegime(
        # name="RegimeMapLiquid",
        regimeMap="regimeMap01",
        regimes=[
            {"name": "regime0", "type": "constant", "value": 1},
            {"name": "regime1", "type": "constant", "value": 0},
        ]
    )
)

# Phase Change Model
thSolver.fluid_fluid.phaseChangeModel = porous.phase_change.HeatDriven(
    mode="conductionLimited",
    correctLatentHeat=False,
    latentHeatModel=porous.phase_change.latent_heat.FinkLeibowitz(adjust=True),
    saturationModel=porous.phase_change.saturation.BrowningPotter()
)

# fvSchemes
thSolver.fvSchemes.laplacianSchemes['default'] = "Gauss linear uncorrected"
thSolver.fvSchemes.divSchemes['default'] = "none"
thSolver.fvSchemes.divSchemes['div(phi,alpha)'] = "Gauss vanLeer"
thSolver.fvSchemes.divSchemes['div(phir,alpha)'] = "Gauss vanLeer"
thSolver.fvSchemes.divSchemes['div\(phi.*,U.*\)'] = "Gauss upwind"
thSolver.fvSchemes.divSchemes['div(alphaRhoPhi,U)'] = "Gauss upwind"
thSolver.fvSchemes.divSchemes['div(alphaRhoPhiNu,U)'] = "Gauss linear"
thSolver.fvSchemes.divSchemes['div(alphaRhoPhi,K)'] = "Gauss upwind"
thSolver.fvSchemes.divSchemes['div\(alphaRhoPhi.*,(h|e).*\)'] = "Gauss upwind"

# fvSolution
thSolution = ffn.numerics.fvSolution()
thSolution.append('"p_rgh.*"', ffn.numerics.fvSolutionSolver(
    solver='GAMG',
    smoother='DIC',
    tolerance=1e-6,
    relTol=0
))
smoothSolver = ffn.numerics.fvSolutionSolver(
    solver='smoothSolver',
    smoother='symGaussSeidel',
    tolerance=1e-6,
    relTol=0,
    minIter=1
)
thSolution.append('"e.*"', smoothSolver)
thSolution.append('"h.*"', smoothSolver)
thSolution.append('alpha', ffn.numerics.fvSolutionSolver(
    solver='MULES',
    adjustSubCycles=True,
    alphaMaxCo=0.25,
    solverPhase="vapour"
))
thSolution.append('".*"', ffn.numerics.fvSolutionSolver(
    solver='PBiCGStab',
    preconditioner='diagonal',
    tolerance=1e-6,
    relTol=0.001
))

thSolver.fvSolution = thSolution


# Pimple Control
thSolver.pimpleOptions.nCorrectors = 12
thSolver.pimpleOptions.nOuterCorrectors = 12
thSolver.pimpleOptions.correctUntilConvergence = True
thSolver.pimpleOptions.partialEliminationMode = "implicit"
thSolver.pimpleOptions.momentumMode = "faceCentered"
thSolver.pimpleOptions.massTransferSafetyFactor = 0.9
thSolver.pimpleOptions.enthalpyStabilizationMode = "source"

thSolver.pimpleOptions.add_residual_control_on_field(
    fieldName="p_rgh",
    tolerance=1e-5,
    relTol=1,
    useFirstPISOInitialResidual=True
)
thSolver.pimpleOptions.add_residual_control_on_field(
    fieldName="h.liquid",
    tolerance=1e-6,
    relTol=0
)

thSolver.add_relaxation_on_equation('".*"', 1)
thSolver.add_relaxation_on_equation("h.vapour", 0.25)
thSolver.add_relaxation_on_equation("h.vapourFinal", 0.25)
thSolver.add_relaxation_on_field("T.interface", 0.25)
thSolver.add_relaxation_on_field("T.interfaceFinal", 0.25)
thSolver.add_relaxation_on_field('"dmdt.*"', 0.125)


#==============================================================================*
# Settings

model = ffn.case.Case()

model.solvers.append(thSolver)
model.timeFolders = [timeFolder0]

settings: ffn.control.ControlDict = model.settings

settings.application = 'GeN-Foam'
settings.endTime = 15
settings.deltaT = 0.0001
settings.writeControl = 'adjustableRunTime'
settings.writeInterval = 1
settings.writePrecision = 6
settings.timePrecision = 6
settings.runTimeModifiable = True
settings.adjustTimeStep = True
settings.maxDeltaT = 0.01
settings.minDeltaT = 0.0001
settings.maxCo = 10
settings.maxCoTwoPhase = 0.5
settings.writeRestartFields = True
settings.writeContinuityErrors = True
settings.includeKineticEnergy = False

funcObjs = {}
for patchName in ["inletP", "inletS", "outletP", "outletS"]:
    funcObjs[patchName] = ffn.functions.SurfaceFieldValue(
        name=f'{patchName}Values',
        fields=["alphaRhoPhi.liquid", "alphaRhoPhi.vapour", "T.liquid", "T.vapour"],
        operation="average",
        region=thMesh.region,
        regionType="patch",
        regionName=patchName,
        log=True,
        writeFields=False
    )

    model.add_function_object(funcObjs[patchName])

print(model)

# Export to OpenFOAM
ffn.allclean()
model.export_to_openfoam()


#==============================================================================*
# Run

ffn.run(case=model, is_preprocessing=True)


#==============================================================================*
# Post-processing

for fieldName, unit in [
    ('T.liquid', 'K'),
    ('T.vapour', 'K'),
    ('alpha.liquid', '-'),
    ('alpha.vapour', '-'),
]:
    model.plot_slice(
        region=thMesh.region,
        time=settings.endTime,
        fieldName=fieldName,
        cmap='RdBu_r',
        unit=unit
    )
    model.plot_animation(
        region=thMesh,
        fieldName=fieldName,
        cmap='RdBu_r',
        normal='y',
        unit=unit
    )


inletP = funcObjs['inletP'].read_from_case(startTime=0)
inletS = funcObjs['inletS'].read_from_case(startTime=0)
outletP = funcObjs['outletP'].read_from_case(startTime=0)
outletS = funcObjs['outletS'].read_from_case(startTime=0)

fig, ax = plt.subplots(figsize=(5, 4), dpi=200)

ax.plot(inletP["Time"], inletP["average(T.liquid)"], label="Inlet Primary")
ax.plot(inletS["Time"], inletS["average(T.liquid)"], label="Inlet Secondary")
ax.plot(outletP["Time"], outletP["average(T.liquid)"], label="Outlet Primary")
ax.plot(outletS["Time"], outletS["average(T.liquid)"], label="Outlet Secondary")
ax.set_xlabel("Time [s]")
ax.set_ylabel("Temperature [K]")
ax.legend()
fig.tight_layout()
fig.savefig("fig_results_T.png")
plt.close()


#==============================================================================*
