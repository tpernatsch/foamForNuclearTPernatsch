"""

"""
#==============================================================================*
# Imports

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh

import matplotlib.pyplot as plt

#==============================================================================*
# Mesh

thMesh = mesh.BlockMesh(region="fluidRegion")

pPipe1 = thMesh.createCube("pPipe", 0, 0, 0, 0.1, 0.1, 0.5, nz=10)
heater = thMesh.extrudeTop([pPipe1], "heater", dz=1, nz=20)
pHX = thMesh.extrudeTop([heater], "pHX", dz=1, nz=20)
pPipe2 = thMesh.extrudeTop([pHX], "pPipe", dz=0.5, nz=10)

sPipe1 = thMesh.createCube("sPipe", 0.2, 0, 0, 0.3, 0.1, 0.5, nz=10)
sPipe2 = thMesh.extrudeTop([sPipe1], "sPipe", dz=1, nz=20)
sHX = thMesh.extrudeTop([sPipe2], "sHX", dz=1, nz=20)
sPipe3 = thMesh.extrudeTop([sHX], "sPipe", dz=0.5, nz=10)

outletP = ffn.Face("outletP")
outletP.addSubFace(pPipe2.topFace())

outletS = ffn.Face("outletS")
outletS.addSubFace(sPipe3.topFace())

inletP = ffn.Face("inletP")
inletP.addSubFace(pPipe1.bottomFace())

inletS = ffn.Face("inletS")
inletS.addSubFace(sPipe1.bottomFace())

walls = ffn.Face("walls", boundaryType="wall")
for face in [pPipe1, heater, pHX, pPipe2, sPipe1, sPipe2, sHX, sPipe3]:
    walls.addSubFace(face.leftFace())
    walls.addSubFace(face.rightFace())
    walls.addSubFace(face.frontFace())
    walls.addSubFace(face.backFace())

thMesh.addBoundary(outletP)
thMesh.addBoundary(outletS)
thMesh.addBoundary(inletP)
thMesh.addBoundary(inletS)
thMesh.addBoundary(walls)


#==============================================================================*
# Fields

timeFolder0 = ffn.TimeFolder(0)


Tliquid = ffn.Field("T.liquid", region="fluidRegion")
Tliquid.dimensions = ffn.Dimension(default='T')
Tliquid.internalField = 653.15
Tliquid.set_boundary_condition("inletP", bc.FixedValue(Tliquid.internalField))
Tliquid.set_boundary_condition("inletS", bc.FixedValue(Tliquid.internalField))
Tliquid.set_boundary_condition("outletP", bc.ZeroGradient())
Tliquid.set_boundary_condition("outletS", bc.ZeroGradient())
Tliquid.set_boundary_condition("walls", bc.ZeroGradient())

Tvapour = ffn.Field("T.vapour", region="fluidRegion")
Tvapour.dimensions = ffn.Dimension(default='T')
Tvapour.internalField = 653.15
for boundary in thMesh.faces:
    Tvapour.set_boundary_condition(boundary, bc.ZeroGradient())

Uliquid = ffn.Field("U.liquid", region="fluidRegion")
Uliquid.dimensions = ffn.Dimension(default='U')
Uliquid.internalField = ffn.Vector(0, 0, 2)
Uliquid.set_boundary_condition("inletP", bc.FixedValue(Uliquid.internalField))
Uliquid.set_boundary_condition("inletS", bc.FixedValue(Uliquid.internalField))
Uliquid.set_boundary_condition("outletP", bc.ZeroGradient())
Uliquid.set_boundary_condition("outletS", bc.ZeroGradient())
Uliquid.set_boundary_condition("walls", bc.Slip())

Uvapour = ffn.Field("U.vapour", region="fluidRegion")
Uvapour.dimensions = ffn.Dimension(default='U')
Uvapour.internalField = ffn.Vector(0, 0, 0.01)
Uvapour.set_boundary_condition("inletP", bc.ZeroGradient())
Uvapour.set_boundary_condition("inletS", bc.ZeroGradient())
Uvapour.set_boundary_condition("outletP", bc.ZeroGradient())
Uvapour.set_boundary_condition("outletS", bc.ZeroGradient())
Uvapour.set_boundary_condition("walls", bc.Slip())

alphaVapour = ffn.Field("alpha.vapour", region="fluidRegion")
alphaVapour.dimensions = ffn.Dimension()
alphaVapour.internalField = 0
for boundary in thMesh.faces:
    alphaVapour.set_boundary_condition(boundary, bc.ZeroGradient())

p = ffn.Field("p", region="fluidRegion")
p.dimensions = ffn.Dimension(default="p")
p.internalField = 1e5
for boundary in thMesh.faces:
    p.set_boundary_condition(boundary, bc.Calculated(p.internalField))

p_rgh = ffn.Field("p_rgh", region="fluidRegion")
p_rgh.dimensions = ffn.Dimension(default="p")
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

thSolver = ffn.ThermalHydraulicsSolver(
    region="fluidRegion",
    solver="twoPhase",
    removeBaffles=False,
    mesh=thMesh
)
thSolver.turbulenceProperties.simulationType = "laminar"
thSolver.turbulencePropertiesSecondFluid.simulationType = "laminar"
thSolver.thermophysicalProperties = ffn.thermophysicalProperty.SodiumPolynomial()
thSolver.thermophysicalPropertiesSecondFluid = ffn.thermophysicalProperty.SodiumVapourPerfectGas()

phaseProperties = thSolver.phaseProperties
phaseProperties.phaseNames = ['liquid', 'vapour']

# Structure properties
pipes = ffn.StructureProperty(
    zones=["pPipe", "sPipe", "pHX", "sHX"],
    volumeFraction=0.5,
    Dh=0.005
)

heater = ffn.StructureProperty(zones=['heater'], volumeFraction=0.5, Dh=0.005)
heater.powerModel = ffn.FixedPower(
    volumetricArea=300,
    T=653.15,
    Cp=500,
    rho=7700,
    powerDensity=1.1e9,
)

heatExchanger = ffn.HeatExchangerModel(
    name="HeatExchanger1",
    primary="pHX",
    secondary="sHX",
    volumetricArea=250,
    wallConductance=2.5e4
)

phaseProperties.structureProperties.append(pipes)
phaseProperties.structureProperties.append(heater)
phaseProperties.structureProperties.append(heatExchanger)

# Fluid Properties
liquidProperty = ffn.FluidProperty(
    thermoResidualAlpha=1e-5,
    stateOfMatter="liquid",
    dispersedDiameterModel=ffn.ConstantDispersedDiameterModel(0.005)
)
vapourProperty = ffn.FluidProperty(
    thermoResidualAlpha=0.01,
    stateOfMatter="gas",
    dispersedDiameterModel=ffn.ConstantDispersedDiameterModel(0.005)
)
phaseProperties.add_fluid_properties(liquidProperty)
phaseProperties.add_fluid_properties(vapourProperty)

# Regime map
regimeMap = ffn.RegimeMapOneParameter(
    name="regimeMap01",
    interpolationMode="quadratic",
    parameter="normalized.alpha.vapour",
    regimeBounds=ffn.OpenFOAMDict({
        "regime0": ffn.List([0, 0.79]),
        "regime1": ffn.List([0.99, 1]),
    })
)
phaseProperties.regimeMapModels.append(regimeMap)

# Drag Models
phaseProperties.dragModels.append(
    ffn.SchillerNaumann(zones=["liquid.vapour"]),
    interactionType="fluid-fluid"
)
phaseProperties.dragModels.append(
    ffn.ReynoldsPower(coeff=0.316, exp=-0.25, zones=["pPipe", "sPipe", "pHX", "sHX", "heater"]),
    phaseName="liquid"
)

# Heat Transfer Models
phaseProperties.heatTransferModels.append(
    ffn.NusseltReynoldsPrandtlPower(const=4, coeff=0.015, expRe=0.7, expPr=0.7, zones=["heater", "pHX", "sHX"]),
    interactionType="fluid-structure",
    phaseName="liquid"
)
phaseProperties.heatTransferModels.append(
    ffn.NusseltReynoldsPrandtlPowerFluidFluid(const=10, coeff=0, expRe=0, expPr=0, phaseName="liquid"),
    interactionType="fluid-fluid"
)
phaseProperties.heatTransferModels.append(
    ffn.NusseltReynoldsPrandtlPowerFluidFluid(const=10, coeff=0, expRe=0, expPr=0, phaseName="vapour"),
    interactionType="fluid-fluid"
)

# Two Phase Drag Multiplier
phaseProperties.twoPhaseDragMultiplierModel = ffn.LottesFlinn("liquid")

# Pair Geometry Models
phaseProperties.pairGeometryModels.append(
    ffn.PairGeometryModel(
        dispersionModel=ffn.ByRegimeDispersionModel(
            regimeMap="regimeMap01",
            regimes={
                "regime0": ffn.OpenFOAMDict({"type": "constant", "dispersedPhase": "vapour"}),
                "regime1": ffn.OpenFOAMDict({"type": "constant", "dispersedPhase": "liquid"}),
            }
        ),
        interfacialAreaDensityModel=ffn.SphericalInterfacialAreaDensityModel()
    ),
    interactionType="fluid-fluid"
)
phaseProperties.pairGeometryModels.append(
    ffn.PairGeometryModel(
        contactPartitionModel=ffn.ByRegimeContactPartitionModel(
            regimeMap="regimeMap01",
            regimes={
                "regime0": ffn.OpenFOAMDict({"type": "constant", "value": 1}),
                "regime1": ffn.OpenFOAMDict({"type": "constant", "value": 0}),
            }
        )
    ),
    phaseName="liquid"
)

# Phase Change Model
phaseProperties.phaseChangeModel = ffn.HeatDriven(
    mode="conductionLimited",
    latentHeatModel=ffn.FinkLeibowitzLatentHeat(adjust=True),
    saturationModel=ffn.BrowningPotterSaturationModel()
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
thSolution = ffn.fvSolution()
thSolution.append('"p_rgh.*"', ffn.fvSolutionSolver(
    solver='GAMG',
    smoother='DIC',
    tolerance=1e-6,
    relTol=0
))
smoothSolver = ffn.fvSolutionSolver(
    solver='smoothSolver',
    smoother='symGaussSeidel',
    tolerance=1e-6,
    relTol=0,
    minIter=1
)
thSolution.append('"e.*"', smoothSolver)
thSolution.append('"h.*"', smoothSolver)
thSolution.append('alpha', ffn.fvSolutionSolver(
    solver='MULES',
    adjustSubCycles=True,
    alphaMaxCo=0.25,
    solverPhase="vapour"
))
thSolution.append('".*"', ffn.fvSolutionSolver(
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

model = ffn.Model()

model.solvers.append(thSolver)
model.timeFolders = [timeFolder0]

settings: ffn.ControlDict = model.settings

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
    funcObjs[patchName] = ffn.SurfaceFieldValue(
        name=f'{patchName}Values',
        fields=["alphaRhoPhi.liquid", "alphaRhoPhi.vapour", "T.liquid", "T.vapour"],
        operation="average",
        region=thMesh.region,
        regionType="patch",
        regionName=patchName,
        log=True,
        writeFields=False
    )

    settings.add_function_object(funcObjs[patchName])

print(model)

# Export to OpenFOAM
ffn.allclean()
model.export_to_openfoam()


#==============================================================================*
# Run

ffn.run(model=model, is_preprocessing=True)


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
