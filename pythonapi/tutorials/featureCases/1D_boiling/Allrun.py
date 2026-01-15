"""

"""
#==============================================================================*
# Imports

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh


#==============================================================================*

ffn.allclean()

#==============================================================================*
# Mesh

thMesh = mesh.BlockMesh(region="fluidRegion")

lowBlock = thMesh.create_cube("low", 0, 0, 0, 0.1, 0.01, 0.5, nz=10)
midBlock = thMesh.extrude_top([lowBlock], "mid", dz=1, nz=20)
topBlock = thMesh.extrude_top([midBlock], "top", dz=0.5, nz=10)

outlet = ffn.Face("outlet", boundaryType="wall")
outlet.add_sub_face(topBlock.topFace())

inlet = ffn.Face("inlet", boundaryType="wall")
inlet.add_sub_face(lowBlock.bottomFace())

walls = ffn.Face("walls", boundaryType="empty")
for face in [lowBlock, midBlock, topBlock]:
    walls.add_sub_face(face.leftFace())
    walls.add_sub_face(face.rightFace())
    walls.add_sub_face(face.frontFace())
    walls.add_sub_face(face.backFace())

thMesh.add_boundary(outlet)
thMesh.add_boundary(inlet)
thMesh.add_boundary(walls)


#==============================================================================*
# Fields

timeFolder0 = ffn.TimeFolder(0)


Tliquid = ffn.Field("T.liquid", region="fluidRegion")
Tliquid.dimensions = ffn.Dimension(default='T')
Tliquid.internalField = 652.15
Tliquid.set_boundary_condition("inlet", bc.FixedValue(Tliquid.internalField))
Tliquid.set_boundary_condition("outlet", bc.ZeroGradient())
Tliquid.set_boundary_condition("walls", bc.Empty())

Tstructure = ffn.Field("T.structure", region="fluidRegion")
Tstructure.dimensions = ffn.Dimension(default='T')
Tstructure.internalField = 652.15
Tstructure.set_boundary_condition("inlet", bc.FixedValue(Tliquid.internalField))
Tstructure.set_boundary_condition("outlet", bc.ZeroGradient())
Tstructure.set_boundary_condition("walls", bc.Empty())

Tvapour = ffn.Field("T.vapour", region="fluidRegion")
Tvapour.dimensions = ffn.Dimension(default='T')
Tvapour.internalField = 652.15
Tvapour.set_boundary_condition("inlet", bc.ZeroGradient())
Tvapour.set_boundary_condition("outlet", bc.ZeroGradient())
Tvapour.set_boundary_condition("walls", bc.Empty())

Uliquid = ffn.Field("U.liquid", region="fluidRegion")
Uliquid.dimensions = ffn.Dimension(default='U')
Uliquid.internalField = ffn.Vector(0, 0, 1)
Uliquid.set_boundary_condition("inlet", bc.ZeroGradient())
Uliquid.set_boundary_condition("outlet", bc.ZeroGradient())
Uliquid.set_boundary_condition("walls", bc.Empty())

Uvapour = ffn.Field("U.vapour", region="fluidRegion")
Uvapour.dimensions = ffn.Dimension(default='U')
Uvapour.internalField = ffn.Vector(0, 0, 1)
Uvapour.set_boundary_condition("inlet", bc.ZeroGradient())
Uvapour.set_boundary_condition("outlet", bc.ZeroGradient())
Uvapour.set_boundary_condition("walls", bc.Empty())

alphaLiquid = ffn.Field("alpha.liquid", region="fluidRegion")
alphaLiquid.dimensions = ffn.Dimension()
alphaLiquid.internalField = 1
alphaLiquid.set_boundary_condition("inlet", bc.Calculated(alphaLiquid.internalField))
alphaLiquid.set_boundary_condition("outlet", bc.Calculated(alphaLiquid.internalField))
alphaLiquid.set_boundary_condition("walls", bc.Empty())

alphaVapour = ffn.Field("alpha.vapour", region="fluidRegion")
alphaVapour.dimensions = ffn.Dimension()
alphaVapour.internalField = 0
alphaVapour.set_boundary_condition("inlet", bc.ZeroGradient())
alphaVapour.set_boundary_condition("outlet", bc.ZeroGradient())
alphaVapour.set_boundary_condition("walls", bc.Empty())

p = ffn.Field("p", region="fluidRegion")
p.dimensions = ffn.Dimension(default="p")
p.internalField = 1e5
p.set_boundary_condition("inlet", bc.Calculated(p.internalField))
p.set_boundary_condition("outlet", bc.Calculated(p.internalField))
p.set_boundary_condition("walls", bc.Empty())

p_rgh = ffn.Field("p_rgh", region="fluidRegion")
p_rgh.dimensions = ffn.Dimension(default="p")
p_rgh.internalField = 1e5
p_rgh.set_boundary_condition("inlet", bc.FixedValue(1.4e5))
p_rgh.set_boundary_condition("outlet", bc.FixedValue(p_rgh.internalField))
p_rgh.set_boundary_condition("walls", bc.Empty())


timeFolder0.append(Tliquid)
timeFolder0.append(Tstructure)
timeFolder0.append(Tvapour)
timeFolder0.append(Uliquid)
timeFolder0.append(Uvapour)
timeFolder0.append(alphaLiquid)
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
phaseProperties.residualKd = 10

# Structure properties
lowTopStruct = ffn.StructureProperty(
    zones=["low", "top"],
    volumeFraction=0.5,
    Dh=0.01
)

midStruct = ffn.StructureProperty(zones=['mid'], volumeFraction=0.5, Dh=0.005)
midStruct.powerModel = ffn.FixedPower(
    volumetricArea=100,
    T=652.15,
    Cp=500,
    rho=7700,
    powerDensity=4.5e8,
)

phaseProperties.structureProperties.append(lowTopStruct)
phaseProperties.structureProperties.append(midStruct)

# Fluid Properties
liquidProperty = ffn.FluidProperty(
    stateOfMatter="liquid",
    dispersedDiameterModel=ffn.ConstantDispersedDiameterModel(0.01)
)
vapourProperty = ffn.FluidProperty(
    stateOfMatter="gas",
    thermoResidualAlpha=0.1,
    dispersedDiameterModel=ffn.ConstantDispersedDiameterModel(0.01)
)
phaseProperties.add_fluid_properties(liquidProperty)
phaseProperties.add_fluid_properties(vapourProperty)

# Power Off Criterion
phaseProperties.structureProperties.powerOffCriterionModel = ffn.TimerPowerOffCriterionModel(
    time=7
)

# Regime map
regimeMap = ffn.RegimeMapOneParameter(
    name="slugMist",
    parameter="normalized.alpha.vapour",
    regimeBounds=ffn.OpenFOAMDict({
        "slug": ffn.List([0, 0.85]),
        "mist": ffn.List([0.95, 1]),
    })
)
phaseProperties.regimeMapModels.append(regimeMap)

# Drag Models
phaseProperties.dragModels.append(
    ffn.SchillerNaumann(zones=["liquid.vapour"]),
    interactionType="fluid-fluid"
)
phaseProperties.dragModels.append(
    ffn.ReynoldsPower(coeff=2, exp=-0.125, zones=["low", "mid", "top"]),
    phaseName="liquid"
)

# Heat Transfer Models
phaseProperties.heatTransferModels.append(
    ffn.SuperpositionNucleateBoiling(
        zones=['mid'],
        forcedConvectionModel=ffn.NusseltReynoldsPrandtlPower(
            const=7.48467, coeff=0.02994, expRe=0.77, expPr=0.77
        ),
        poolBoilingModel=ffn.Shah(useExplicitHeatFlux=False),
        flowEnhancementFactorModel=ffn.CobraTfFlowEnhancementFactor(),
        suppressionFactorModel=ffn.CobraTfSuppressionFactor()
    ),
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
phaseProperties.twoPhaseDragMultiplierModel = ffn.Kaiser88("liquid")

# Pair Geometry Models
phaseProperties.pairGeometryModels.append(
    ffn.PairGeometryModel(
        dispersionModel=ffn.ByRegimeDispersionModel(
            regimeMap="slugMist",
            regimes={
                "slug": ffn.OpenFOAMDict({"type": "constant", "dispersedPhase": "vapour"}),
                "mist": ffn.OpenFOAMDict({"type": "constant", "dispersedPhase": "liquid"}),
            }
        ),
        interfacialAreaDensityModel=ffn.SphericalInterfacialAreaDensityModel()
    ),
    interactionType="fluid-fluid"
)
phaseProperties.pairGeometryModels.append(
    ffn.PairGeometryModel(
        contactPartitionModel=ffn.ByRegimeContactPartitionModel(
            regimeMap="slugMist",
            regimes={
                "slug": ffn.OpenFOAMDict({"type": "constant", "value": 1.0}),
                "mist": ffn.OpenFOAMDict({"type": "constant", "value": 0.1}),
            }
        )
    ),
    phaseName="liquid"
)

# Phase Change Model
phaseProperties.phaseChangeModel = ffn.HeatDriven(
    mode="conductionLimited",
    correctLatentHeat=False,
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
    alphaMaxCo=0.25
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
thSolver.pimpleOptions.nOuterCorrectors = 24
thSolver.pimpleOptions.partialEliminationMode = "implicit"
thSolver.pimpleOptions.momentumMode = "faceCentered"
thSolver.pimpleOptions.oscillationLimiterFraction = 0
thSolver.pimpleOptions.correctUntilConvergence = True
thSolver.pimpleOptions.massTransferSafetyFactor = 0.1
thSolver.pimpleOptions.maxTInterfaceDdt = 10000
thSolver.pimpleOptions.enthalpyStabilizationMode = "source"

thSolver.pimpleOptions.add_residual_control_on_field(
    fieldName="p_rgh",
    tolerance=1e-5,
    relTol=0,
    useFirstPISOInitialResidual=True
)
thSolver.pimpleOptions.add_residual_control_on_field(
    fieldName="h.liquid",
    tolerance=1e-5,
    relTol=0
)

thSolver.add_relaxation_on_equation('".*"', 1)
thSolver.add_relaxation_on_field('"alpha.*"', 0.5)
thSolver.add_relaxation_on_field('"dmdt.*"', 0.125)
thSolver.add_relaxation_on_field("T.interface", 0.125)
thSolver.add_relaxation_on_field("T.interfaceFinal", 0.125)
thSolver.add_relaxation_on_field("h.vapour", 0.0125)
thSolver.add_relaxation_on_field("h.vapourFinal", 0.0125)


#==============================================================================*
# Settings

model = ffn.Model()

model.solvers.append(thSolver)
model.timeFolders = [timeFolder0]

settings: ffn.ControlDict = model.settings

settings.application = 'GeN-Foam'
settings.endTime = 12
settings.deltaT = 0.001
settings.writeControl = 'adjustableRunTime'
settings.writeInterval = 0.1
settings.writePrecision = 6
settings.timePrecision = 6
settings.runTimeModifiable = True
settings.adjustTimeStep = True
settings.maxDeltaT = 0.01
settings.maxCo = 2
settings.maxCoTwoPhase = 0.25
settings.writeContinuityErrors = True

print(model)

# Export to OpenFOAM
model.export_to_openfoam()


#==============================================================================*
# Run

ffn.run(model=model, is_preprocessing=True)


#==============================================================================*
# Post-processing

model.plot_slice(
    region=thMesh.region,
    time=settings.endTime,
    fieldName='T.liquid',
    cmap='RdBu_r',
    unit='K'
)
model.plot_slice(
    region=thMesh.region,
    time=settings.endTime,
    fieldName='T.vapour',
    cmap='RdBu_r',
    unit='K'
)
model.plot_slice(
    region=thMesh.region,
    time=settings.endTime,
    fieldName='alpha.liquid',
    cmap='RdBu_r',
    unit='K'
)
model.plot_slice(
    region=thMesh.region,
    time=settings.endTime,
    fieldName='alpha.vapour',
    cmap='RdBu_r',
    unit='K'
)
model.plot_animation(
    region=thMesh,
    fieldName='alpha.liquid',
    normal='y',
    cmap='RdBu_r',
    unit='K'
)
model.plot_animation(
    region=thMesh,
    fieldName='T.liquid',
    normal='y',
    cmap='RdBu_r',
    unit='K'
)
model.plot_animation(
    region=thMesh,
    fieldName='T.vapour',
    normal='y',
    cmap='RdBu_r',
    unit='K'
)
model.plot_residuals(
    parameters=["h.liquid", "h.vapour", "p_rgh"],
    title="Thermal-hydraulics"
)


#==============================================================================*
