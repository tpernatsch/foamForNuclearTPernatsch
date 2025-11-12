#==============================================================================*
# Imports

# import numpy as np
import foamForNuclear as ffn
import foamForNuclear.mesh as mesh
# import foamForNuclear.boundaryConditions as bc


#==============================================================================*
# Mesh

thMesh = mesh.BlockMesh(region='fluidRegion')

thMesh.add_pipes_from_json(
    filename="diagram.json",
    firstBlockName="core",
    # firstBlockName="PipeOutlet",
    elbowRadius=0.1,
    isAddBoundaryConditions=True
)


#==============================================================================*
# Time folder

timeFolder0 = ffn.TimeFolder(0)

T = ffn.Field("T", region=thMesh.region)
T.dimensions = ffn.Dimension(default='T')
T.internalField = 600
# T.set_boundary_condition("pipeWall", bc.ZeroGradient())

U = ffn.Field("U", region=thMesh.region)
U.dimensions = ffn.Dimension(default='U')
U.internalField = ffn.Vector(0, 0, 0)
# U.set_boundary_condition("pipeWall", bc.Slip())

p_rgh = ffn.Field("p_rgh", region=thMesh.region)
p_rgh.dimensions = ffn.Dimension(default='p')
p_rgh.internalField = 1e5
# p_rgh.set_boundary_condition("pipeWall", bc.ZeroGradient())

# for face in [
#     "core_outlet", "hotLeg_inlet", "hotLeg_outlet", "heatExchanger_inlet",
#     "heatExchanger_outlet", "pump_inlet", "pump_outlet", "coldLeg_inlet",
#     "coldLeg_outlet", "core_inlet"
# ]:
#     T.set_boundary_condition(face, bc.Cyclic())
#     U.set_boundary_condition(face, bc.Cyclic())
#     p_rgh.set_boundary_condition(face, bc.Cyclic())
    # p_rgh.set_boundary_condition(face, bc.FixedJump(value=p_rgh.internalField, jump=0, rho='"thermo:rho"'))

timeFolder0.append(T)
timeFolder0.append(U)
timeFolder0.append(p_rgh)


#==============================================================================*
# Solver

thSolver = ffn.ThermalHydraulicsSolver(
    region=thMesh.region,
    mesh=thMesh,
    solver='onePhase',
    isSetFvSchemesToDefault=False
)

# thSolver.thermophysicalProperties = ffn.thermophysicalProperty.SodiumConst()

# thSolver.turbulenceProperties.simulationType = 'laminar'

# corePowerModel = ffn.StructureProperty(['core'], volumeFraction=0.5, Dh=0.01)

# corePowerModel.powerModel = ffn.NuclearFuelPin(
#     fuelInnerRadius=0.0012,
#     fuelOuterRadius=0.004715,
#     cladInnerRadius=0.004865,
#     cladOuterRadius=0.005365,
#     fuelMeshSize=16,
#     cladMeshSize=4,
#     fuelRho=10480,
#     fuelCp=250,
#     fuelK=3,
#     cladRho=7500,
#     cladCp=500,
#     cladK=20,
#     fuelT=T.internalField,
#     cladT=T.internalField,
#     gapH=3000,
#     powerDensity=100e6
# )

# corePowerModel.add_passive_structure(
#     volumetricArea=2,
#     rho=7700,
#     Cp=500,
#     T=T.internalField
# )

# heatExchangerModel = ffn.StructureProperty(['heatExchanger'], volumeFraction=0.5, Dh=0.01)
# heatExchangerModel.powerModel = ffn.FixedTemperature(
#     volumetricArea=20,
#     T=T.internalField,
#     Cp=500,
#     rho=7700,
#     powerDensity=0
# )

# pumpModel = ffn.Pump(
#     zones=['pump'],
#     volumeFraction=0.2,
#     Dh=1,
#     momentumSource=ffn.Vector(-100000, 0, 0)
# )

# thSolver.add_structure_property(corePowerModel)
# thSolver.add_structure_property(heatExchangerModel)
# thSolver.add_structure_property(pumpModel)

# thSolver.add_drag_model(
#     ffn.ReynoldsPower(0.687, -0.25, zones=['core', 'heatExchanger'])
# )
# thSolver.add_heat_transfer_model(
#     ffn.NusseltReynoldsPrandtlPower(4.82, 0.0185, 0.827, 0.827, zones=['core', 'heatExchanger'])
# )

# thSolver.pimpleOptions.nOuterCorrectors = 2
# thSolver.pimpleOptions.nCorrectors = 1
# thSolver.pimpleOptions.nNonOrthogonalCorrectors = 0
# thSolver.pimpleOptions.correctPhi = False

# thSolver.fvSchemes.divSchemes['default'] = "Gauss upwind"
# thSolver.fvSchemes.divSchemes['div(alphaRhoPhiNu,U)'] = "Gauss linear"

# thSolver.add_relaxation_on_equation('"U.*"', 1)


#==============================================================================*
# Model

model = ffn.Model(timeFolders=[timeFolder0])
model.settings.application = "GeN-Foam"
model.settings.endTime = 60
model.settings.deltaT = 0.1
model.settings.writeInterval = 1
model.settings.writeControl = "adjustableRunTime"
model.settings.adjustTimeStep = True
model.settings.runTimeModifiable = True
model.settings.maxCo = 0.1
model.settings.maxDeltaT = 1

model.solvers.append(thSolver)


#==============================================================================*
# Run

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=thMesh, show_edges=True)
model.plot_mesh(region=thMesh, show_edges=True, normal='y')
model.plot_mesh(region=thMesh, show_edges=True, normal='z')


#==============================================================================*
