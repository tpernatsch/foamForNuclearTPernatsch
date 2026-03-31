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

ffn.allclean()

#==============================================================================*
# Mesh

thMesh = mesh.BlockMesh(region="fluidRegion")

pPipe1 = thMesh.create_cube("pPipe", 0, 0, 0, 0.1, 0.1, 0.5, nz=10)
heater = thMesh.extrude_top(pPipe1, "heater", dz=1, nz=20)
pHX = thMesh.extrude_top(heater, "pHX", dz=1, nz=20)
pPipe2 = thMesh.extrude_top(pHX, "pPipe", dz=0.5, nz=10)

sPipe1 = thMesh.create_cube("sPipe", 0.2, 0, 0, 0.3, 0.1, 0.5, nz=10)
sPipe2 = thMesh.extrude_top(sPipe1, "sPipe", dz=1, nz=20)
sHX = thMesh.extrude_top(sPipe2, "sHX", dz=1, nz=20)
sPipe3 = thMesh.extrude_top(sHX, "sPipe", dz=0.5, nz=10)

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

timeFolder0 = ffn.TimeFolder(0)


T = ffn.fields.Field("T", region="fluidRegion")
T.dimensions = ffn.fields.Dimension(default='T')
T.internalField = 653.15
T.set_boundary_condition("inletP", bc.FixedValue(T.internalField))
T.set_boundary_condition("inletS", bc.FixedValue(T.internalField))
T.set_boundary_condition("outletP", bc.ZeroGradient())
T.set_boundary_condition("outletS", bc.ZeroGradient())
T.set_boundary_condition("walls", bc.ZeroGradient())

U = ffn.fields.Field("U", region="fluidRegion")
U.dimensions = ffn.fields.Dimension(default='U')
U.internalField = ffn.common.Vector(0, 0, 2)
U.set_boundary_condition("inletP", bc.FixedValue(U.internalField))
U.set_boundary_condition("inletS", bc.FixedValue(U.internalField))
U.set_boundary_condition("outletP", bc.ZeroGradient())
U.set_boundary_condition("outletS", bc.ZeroGradient())
U.set_boundary_condition("walls", bc.Slip())

p = ffn.fields.Field("p", region="fluidRegion")
p.dimensions = ffn.fields.Dimension(default="p")
p.internalField = 1e5
p.set_boundary_condition("inletP", bc.Calculated(p.internalField))
p.set_boundary_condition("inletS", bc.Calculated(p.internalField))
p.set_boundary_condition("outletP", bc.Calculated(p.internalField))
p.set_boundary_condition("outletS", bc.Calculated(p.internalField))
p.set_boundary_condition("walls", bc.Calculated(p.internalField))

p_rgh = ffn.fields.Field("p_rgh", region="fluidRegion")
p_rgh.dimensions = ffn.fields.Dimension(default="p")
p_rgh.internalField = 1e5
p_rgh.set_boundary_condition("inletP", bc.FixedFluxPressure(p_rgh.internalField))
p_rgh.set_boundary_condition("inletS", bc.FixedFluxPressure(p_rgh.internalField))
p_rgh.set_boundary_condition("outletP", bc.FixedFluxPressure(p_rgh.internalField))
p_rgh.set_boundary_condition("outletS", bc.FixedFluxPressure(p_rgh.internalField))
p_rgh.set_boundary_condition("walls", bc.FixedFluxPressure(p_rgh.internalField))


timeFolder0.append(T)
timeFolder0.append(U)
timeFolder0.append(p)
timeFolder0.append(p_rgh)


#==============================================================================*
# Solvers

thSolver = ffn.solvers.thermal_hydraulics.OnePhase(
    region=thMesh.region,
    removeBaffles=False,
    mesh=thMesh
)
thSolver.fluid.turbulenceProperties.simulationType = "laminar"
# GeN-Foam crashes with Polynomial thermophysical properties
# thSolver.fluid.thermophysicalProperties = ffn.thermo.SodiumPolynomial()
thSolver.fluid.thermophysicalProperties = ffn.thermo.SodiumConst()


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

# pump = ffn.Pump(
#     zones=['pump'],
#     Dh=0.5,
#     volumeFraction=0,
#     momentumSource=ffn.Vector(0, 0, -166000),
#     momentumSourceTimeProfile=ffn.TimeProfile(
#         type='table',
#         table=[
#             (   0,   1   ),
#             (   1,   0.8490104619    ),
#             (   2,   0.7126676467    ),
#             (   3,   0.5704819773    ),
#             (   4,   0.4544049355    ),
#             (   5,   0.4067212859    ),
#         ]
#     )
# )


thSolver.structures.append(pipes)
thSolver.structures.append(heater)
thSolver.structures.append(heatExchanger)

thSolver.fluid_structure.dragModels.append(
    porous.drag.ReynoldsPower(
        coeff=0.316,
        exp=-0.25,
        zones=["pPipe", "sPipe", "pHX", "sHX", "heater"]
    )
)
thSolver.fluid_structure.heatTransferModels.append(
    porous.heat_transfer.NusseltReynoldsPrandtlPower(
        const=4,
        coeff=0.02,
        expRe=0.8,
        expPr=0.8,
        zones=["heater", "pHX", "sHX"]
    )
)

thSolver.fvSchemes.laplacianSchemes['default'] = "Gauss linear uncorrected"
thSolver.fvSchemes.divSchemes['default'] = "none"
thSolver.fvSchemes.divSchemes['div(phi,alpha)'] = "Gauss vanLeer"
thSolver.fvSchemes.divSchemes['div(phir,alpha)'] = "Gauss vanLeer"
thSolver.fvSchemes.divSchemes['div\(phi.*,U.*\)'] = "Gauss upwind"
thSolver.fvSchemes.divSchemes['div(alphaRhoPhi,U)'] = "Gauss upwind"
thSolver.fvSchemes.divSchemes['div(alphaRhoPhiNu,U)'] = "Gauss linear"
thSolver.fvSchemes.divSchemes['div(alphaRhoPhi,K)'] = "Gauss upwind"
thSolver.fvSchemes.divSchemes['div\(alphaRhoPhi.*,(h|e).*\)'] = "Gauss upwind"

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
thSolution.append('".*"', ffn.numerics.fvSolutionSolver(
    solver='PBiCGStab',
    preconditioner='diagonal',
    tolerance=1e-6,
    relTol=0.001
))

thSolver.fvSolution = thSolution


thSolver.pimpleOptions.nCorrectors = 6
thSolver.pimpleOptions.nOuterCorrectors = 6
thSolver.pimpleOptions.correctUntilConvergence = True
thSolver.pimpleOptions.momentumMode = "faceCentered"
thSolver.pimpleOptions.continuityErrorCompensationMode = "Su"
thSolver.pimpleOptions.continuityErrorScaleFactor = 0.9

thSolver.add_relaxation_on_equation('".*"', 1)


#==============================================================================*
# Function objects

massFlow = ffn.functions.MassFlow(
    'inletMassFlow',
    log=True,
    writeFields=False,
    writeControl="adjustableRunTime",
    writeInterval=0.01,
    region=thMesh.region,
    regionName='inletP'
)
TBulk = ffn.functions.TBulk(
    'outletTBulkP',
    log=True,
    writeFields=False,
    writeControl="adjustableRunTime",
    writeInterval=0.01,
    region=thMesh.region,
    regionName='outletP'
)

#==============================================================================*
# Settings

model = ffn.Case()

model.solvers.append(thSolver)
model.timeFolders = [timeFolder0]

model.add_function_object(massFlow)
model.add_function_object(TBulk)

settings: ffn.control.ControlDict = model.settings

settings.application = 'GeN-Foam'
settings.endTime = 20
settings.deltaT = 0.001
settings.writeControl = 'adjustableRunTime'
settings.writeInterval = 1
settings.writePrecision = 6
settings.timePrecision = 6
settings.runTimeModifiable = True
settings.adjustTimeStep = True
settings.maxDeltaT = 0.01

# settings.functionObjects.append(ffn.FunctionObject(
#     'vOutlet',
#     "surfaceFieldValue",
#     region=thMesh.region,
#     writeFields=False,
# ))

print(model)

# Export to OpenFOAM
model.export_to_openfoam()


#==============================================================================*
# Run

ffn.run(case=model, is_preprocessing=True)


#==============================================================================*
# Post-processing

massFlowData = massFlow.read_from_case(0)
TBulkData = TBulk.read_from_case(0)

fig, ax = plt.subplots()
ax.plot(TBulkData['Time'], TBulkData['TBulk'], label="Outlet Primary")
ax.set_xlabel('Time [s]')
ax.set_ylabel('Temperature [K]')
ax.legend()
fig.tight_layout()
fig.savefig("fig_results_TBulk.png")

model.plot_slice(
    region=thMesh,
    time=20,
    fieldName='T',
    cmap='RdBu_r',
    unit='K',
    offset=(0, 0.05, 0)
)
model.plot_animation(
    region=thMesh,
    fieldName='T',
    cmap='RdBu_r',
    unit='K',
    fps=2,
    normal='y'
)

#==============================================================================*
