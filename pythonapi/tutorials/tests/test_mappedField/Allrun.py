import numpy as np
import foamForNuclear as ffn
import foamForNuclear.mesh as mesh
import foamForNuclear.boundaryConditions as bc


#==============================================================================*
# Mesh

# --- Region 1

th1Mesh = mesh.BlockMesh(region="fluid1Region")

core = th1Mesh.create_cube("core", 0, 0, 0, 1, 0.1, 1 - 0.2, nx=25, nz=25)

# Add clear fluid for proper coupling, specificaly for U field that may be
# affected by a porous structure
plenum1 = th1Mesh.extrude_top([core], 'plenum', dz=0.2, nz=6)

inlet1 = mesh.Face("inlet1")
inlet1.add_sub_face(core.bottomFace())

outlet1 = mesh.FaceMappedPatch(
    name="outlet1",
    sampleRegion='fluid2Region',
    sampleMode='nearestPatchFace',
    samplePatch='inlet2'
)
outlet1.add_sub_face(plenum1.topFace())

fixedWalls1 = mesh.Face("fixedWalls1", boundaryType="wall")
fixedWalls1.add_sub_face(core.leftFace())
fixedWalls1.add_sub_face(core.rightFace())
fixedWalls1.add_sub_face(plenum1.leftFace())
fixedWalls1.add_sub_face(plenum1.rightFace())

frontAndBack1 = mesh.Face("frontAndBack1", boundaryType="empty")
frontAndBack1.add_sub_face(core.frontFace())
frontAndBack1.add_sub_face(core.backFace())
frontAndBack1.add_sub_face(plenum1.frontFace())
frontAndBack1.add_sub_face(plenum1.backFace())

th1Mesh.add_boundary(inlet1)
th1Mesh.add_boundary(outlet1)
th1Mesh.add_boundary(fixedWalls1)
th1Mesh.add_boundary(frontAndBack1)


# --- Region 2

th2Mesh = mesh.BlockMesh(region="fluid2Region")

# Add clear fluid for proper coupling, specificaly for U field that may be
# affected by a porous structure
plenum2 = th2Mesh.create_cube("plenum", 0, 0, 1, 1, 0.1, 1.2, nx=25, nz=6)

hx = th2Mesh.extrude_top([plenum2], "hx", dz=1 - 0.2, nz=25)

inlet2 = mesh.FaceMappedPatch(
    name="inlet2",
    sampleRegion='fluid2Region',
    sampleMode='nearestPatchFace',
    samplePatch='outlet1'
)
inlet2.add_sub_face(plenum2.bottomFace())

outlet2 = mesh.Face("outlet2")
outlet2.add_sub_face(hx.topFace())

fixedWalls2 = mesh.Face("fixedWalls2", boundaryType="wall")
fixedWalls2.add_sub_face(hx.leftFace())
fixedWalls2.add_sub_face(hx.rightFace())
fixedWalls2.add_sub_face(plenum2.leftFace())
fixedWalls2.add_sub_face(plenum2.rightFace())

frontAndBack2 = mesh.Face("frontAndBack2", boundaryType="empty")
frontAndBack2.add_sub_face(hx.frontFace())
frontAndBack2.add_sub_face(hx.backFace())
frontAndBack2.add_sub_face(plenum2.frontFace())
frontAndBack2.add_sub_face(plenum2.backFace())

th2Mesh.add_boundary(inlet2)
th2Mesh.add_boundary(outlet2)
th2Mesh.add_boundary(fixedWalls2)
th2Mesh.add_boundary(frontAndBack2)


#==============================================================================*
# Time folder

# Coupling
# | Region 1 | Region 2 |
#
# - 1->2: zeroGrad U and   mapped p_rgh
# - 2->1:   mapped U and zeroGrad p_rgh

timeFolder0 = ffn.TimeFolder(0)

# --- Region 1

T1 = ffn.Field("T", region=th1Mesh.region)
T1.dimensions = ffn.Dimension(default='T')
T1.internalField = 600
T1.set_boundary_condition(fixedWalls1, bc.ZeroGradient())
T1.set_boundary_condition(inlet1, bc.FixedValue(T1.internalField))
T1.set_boundary_condition(outlet1, bc.ZeroGradient())
T1.set_boundary_condition(frontAndBack1, bc.Empty())

U1 = ffn.Field("U", region=th1Mesh.region)
U1.dimensions = ffn.Dimension(default='U')
U1.internalField = ffn.Vector(0, 0, 1)
U1.set_boundary_condition(fixedWalls1, bc.Slip())
U1.set_boundary_condition(inlet1, bc.UniformFixedValue(
    ffn.Table([
        (20, ffn.Vector(0, 0, 1)),
        (50, ffn.Vector(0, 0, 0.1)),
    ])
))
U1.set_boundary_condition(outlet1, bc.ZeroGradient())
U1.set_boundary_condition(frontAndBack1, bc.Empty())

p1 = ffn.Field("p", region=th1Mesh.region)
p1.dimensions = ffn.Dimension(default='p')
p1.internalField = 1e5
p1.set_boundary_condition(fixedWalls1, bc.Calculated(p1.internalField))
p1.set_boundary_condition(inlet1, bc.Calculated(p1.internalField))
p1.set_boundary_condition(outlet1, bc.Calculated(p1.internalField))
p1.set_boundary_condition(frontAndBack1, bc.Empty())

p_rgh1 = ffn.Field("p_rgh", region=th1Mesh.region)
p_rgh1.dimensions = ffn.Dimension(default='p')
p_rgh1.internalField = 1e5
p_rgh1.set_boundary_condition(fixedWalls1, bc.ZeroGradient())
p_rgh1.set_boundary_condition(inlet1, bc.FixedFluxPressure(p_rgh1.internalField))
p_rgh1.set_boundary_condition(outlet1, bc.MappedField(
    value=p_rgh1.internalField,
    sampleRegion=th2Mesh.region,
    sampleMode='nearestPatchFace',
    samplePatch='inlet2'
))
p_rgh1.set_boundary_condition(frontAndBack1, bc.Empty())


# --- Region 2

T2 = ffn.Field("T", region=th2Mesh.region)
T2.dimensions = ffn.Dimension(default='T')
T2.internalField = 600
T2.set_boundary_condition(fixedWalls2, bc.ZeroGradient())
T2.set_boundary_condition(inlet2, bc.MappedField(
    value=T2.internalField,
    sampleRegion=th1Mesh.region,
    sampleMode='nearestPatchFace',
    samplePatch='outlet1'
))
T2.set_boundary_condition(outlet2, bc.ZeroGradient())
T2.set_boundary_condition(frontAndBack2, bc.Empty())

U2 = ffn.Field("U", region=th2Mesh.region)
U2.dimensions = ffn.Dimension(default='U')
U2.internalField = ffn.Vector(0, 0, 1)
U2.set_boundary_condition(fixedWalls2, bc.Slip())
U2.set_boundary_condition(inlet2, bc.MappedField(
    value=U2.internalField,
    sampleRegion=th1Mesh.region,
    sampleMode='nearestPatchFace',
    samplePatch='outlet1',
))
U2.set_boundary_condition(outlet2, bc.ZeroGradient())
U2.set_boundary_condition(frontAndBack2, bc.Empty())

p2 = ffn.Field("p", region=th2Mesh.region)
p2.dimensions = ffn.Dimension(default='p')
p2.internalField = 1e5
p2.set_boundary_condition(fixedWalls2, bc.Calculated(p2.internalField))
p2.set_boundary_condition(inlet2, bc.Calculated(p2.internalField))
p2.set_boundary_condition(outlet2, bc.Calculated(p2.internalField))
p2.set_boundary_condition(frontAndBack2, bc.Empty())

p_rgh2 = ffn.Field("p_rgh", region=th2Mesh.region)
p_rgh2.dimensions = ffn.Dimension(default='p')
p_rgh2.internalField = 1e5
p_rgh2.set_boundary_condition(fixedWalls2, bc.ZeroGradient())
p_rgh2.set_boundary_condition(inlet2, bc.ZeroGradient())
p_rgh2.set_boundary_condition(outlet2, bc.FixedValue(p_rgh2.internalField))
p_rgh2.set_boundary_condition(frontAndBack2, bc.Empty())

timeFolder0.append(T1)
timeFolder0.append(U1)
timeFolder0.append(p1)
timeFolder0.append(p_rgh1)
timeFolder0.append(T2)
timeFolder0.append(U2)
timeFolder0.append(p2)
timeFolder0.append(p_rgh2)


#==============================================================================*
# Thermal-hydraulics solver

# --- Region 1

thSolver1 = ffn.ThermalHydraulicsSolver(
    th1Mesh.region, "onePhase",
    mesh=th1Mesh,
    removeBaffles=True,
    isSetFvSolutionToDefault=False
)

thSolver1.thermophysicalProperties = ffn.thermophysicalProperty.SodiumConst()

thSolver1.turbulenceProperties.simulationType = 'laminar'

core = ffn.StructureProperty(['core'], volumeFraction=0.5, Dh=0.01)

core.powerModel = ffn.NuclearFuelPin(
    powerDensity=100e6,
    fuelInnerRadius=0.0012,
    fuelOuterRadius=0.004715,
    cladInnerRadius=0.004865,
    cladOuterRadius=0.005365,
    fuelMeshSize=16,
    cladMeshSize=4,
    fuelRho=10480,
    fuelCp=250,
    fuelK=3,
    cladRho=7500,
    cladCp=500,
    cladK=20,
    fuelT=600,
    cladT=600,
    gapH=3000
)

core.add_passive_structure(
    volumetricArea=2,
    rho=7700,
    Cp=500,
    T=600
)

thSolver1.add_structure_property(core)

thSolver1.add_drag_model(
    ffn.ReynoldsPower(0.687, -0.25, zones=['core'])
)
thSolver1.add_heat_transfer_model(
    ffn.NusseltReynoldsPrandtlPower(4.82, 0.0185, 0.827, 0.827, zones=['core'])
)


thSolution = ffn.fvSolution()

thSolution.append('"p_rgh.*"', ffn.fvSolutionSolver(
    solver='GAMG', smoother='DIC', tolerance=1e-8, relTol=0
))
thSolution.append('"e.*"', ffn.fvSolutionSolver(
    solver='smoothSolver',
    smoother='symGaussSeidel',
    tolerance=1e-8, relTol=0, minIter=0
))
thSolution.append('"h.*"', ffn.fvSolutionSolver(
    solver='smoothSolver',
    smoother='symGaussSeidel',
    tolerance=1e-8, relTol=0, minIter=0
))
thSolution.append('".*"', ffn.fvSolutionSolver(
    solver='PBiCGStab',
    preconditioner='diagonal',
    tolerance=1e-6, relTol=0.001
))

thSolver1.fvSolution = thSolution

thSolver1.pimpleOptions.nOuterCorrectors = 3


thSolver1.fvSchemes.divSchemes['default'] = 'none'
thSolver1.fvSchemes.divSchemes['div(phi,alpha)'] = 'Gauss vanLeer'
thSolver1.fvSchemes.divSchemes['div(phir,alpha)'] = 'Gauss vanLeer'
thSolver1.fvSchemes.divSchemes['div\(phi.*,U.*\)'] = 'Gauss upwind'
thSolver1.fvSchemes.divSchemes['div(alphaRhoPhi,U)'] = 'Gauss upwind'
thSolver1.fvSchemes.divSchemes['div(alphaRhoPhiNu,U)'] = 'Gauss linear'
thSolver1.fvSchemes.divSchemes['div(alphaRhoPhi,K)'] = 'Gauss upwind'
thSolver1.fvSchemes.divSchemes['div\(alphaRhoPhi.*,k.*\)'] = 'Gauss upwind'
thSolver1.fvSchemes.divSchemes['div\(alphaRhoPhi.*,epsilon.*\)'] = 'Gauss upwind'
thSolver1.fvSchemes.divSchemes['div\(alphaRhoPhi.*,(h|e).*\)'] = 'Gauss upwind'

thSolver1.fvSchemes.laplacianSchemes['default'] = 'Gauss linear uncorrected'

thSolver1.fvSchemes.snGradSchemes['default'] = 'uncorrected'


# --- Region 2

thSolver2 = ffn.ThermalHydraulicsSolver(
    th2Mesh.region, "onePhase",
    mesh=th2Mesh,
    removeBaffles=True,
    isSetFvSolutionToDefault=False
)

thSolver2.thermophysicalProperties = ffn.thermophysicalProperty.SodiumConst()

thSolver2.turbulenceProperties.simulationType = 'laminar'

hx = ffn.StructureProperty(['hx'], volumeFraction=0.5, Dh=0.01)

thSolver2.add_structure_property(hx)

thSolver2.add_drag_model(
    ffn.ReynoldsPower(0.687, -0.25, zones=['hx'])
)
thSolver2.add_heat_transfer_model(
    ffn.NusseltReynoldsPrandtlPower(4.82, 0.0185, 0.827, 0.827, zones=['hx'])
)


thSolver2.fvSolution = thSolution

thSolver2.pimpleOptions.nOuterCorrectors = 3


thSolver2.fvSchemes.divSchemes['default'] = 'none'
thSolver2.fvSchemes.divSchemes['div(phi,alpha)'] = 'Gauss vanLeer'
thSolver2.fvSchemes.divSchemes['div(phir,alpha)'] = 'Gauss vanLeer'
thSolver2.fvSchemes.divSchemes['div\(phi.*,U.*\)'] = 'Gauss upwind'
thSolver2.fvSchemes.divSchemes['div(alphaRhoPhi,U)'] = 'Gauss upwind'
thSolver2.fvSchemes.divSchemes['div(alphaRhoPhiNu,U)'] = 'Gauss linear'
thSolver2.fvSchemes.divSchemes['div(alphaRhoPhi,K)'] = 'Gauss upwind'
thSolver2.fvSchemes.divSchemes['div\(alphaRhoPhi.*,k.*\)'] = 'Gauss upwind'
thSolver2.fvSchemes.divSchemes['div\(alphaRhoPhi.*,epsilon.*\)'] = 'Gauss upwind'
thSolver2.fvSchemes.divSchemes['div\(alphaRhoPhi.*,(h|e).*\)'] = 'Gauss upwind'

thSolver2.fvSchemes.laplacianSchemes['default'] = 'Gauss linear uncorrected'

thSolver2.fvSchemes.snGradSchemes['default'] = 'uncorrected'


#==============================================================================*
# Solvers

solvers = ffn.Solvers([thSolver1, thSolver2])


#==============================================================================*
# Coupling

fluidLoop = ffn.MultiPhysicsLoop(
    solver="picardLoop",
    region="fluidLoop",
    maxResidual=1e-6,
    maxIterations=1,
    solvers=[thSolver1, thSolver2]
)

coupling = ffn.Coupling(solvers=[fluidLoop])

coupling.plot_solving_flowchart()
coupling.plot_solving_graph()


#==============================================================================*
# Model

model = ffn.Model(
    solvers=solvers,
    coupling=coupling,
    timeFolders=[timeFolder0]
)

settings = model.settings

settings.application = "GeN-Foam"
settings.startTime = 0
settings.endTime = 100
settings.deltaT = 1
settings.writeControl = 'adjustableRunTime'
settings.writeInterval = 2
settings.adjustTimeStep = True

print(model)

ffn.allclean()
model.export_to_openfoam()


#==============================================================================*
# Run

ffn.run(model, is_preprocessing=True)


#==============================================================================*
# Post-processing

for mesh_ in [th1Mesh, th2Mesh]:
    model.plot_slice(
        region=mesh_.region,
        time=settings.endTime,
        fieldName='U',
        cmap='RdBu_r',
        show_edges=True,
        unit='m/s'
    )
    model.plot_slice(
        region=mesh_.region,
        time=settings.endTime,
        fieldName='T',
        cmap='RdBu_r',
        show_edges=True,
        unit='K'
    )
    model.plot_slice(
        region=mesh_.region,
        time=settings.endTime,
        fieldName='p_rgh',
        cmap='RdBu_r',
        show_edges=True,
        unit='Pa'
    )


#==============================================================================*
