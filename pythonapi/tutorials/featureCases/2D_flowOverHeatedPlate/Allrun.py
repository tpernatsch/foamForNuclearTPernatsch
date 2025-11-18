"""

"""
#==============================================================================*
# Imports

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh

import matplotlib.pyplot as plt
import pandas as pd


#==============================================================================*
# Mesh

fluidMesh = mesh.BlockMesh(region="fluid")

fluidBlock2 = fluidMesh.create_cube("zone0", 0, 0, 0, 1, 0.5, 0.4, nx=200, ny=41, nz=1, gradx=5, grady=16)
fluidBlock3 = fluidMesh.extrude_right([fluidBlock2], "zone0", dx=2, nx=51)
fluidBlock1 = fluidMesh.extrude_left([fluidBlock2], "zone0", dx=0.5, nx=81, gradx=0.2)

fluidInlet = ffn.Face("inlet")
fluidInlet.add_sub_face(fluidBlock1.leftFace())

fluidOutlet = ffn.Face("outlet")
fluidOutlet.add_sub_face(fluidBlock3.rightFace())

fluidTop = ffn.Face("top", boundaryType="wall")
for face in [fluidBlock1, fluidBlock2, fluidBlock3]:
    fluidTop.add_sub_face(face.backFace())

fluidBottom = ffn.Face("bottom")
fluidBottom.add_sub_face(fluidBlock3.frontFace())

fluidSlip_bottom = ffn.Face("slip-bottom")
fluidSlip_bottom.add_sub_face(fluidBlock1.frontFace())

fluidInterface = ffn.Face(
    "interface",
    boundaryType="mappedWall",
    inGroups=["wall"],
    extraParameters={
        "sampleRegion": "solid",
        "samplePatch": "top",
        "sampleMode": "nearestPatchFace"
    }
)
fluidInterface.add_sub_face(fluidBlock2.frontFace())

fluidMesh.add_boundary(fluidInlet)
fluidMesh.add_boundary(fluidOutlet)
fluidMesh.add_boundary(fluidTop)
fluidMesh.add_boundary(fluidBottom)
fluidMesh.add_boundary(fluidSlip_bottom)
fluidMesh.add_boundary(fluidInterface)


solidMesh = mesh.BlockMesh(region="solid")

solidBlock = solidMesh.create_cube("solid", 0, -0.25, 0, 1, 0, 0.4, nx=200, ny=41, nz=1, gradx=5, grady=0.0625)

solidLeft = ffn.Face("left", boundaryType="wall")
solidLeft.add_sub_face(solidBlock.leftFace())

solidRight = ffn.Face("right", boundaryType="wall")
solidRight.add_sub_face(solidBlock.rightFace())

solidTop = ffn.Face(
    "top",
    boundaryType="mappedWall",
    inGroups=["wall"],
    extraParameters={
        "sampleRegion": "fluid",
        "samplePatch": "interface",
        "sampleMode": "nearestPatchFace"
    }
)
solidTop.add_sub_face(solidBlock.backFace())

solidBottom = ffn.Face("bottom", boundaryType="wall")
solidBottom.add_sub_face(solidBlock.frontFace())

solidMesh.add_boundary(solidLeft)
solidMesh.add_boundary(solidRight)
solidMesh.add_boundary(solidTop)
solidMesh.add_boundary(solidBottom)


#==============================================================================*
# Fields

timeFolder0 = ffn.TimeFolder(0)


Tfluid = ffn.Field("T", region=fluidMesh.region)
Tfluid.dimensions = ffn.Dimension(default='T')
Tfluid.internalField = 300
Tfluid.set_boundary_condition(fluidInterface, bc.Mixed(value=Tfluid.internalField, refValue=Tfluid.internalField, refGradient=0, valueFraction=0))
Tfluid.set_boundary_condition(fluidInlet, bc.FixedValue(Tfluid.internalField))
for face in [fluidOutlet, fluidTop, fluidBottom, fluidSlip_bottom]: # , outerWall]:
    Tfluid.set_boundary_condition(face, bc.ZeroGradient())
Tfluid.set_boundary_condition("defaultFaces", bc.Empty())

U = ffn.Field("U", region=fluidMesh.region)
U.dimensions = ffn.Dimension(default='U')
U.internalField = ffn.Vector(1, 0, 0)
U.set_boundary_condition(fluidInterface, bc.FixedValue(ffn.Vector(0, 0, 0)))
U.set_boundary_condition(fluidInlet, bc.FixedValue(U.internalField))
U.set_boundary_condition(fluidBottom, bc.FixedValue(U.internalField))
for face in [fluidOutlet, fluidTop, fluidSlip_bottom]:
    U.set_boundary_condition(face, bc.ZeroGradient())
U.set_boundary_condition("defaultFaces", bc.Empty())

p_rgh = ffn.Field("p_rgh", region=fluidMesh.region)
p_rgh.dimensions = ffn.Dimension(default="p")
p_rgh.internalField = 0
p_rgh.set_boundary_condition(fluidOutlet, bc.FixedValue(0))
for face in [fluidInterface, fluidInlet, fluidTop, fluidBottom, fluidSlip_bottom]:
    p_rgh.set_boundary_condition(face, bc.ZeroGradient())
p_rgh.set_boundary_condition("defaultFaces", bc.Empty())


Tsolid = ffn.Field("T", region=solidMesh.region)
Tsolid.dimensions = ffn.Dimension(default='T')
Tsolid.internalField = 310
Tsolid.set_boundary_condition(solidTop, bc.Mixed(value=Tsolid.internalField, refValue=Tsolid.internalField, refGradient=0, valueFraction=0))
Tsolid.set_boundary_condition(solidBottom, bc.FixedValue(Tsolid.internalField))
Tsolid.set_boundary_condition(solidLeft, bc.ZeroGradient())
Tsolid.set_boundary_condition(solidRight, bc.ZeroGradient())
Tsolid.set_boundary_condition("defaultFaces", bc.Empty())


timeFolder0.append(Tfluid)
timeFolder0.append(U)
timeFolder0.append(p_rgh)
timeFolder0.append(Tsolid)


#==============================================================================*
# Solvers

fluidSolver = ffn.ThermalHydraulicsSolver(
    region=fluidMesh.region,
    solver="onePhase",
    removeBaffles=False,
    mesh=fluidMesh,
    isSetFvSchemesToDefault=False,
    isSetFvSolutionToDefault=False
)
fluidSolver.g = ffn.Vector(0, 0, 0)

fluidSolver.turbulenceProperties.simulationType = "laminar"

thermo = ffn.thermophysicalProperty.SodiumBoussinesq()
thermo.description = ""
thermo.molWeight = 1
thermo.rho0 = 1
thermo.T0 = 303.0
thermo.beta = 0
thermo.Cp = 2000 # from Cp
thermo.Hf = 0
thermo.mu = 2.0e-4 # nu × rho
thermo.Pr = 0.01 # from Cp, lambda, mu

fluidSolver.thermophysicalProperties = thermo

fluidSolver.phaseProperties.pMin = -1e6
fluidSolver.phaseProperties.pRefCell = 0
fluidSolver.phaseProperties.pRefValue = 0

# fvSchemes
fluidSolver.fvSchemes.ddtSchemes['default'] = "Euler"
fluidSolver.fvSchemes.gradSchemes['default'] = "Gauss linear"
fluidSolver.fvSchemes.divSchemes['default'] = "none"
fluidSolver.fvSchemes.divSchemes['div(alphaRhoPhi,U)'] = "Gauss vanLeer"
fluidSolver.fvSchemes.divSchemes['div(alphaRhoPhi,he)'] = "Gauss vanLeer"
fluidSolver.fvSchemes.divSchemes['div((nuEff*dev2(T(grad(U)))))'] = "Gauss linear"
fluidSolver.fvSchemes.divSchemes['div((nuEff*dev(T(grad(U)))))'] = "Gauss linear"
fluidSolver.fvSchemes.divSchemes['div(alphaRhoPhiNu,U)'] = "Gauss linear"
fluidSolver.fvSchemes.divSchemes['div(alphaRhoPhi,K)'] = "Gauss upwind"
fluidSolver.fvSchemes.laplacianSchemes['default'] = "none"
fluidSolver.fvSchemes.laplacianSchemes['laplacian(((alpha*)*nuEff),U)'] = "Gauss linear corrected"
fluidSolver.fvSchemes.laplacianSchemes['laplacian(rAUf,p_rgh)'] = "Gauss linear corrected"
fluidSolver.fvSchemes.laplacianSchemes['laplacian((1|A(U)),p_rgh)'] = "Gauss linear corrected"
fluidSolver.fvSchemes.laplacianSchemes['laplacian(diffusivity,cellMotionU)'] = "Gauss linear corrected"
fluidSolver.fvSchemes.laplacianSchemes['laplacian((interpolate(alpha)*interpolate(alphaEff)),h)'] = "Gauss linear corrected"
fluidSolver.fvSchemes.snGradSchemes['default'] = "corrected"
fluidSolver.fvSchemes.fluxRequired['default'] = "no"
fluidSolver.fvSchemes.fluxRequired['p_rgh'] = ""
fluidSolver.fvSchemes.interpolationSchemes['default'] = "linear"

# fvSolution
fluidSolution = ffn.fvSolution()
fluidSolution.append('"p_rgh|p_rghFinal"', ffn.fvSolutionSolver(
    solver='GAMG',
    tolerance=1e-8,
    relTol=1e-3,
    minIter=1,
    maxIter=1000,
    smoother='GaussSeidel',
    nPreSweeps=0,
    nPostSweeps=2,
    nFinestSweeps=2,
    scaleCorrection=True,
    directSolveCoarsest=False,
    cacheAgglomeration=True,
    nCellsInCoarsestLevel=20,
    agglomerator="faceAreaPair",
    mergeLevels=1
))
fluidSolution.append('cellMotionU', ffn.fvSolutionSolver(
    solver='GAMG',
    tolerance=1e-6,
    relTol=1e-3,
    minIter=1,
    maxIter=1000,
    smoother='GaussSeidel',
    nPreSweeps=0,
    nPostSweeps=2,
    nFinestSweeps=2,
    scaleCorrection=True,
    directSolveCoarsest=False,
    cacheAgglomeration=True,
    nCellsInCoarsestLevel=20,
    agglomerator="faceAreaPair",
    mergeLevels=1
))
fluidSolution.append('"U|UFinal"', ffn.fvSolutionSolver(
    solver='PBiCG',
    preconditioner="DILU",
    tolerance=1e-8,
    relTol=1e-3,
    minIter=1
))
fluidSolution.append('"h|hFinal"', ffn.fvSolutionSolver(
    solver='PBiCG',
    preconditioner="DILU",
    tolerance=1e-12,
    relTol=1e-3,
    minIter=1
))

fluidSolver.fvSolution = fluidSolution


# Pimple Control
fluidSolver.pimpleOptions.nOuterCorrectors = 1
fluidSolver.pimpleOptions.nCorrectors = 3
fluidSolver.pimpleOptions.nNonOrthogonalCorrectors = 1
fluidSolver.pimpleOptions.solveEnergy = True
fluidSolver.pimpleOptions.solveFluidMechanics = True

fluidSolver.pimpleOptions.add_residual_control_on_field(
    fieldName="U",
    relTol=1e-6,
    tolerance=1e-6,
)

fluidSolver.add_relaxation_on_equation('"U|UFinal"', 0.9)


#==============================================================================*
# Thermal solver

solidSolver = ffn.OffbeatSolver(
    region=solidMesh.region,
    solver="extendedThermoMechanics",
    mesh=solidMesh,
    thermalSolverOptions=ffn.SolidConductionThermalSolverOptions(),
    mechanicsSolverOptions=ffn.MechanicsSolverOptions(),
    isSetFvSchemesToDefault=False,
    isSetFvSolutionToDefault=False
)
solidSolver.couplingOptions.correctTFromTH = False
solidSolver.couplingOptions.correctDispForNeutro = False

solidSolver.add_material(
    ffn.thermomechanicalMaterial.ConstantMaterial(
        name="solid",
        rho=1,
        Cp=100,
        k=100,
        emissivity=0,
        E=2e+11,
        nu=0.3,
        G=0,
        alpha=1e-5,
        Tref=0
    )
)

# fvSchemes
solidSolver.fvSchemes.d2dt2Schemes['default'] = "backward"
solidSolver.fvSchemes.ddtSchemes['default'] = "backward"
solidSolver.fvSchemes.gradSchemes['default'] = "Gauss linear"
solidSolver.fvSchemes.divSchemes['default'] = "Gauss linear"
solidSolver.fvSchemes.laplacianSchemes['default'] = "Gauss linear uncorrected"
solidSolver.fvSchemes.snGradSchemes['default'] = "uncorrected"
solidSolver.fvSchemes.interpolationSchemes['default'] = "linear"
solidSolver.fvSchemes.fluxRequired['default'] = "true"

# fvSolution
solidSolution = ffn.fvSolution()
solidSolution.append('T', ffn.fvSolutionSolver(
    solver='PCG',
    preconditioner='DIC',
    tolerance=1e-10,
    relTol=0.01,
    minIter=1,
    maxIter=100
))

solidSolver.fvSolution = solidSolution

solidSolver.stressAnalysis.nCorrectors = 5
solidSolver.stressAnalysis.maxOuterIter = 1
solidSolver.stressAnalysis.relT = 1e-6


#==============================================================================*
# Coupling

solvers = ffn.Solvers([fluidSolver, solidSolver])

chtLoop = ffn.CHTLoop(
    region="Level_1",
    solvers=[fluidSolver, solidSolver],
    minResidual=1e-6,
    maxIterations=100,
    fluidRegionName=fluidSolver.region,
    solidRegionName=solidSolver.region,
    fluidPatches=[fluidInterface.name],
    solidPatches=[solidTop.name],
    useHTC=False
)

coupling = ffn.Coupling(solvers=[chtLoop])


#==============================================================================*
# Settings

model = ffn.Model(
    timeFolders=[timeFolder0],
    solvers=solvers,
    coupling=coupling
)

settings: ffn.ControlDict = model.settings

settings.application = 'GeN-Foam'
settings.endTime = 8
settings.deltaT = 0.01
settings.writeControl = 'runTime'
settings.writeInterval = 0.5
settings.writePrecision = 6
settings.timePrecision = 6
settings.runTimeModifiable = True
settings.adjustTimeStep = False
settings.maxDeltaT = 1
settings.maxCo = 0.5

probesFunctionObject = ffn.Probes(
    name="probes",
    fields=["T"],
    enabled=True,
    writeControl="writeTime",
    writeInterval=1,
    region=fluidSolver.region,
    probeLocations=[
        (0.01, 0, 0),
        (0.06, 0, 0),
        (0.11, 0, 0),
        (0.16, 0, 0),
        (0.21, 0, 0),
        (0.26, 0, 0),
        (0.31, 0, 0),
        (0.36, 0, 0),
        (0.41, 0, 0),
        (0.46, 0, 0),
        (0.51, 0, 0),
        (0.56, 0, 0),
        (0.61, 0, 0),
        (0.66, 0, 0),
        (0.71, 0, 0),
        (0.76, 0, 0),
        (0.81, 0, 0),
        (0.86, 0, 0),
        (0.91, 0, 0),
        (0.96, 0, 0),
    ]
)

model.add_function_object(probesFunctionObject)

print(model)


#==============================================================================*
# Preprocessing

ffn.allclean()

# Export to OpenFOAM
model.export_to_openfoam()

coupling.plot_solving_flowchart()
coupling.plot_solving_graph()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=[fluidMesh, solidMesh], show_edges=True, normal="z")
model.plot_mesh(region=fluidMesh, show_edges=True, normal="z")
model.plot_mesh(region=solidMesh, show_edges=True, normal="z")


#==============================================================================*
# Run

ffn.run(model=model, is_preprocessing=False)


#==============================================================================*
# Post-processing

model.plot_residuals(
    parameters=["h", "p_rgh"],
    title="Thermal-hydraulics"
)

for fieldName, unit in [('T', 'K'), ("U", "m/s")]:
    model.plot_slice(
        region=fluidMesh,
        time=settings.endTime,
        fieldName=fieldName,
        cmap='RdBu_r',
        unit=unit,
        normal="z"
    )
    model.plot_animation(
        region=fluidMesh,
        fieldName=fieldName,
        cmap="RdBu_r",
        unit=unit,
        fps=2,
        normal="z"
    )

model.plot_slice(
    region=[solidMesh, fluidMesh],
    time=settings.endTime,
    fieldName='T',
    cmap='RdBu_r',
    unit='K',
    normal="z"
)
model.plot_animation(
    region=solidMesh,
    fieldName='T',
    cmap="RdBu_r",
    unit='K',
    fps=2,
    normal="z"
)


data, locations = probesFunctionObject.read_from_case(startTime=0, fieldName="T")

xRange = [loc.x for loc in locations]
lastTime = model.settings.endTime

# Normalize temperatures
T_norm = (data[lastTime] - 300) / 10

# Get expected values for comparison
numerical = pd.read_csv('benchmark/expected_numerical.csv', header=None, names=['x', 'y'])
analytical = pd.read_csv('benchmark/expected_analytical.csv', header=None, names=['x', 'y'])

# Plot
plt.figure(figsize=(8, 4))
plt.plot(xRange, T_norm, marker='o', label='foamForNuclear')
plt.plot(numerical['x'], numerical['y'], label='Expected Numerical', linestyle='--')
plt.plot(analytical['x'], analytical['y'], label='Expected Analytical', linestyle='-.')
plt.xlabel('Probe Index')
plt.ylabel(r'Normalized Temperature $\frac{T - 300}{T_s - 300}$')
plt.title(f'Normalized Temperature at t = {lastTime}s')
plt.xlim((0, 1))
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig("fig_results_temperatureComparison.png")


#==============================================================================*
