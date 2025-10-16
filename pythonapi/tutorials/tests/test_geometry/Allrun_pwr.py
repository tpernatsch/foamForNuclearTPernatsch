"""
This model is a very simplified representation of a 3-loop PWR vessel.
"""

#==============================================================================*
# Imports

import foamForNuclear as ffn
import foamForNuclear.mesh as mesh
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.thermophysicalProperty as thermo

import numpy as np
import matplotlib.pyplot as plt


#==============================================================================*
# Geometric parameters

isParallel: bool = True

nLoops: int = 3

vesselRadius = 2.75
coreRadius = 2.2
gap = 0.1

coreLength = 4.6
legLength = 3

elbowRadius = 0.5
equivalentHydraulicDiameter = 0.6

nCenter: int = 10 # 3
nrDowncomer: int = 2 # 1


#==============================================================================*
# Mesh

primaryMesh = mesh.BlockMesh(region='primary')

primaryMesh.createHalfSphere(
    "center",
    radius=coreRadius,
    nCenter=nCenter, nBorder=nCenter,
    isAddBoundaryConditions=True
)
primaryMesh.createHollowHalfSphere(
    "ringSphere1",
    innerRadius=coreRadius,
    outerRadius=coreRadius+gap,
    nr=1, nt=nCenter,
    isAddBoundaryConditions=True
)
primaryMesh.createHollowHalfSphere(
    "ringSphere2",
    innerRadius=coreRadius+gap,
    outerRadius=vesselRadius,
    nr=nrDowncomer, nt=nCenter,
    isAddBoundaryConditions=True
)
primaryMesh.createCylinderAlongZ(
    name="core",
    radius=coreRadius,
    lowZ=0, highZ=coreLength,
    nx=nCenter, ny=nCenter, nz=3*nCenter,
    isAddBoundaryConditions=True
)
downcomer = primaryMesh.createRingAlongZ(
    name="downcomer",
    innerRadius=coreRadius+gap,
    outerRadius=vesselRadius,
    lowZ=0, highZ=coreLength,
    nr=nrDowncomer, nt=nCenter, nz=3*nCenter,
    isAddBoundaryConditions=True
)

inletManifold = primaryMesh.createPipeCylindricalManifoldAlongZ(
    'inletManifold',
    nEntries=nLoops,
    innerRadius=coreRadius+gap,
    outerRadius=vesselRadius,
    equivalentHydraulicDiameter=equivalentHydraulicDiameter,
    lowZ=coreLength,
    nr=nrDowncomer, nt=nCenter,
    isAddGap=True
)
outletManifold = primaryMesh.createPipeCylindricalManifoldAlongZ(
    'outletManifold',
    nEntries=nLoops,
    innerRadius=0.5,
    outerRadius=coreRadius,
    equivalentHydraulicDiameter=equivalentHydraulicDiameter,
    lowZ=coreLength,
    angleStart=360/(2*nLoops),
    nr=nCenter, nt=2*nCenter,
    isAddGap=False
)

# Add pipes
inletPipes, outletPipes = [], []
for i, inletPipe in enumerate(inletManifold[:nLoops]):
    inletPipes.append(
        primaryMesh.extrudeNormal(inletPipe, 'front', f'inletPipe{i}', length=legLength, n=4)
    )

for i, outletPipe in enumerate(outletManifold[::2]):
    outletPipes.append(
        primaryMesh.extrudeNormal(outletPipe, 'front', f'outletPipe{i}', length=legLength + vesselRadius-coreRadius, n=4)
    )

# Add inlet and outlet of manifold
inletManifoldOutlet = ffn.Face("inletManifoldOutlet_wall")
for block in inletManifold:
    inletManifoldOutlet.addSubFace(block.bottomFace())

outletManifoldInlet = ffn.Face("outletManifoldInlet_wall")
for block in outletManifold:
    outletManifoldInlet.addSubFace(block.bottomFace())

inlet = ffn.Face("inlet")
for block in inletPipes:
    inlet.addSubFace(block.frontFace())

outlet = ffn.Face("outlet")
for block in outletPipes:
    outlet.addSubFace(block.frontFace())


pipeWall = ffn.Face("pipeWall")
for block in outletPipes:
    pipeWall.addSubFace(block.leftFace())
    pipeWall.addSubFace(block.rightFace())
    pipeWall.addSubFace(block.topFace())
    pipeWall.addSubFace(block.bottomFace())
for block in inletPipes:
    pipeWall.addSubFace(block.leftFace())
    pipeWall.addSubFace(block.rightFace())
    pipeWall.addSubFace(block.topFace())
    pipeWall.addSubFace(block.bottomFace())
for block in inletManifold:
    pipeWall.addSubFace(block.topFace())
    pipeWall.addSubFace(block.backFace())
for i, block in enumerate(inletManifold[nLoops:]):
    pipeWall.addSubFace(block.frontFace())
    if (i % 2 == 0):
        pipeWall.addSubFace(block.rightFace())
    else:
        pipeWall.addSubFace(block.leftFace())
for block in outletManifold:
    pipeWall.addSubFace(block.topFace())
    pipeWall.addSubFace(block.backFace())
for block in outletManifold[1::2]:
    pipeWall.addSubFace(block.frontFace())



primaryMesh.addBoundary(inlet)
primaryMesh.addBoundary(outlet)
primaryMesh.addBoundary(inletManifoldOutlet)
primaryMesh.addBoundary(outletManifoldInlet)
primaryMesh.addBoundary(pipeWall)


primaryMesh.addMergePatchPairs()
# Merge
primaryMesh.mergePatchesWithName(name='wall', includeFacename=['Wall', 'ringSphere1Top'], patchType="wall")
# Rename
# primaryMesh.mergePatchesWithName(name='outlet', includeFacename=['coreTop_'], patchType="wall")
# primaryMesh.mergePatchesWithName(name='inlet', includeFacename=['downcomerTop_'], patchType="wall")
primaryMesh.mergePatchesWithName(name='coreTop_wall', includeFacename=['coreTop_'], patchType="wall")
primaryMesh.mergePatchesWithName(name='downcomerTop_wall', includeFacename=['downcomerTop_'], patchType="wall")


primaryMesh.mergePatchPairsByName("downcomerTop_wall", inletManifoldOutlet.name)
primaryMesh.mergePatchPairsByName("coreTop_wall", outletManifoldInlet.name)


#==============================================================================*
# Time folder

T_inlet = 290+273.15
massFlowRate = nLoops * 2500e3/3600

waterConstThermo = thermo.WaterConst(T=T_inlet)

# Sp = np.pi*(vesselRadius**2 - (coreRadius+gap)**2)
# # Pm = 2*np.pi*(vesselRadius + coreRadius+gap)
# # D = 4*Sp/Pm
# D = equivalentHydraulicDiameter


# rho_inlet = waterConstThermo.rho

# u_inlet = massFlowRate / (rho_inlet * Sp)

# Reynolds = u_inlet * rho_inlet * D / waterConstThermo.mu
# I = 0.16 * Reynolds**(-0.125) # 0.055*Reynolds**(-0.041)
# k_inlet = 3/2 * (u_inlet*I)**2

# lm = 0.07 * D
# epsilon_inlet = 0.09**(3/4) * k_inlet**(3/2) / lm

# print(massFlowRate, rho_inlet, Sp, D, u_inlet, epsilon_inlet, k_inlet, Reynolds, waterConstThermo.mu)


timeFolder0 = ffn.TimeFolder(0)

T = ffn.Field("T", region=primaryMesh.region)
T.dimensions = ffn.Dimension(default='T')
T.internalField = T_inlet
T.set_boundary_condition("inlet", bc.FixedValue(T.internalField))
T.set_boundary_condition("outlet", bc.ZeroGradient())
T.set_boundary_condition('".*wall"', bc.ZeroGradient())

U = ffn.Field("U", region=primaryMesh.region)
U.dimensions = ffn.Dimension(default='U')
U.internalField = ffn.Vector(0, 0, 0)
# U.set_boundary_condition("inlet", bc.ZeroGradient())
U.set_boundary_condition("inlet", bc.FlowRateInletVelocity(
    value=U.internalField,
    # volumetricFlowRate=0.001,#100000/3600/3,
    massFlowRate=massFlowRate,
))
U.set_boundary_condition("outlet", bc.ZeroGradient())
U.set_boundary_condition('".*wall"', bc.Slip())

p_rgh = ffn.Field("p_rgh", region=primaryMesh.region)
p_rgh.dimensions = ffn.Dimension(default='p')
p_rgh.internalField = 155e5
# p_rgh.set_boundary_condition("inlet", bc.FixedValue(p_rgh.internalField+1e5))
p_rgh.set_boundary_condition("inlet", bc.ZeroGradient())
p_rgh.set_boundary_condition("outlet", bc.FixedValue(p_rgh.internalField))
p_rgh.set_boundary_condition('".*wall"', bc.ZeroGradient())


# epsilon = ffn.Field("epsilon", region=primaryMesh.region)
# epsilon.dimensions = ffn.Dimension(length=2, time=-3)
# epsilon.internalField = epsilon_inlet
# epsilon.set_boundary_condition("inlet", bc.FixedValue(epsilon_inlet))
# epsilon.set_boundary_condition("outlet", bc.ZeroGradient())
# epsilon.set_boundary_condition('".*wall"', bc.EpsilonWallFunction(0))

# k = ffn.Field("k", region=primaryMesh.region)
# k.dimensions = ffn.Dimension(length=2, time=-2)
# k.internalField = k_inlet
# k.set_boundary_condition("inlet", bc.FixedValue(k_inlet))
# k.set_boundary_condition("outlet", bc.ZeroGradient())
# k.set_boundary_condition('".*wall"', bc.ZeroGradient())

# alphat = ffn.Field("alphat", region=primaryMesh.region)
# alphat.dimensions = ffn.Dimension(mass=1, length=-1, time=-1)
# alphat.internalField = 0
# alphat.set_boundary_condition("inlet", bc.FixedValue(0))
# alphat.set_boundary_condition("outlet", bc.FixedValue(0))
# alphat.set_boundary_condition('".*wall"', bc.ZeroGradient())

# nut = ffn.Field("nut", region=primaryMesh.region)
# nut.dimensions = ffn.Dimension(length=2, time=-1)
# nut.internalField = 0
# nut.set_boundary_condition("inlet", bc.FixedValue(0))
# nut.set_boundary_condition("outlet", bc.FixedValue(0))
# nut.set_boundary_condition('".*wall"', bc.ZeroGradient())

timeFolder0.append(T)
timeFolder0.append(U)
timeFolder0.append(p_rgh)
# timeFolder0.append(epsilon)
# timeFolder0.append(k)
# timeFolder0.append(alphat)
# timeFolder0.append(nut)


#==============================================================================*
# Solver

thSolver = ffn.ThermalHydraulicsSolver(
    region=primaryMesh.region,
    mesh=primaryMesh,
    solver='onePhase',
    removeBaffles=True,
    isSetFvSchemesToDefault=True,
    isSetFvSolutionToDefault=False
)
thSolver.thermophysicalProperties = thermo.WaterPolynomial()
# thSolver.thermophysicalProperties = thermo.WaterPerfectFluid()
# thSolver.thermophysicalProperties = waterConstThermo

thSolver.phaseProperties.pMin = 150e5

corePowerModel = ffn.StructureProperty(
    ['core'],
    elementDiameter=10.54e-3,
    pitch=12.85e-3,
    latticeType='square'
)

corePowerModel.powerModel = ffn.NuclearFuelPin(
    fuelInnerRadius=0,
    fuelOuterRadius=9.18e-3/2,
    cladInnerRadius=9.34e-3/2,
    cladOuterRadius=10.54e-3/2,
    fuelMeshSize=16,
    cladMeshSize=4,
    fuelRho=10480,
    fuelCp=250,
    fuelK=3,
    cladRho=7500,
    cladCp=500,
    cladK=20,
    fuelT=T.internalField,
    cladT=T.internalField,
    gapH=3000,
    powerDensity=100e6
)

corePowerModel.add_passive_structure(
    volumetricArea=2,
    rho=7700,
    Cp=500,
    T=T.internalField
)

thSolver.add_structure_property(corePowerModel)

thSolver.add_drag_model(
    ffn.ReynoldsPower(0.316, -0.25, zones=['core'])
)
thSolver.add_heat_transfer_model(
    ffn.NusseltReynoldsPrandtlPower(0, 0.023, 0.8, 0.4, zones=['core'])
)

# Turbulence properties
turbulenceProp = thSolver.turbulenceProperties
turbulenceProp.simulationType = 'laminar'
# turbulenceProp.simulationType = "RAS"
# turbulenceProp.RASoptions.RASModel = "porousKEpsilon"
# turbulenceProp.porousKEpsilonProperties.append(
#     ffn.PorousKEpsilonPropertiesPerZone(
#         zones=["core"],
#         DhStruct=corePowerModel.Dh
#     )
# )

# Pimple loop options
thSolver.pimpleOptions.nOuterCorrectors = 2
thSolver.pimpleOptions.nCorrectors = 2
thSolver.pimpleOptions.nNonOrthogonalCorrectors = 0
thSolver.pimpleOptions.momentumMode = "faceCentered"
thSolver.pimpleOptions.solveEnergy = True
thSolver.pimpleOptions.correctPhi = False


# Scheme definition
thSolver.fvSchemes.divSchemes['default'] = "Gauss upwind" # Better than Gauss linear
thSolver.fvSchemes.divSchemes['div(alphaRhoPhiNu,U)'] = "Gauss linear"

# Solver definition
thSolver.fvSolution.solvers['"p_rgh.*"'] = ffn.fvSolutionSolver(
    solver="GAMG",
    smoother="DIC",
    tolerance=1e-5,
    relTol=0
)
thSolver.fvSolution.solvers['"e|h"'] = ffn.fvSolutionSolver(
    solver="smoothSolver",
    smoother="symGaussSeidel",
    tolerance=1e-5,
    relTol=0,
    minIter=0
)
thSolver.fvSolution.solvers['".*"'] = ffn.fvSolutionSolver(
    solver="PBiCGStab",
    preconditioner="diagonal",
    tolerance=1e-5,
    relTol=0.01
)

# Run in parallel
if (isParallel):
    thSolver.decomposeParDict.numberOfSubdomains = 8
    thSolver.decomposeParDict.method = "scotch"


#==============================================================================*
# Function objects (inline post-processing)

massFlowInlet = ffn.MassFlow(
    name="massFlowInlet",
    log=True,
    writeFields=False,
    writeControl="adjustableRunTime",
    writeInterval=1,
    region=primaryMesh.region,
    regionName="inlet",
    scaleFactor=3600/1000/nLoops # Convert in t/h/loop
)

massFlowOutlet = ffn.MassFlow(
    name="massFlowOutlet",
    log=True,
    writeFields=False,
    writeControl="adjustableRunTime",
    writeInterval=massFlowInlet.writeInterval,
    region=primaryMesh.region,
    regionName="outlet",
    scaleFactor=3600/1000/nLoops # Convert in t/h/loop
)

TInlet = ffn.TBulk(
    name="TInlet",
    log=True,
    writeFields=False,
    writeControl="adjustableRunTime",
    writeInterval=1,
    region=primaryMesh.region,
    regionName="inlet"
)

TOutlet = ffn.TBulk(
    name="TOutlet",
    log=True,
    writeFields=False,
    writeControl="adjustableRunTime",
    writeInterval=massFlowInlet.writeInterval,
    region=primaryMesh.region,
    regionName="outlet"
)

fuelFraction = np.pi*corePowerModel.powerModel.fuelOuterRadius**2 / (corePowerModel.pitch**2)

totalPower = ffn.VolFieldValue(
    name="totalPower",
    fields=["powerDensityNeutronics"],
    operation="volIntegrate",
    region=thSolver.region,
    regionType="all",
    log=True,
    writeFields=False,
    writeControl="adjustableRunTime",
    writeInterval=massFlowInlet.writeInterval,
    scaleFactor=fuelFraction
)


#==============================================================================*
# Model

model = ffn.Model(timeFolders=[timeFolder0])
model.settings.application = "GeN-Foam"
model.settings.endTime = 200
model.settings.deltaT = 0.1
model.settings.writeInterval = 10
model.settings.writeControl = "adjustableRunTime"
model.settings.adjustTimeStep = True
model.settings.runTimeModifiable = True
model.settings.maxCo = 1
model.settings.maxDeltaT = 10

model.add_solver(thSolver)
model.add_function_object(massFlowInlet)
model.add_function_object(massFlowOutlet)
model.add_function_object(TInlet)
model.add_function_object(TOutlet)
model.add_function_object(totalPower)


#==============================================================================*

print(model)

# Export to OpenFOAM format
ffn.allclean()
model.export_to_openfoam()

# Preprocessing
ffn.run_preprocessing(model=model)

model.plot_mesh(region=primaryMesh, show_edges=True)
model.plot_mesh(region=primaryMesh, show_edges=True, normal="x")
model.plot_mesh(region=primaryMesh, show_edges=True, normal="y")
model.plot_mesh(region=primaryMesh, show_edges=True, normal="z")

model.plot_slice(region=primaryMesh, show_edges=True, normal="x")
model.plot_slice(region=primaryMesh, show_edges=True, normal="y")

for boundaryName in [f.name for f in primaryMesh.faces if not f.is_empty]+["defaultFaces"]:
    try:
        model.plot_boundary(region=primaryMesh, boundaryName=boundaryName, show_edges=True)
    except:
        pass

# Run the simulation
try:
    ffn.run(model=model)
except:
    pass

# Reconstruct the case at the end the simulation
if (model.is_parallel):
    ffn.reconstructParallelCase(model=model)


#==============================================================================*
# Post-processing

# Mass flow rate
dataMassFlowInlet = massFlowInlet.read_from_case(startTime=0)
dataMassFlowOutlet = massFlowOutlet.read_from_case(startTime=0)
dataTInlet = TInlet.read_from_case(startTime=0)
dataTOutlet = TOutlet.read_from_case(startTime=0)
dataPower = totalPower.read_from_case(startTime=0)

minT, maxT = list(dataTInlet["TBulk"])[-1], list(dataTOutlet["TBulk"])[-1]

print(f"Total power = {dataPower['volIntegrate(powerDensityNeutronics)'][0]/1e9:.2f} GW")

fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
ax.plot(dataMassFlowInlet["Time"], dataMassFlowInlet["MassFlow"], label="Inlet")
ax.plot(dataMassFlowOutlet["Time"], dataMassFlowOutlet["MassFlow"], label="Outlet")
ax.set_xlabel("Time [s]")
ax.set_ylabel("Mass flow rate [t/h/loop]")
ax.legend()
fig.tight_layout()
fig.savefig("fig_results_massFlow.png")
plt.close()

# Residuals
model.plot_residuals(
    parameters=['p_rgh', 'h', 'k', 'epsilon'],
    title="Thermal-hydraulics"
)

lastTime = model.get_time_steps()[-1]
print(f"Last time = {lastTime} s")

minP, maxP = model.get_data_range(primaryMesh, lastTime, "p_rgh")
minTfuel, maxTfuel = model.get_data_range(primaryMesh, lastTime, "T.fuelAvForNeutronics", removeZeros=True)
minTclad, maxTclad = model.get_data_range(primaryMesh, lastTime, "T.cladAvForNeutronics", removeZeros=True)

print(f"deltaT = {(maxT-minT):.2f} °C")
print(f"deltaP = {(maxP-minP)/1e5} bar")
print(f"Tfuel range = [{minTfuel:.1f}, {maxTfuel:.1f}] K")
print(f"Tclad range = [{minTclad:.1f}, {maxTclad:.1f}] K")

# Plot mesh and animation
for fieldname, unit in [
    ('T', 'K'),
    ('U', 'm/s'),
    ('p_rgh', 'Pa')
]:
    model.plot_mesh(
        region=primaryMesh,
        time=lastTime,
        fieldName=fieldname,
        unit=unit,
        cmap="RdBu_r",
        show_edges=True
    )
    model.plot_slice(
        region=primaryMesh,
        time=lastTime,
        fieldName=fieldname,
        unit=unit,
        cmap="RdBu_r",
        show_edges=True
    )
    model.plot_animation(
        region=primaryMesh,
        fieldName=fieldname,
        cmap="RdBu_r",
        unit=unit,
        fps=10,
        show_edges=True
    )
    model.plot_animation(
        region=primaryMesh,
        fieldName=fieldname,
        cmap="RdBu_r",
        unit=unit,
        fps=10,
        show_edges=True,
        normal='y'
    )

for fieldname, unit, minT, maxT in [
    ('T.cladAvForNeutronics', 'K', minTclad, maxTclad),
    ('T.fuelAvForNeutronics', 'K', minTfuel, maxTfuel),
]:
    model.plot_slice(
        region=primaryMesh,
        time=lastTime,
        fieldName=fieldname,
        unit=unit,
        cmap="RdBu_r",
        show_edges=True,
        limits=[minT, maxT]
    )


#==============================================================================*
