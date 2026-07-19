"""

"""
#==============================================================================*
# Imports

import matplotlib.pyplot as plt
import numpy as np

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh


#==============================================================================*
# nMesh

# Use https://foam-for-nuclear.gitlab.io/honeycomb/

lattice = """
0 0 0 0 0 F F F F F 0 0 0 0 0
0 0 0 F F F F F F F F F 0 0 0
0 0 F F F F F F F F F F F 0 0
0 F F F F F F F F F F F F F 0
0 F F F F F F F F F F F F F 0
F F F F F F F F F F F F F F F
F F F F F F F F F F F F F F F
F F F F F F F F F F F F F F F
F F F F F F F F F F F F F F F
F F F F F F F F F F F F F F F
0 F F F F F F F F F F F F F 0
0 F F F F F F F F F F F F F 0
0 0 F F F F F F F F F F F 0 0
0 0 0 F F F F F F F F F 0 0 0
0 0 0 0 0 F F F F F 0 0 0 0 0
"""


nXY = len([line for line in lattice.split("\n") if line != ""])

assemblyPitch = 0.215
fuelOD = 8.2e-3
cladOD = 9.5e-3
fuelPinPitch = 12.6e-3
coreHeight = 4.270

coreNodes = 15
nxyMesh = 1


def createMesh(region: str, isAddRing: bool=False):
    nMesh = mesh.BlockMesh(region=region)

    nMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: nMesh.create_cube(
            name="core",
            lowX=x-assemblyPitch/2,
            highX=x+assemblyPitch/2,
            lowY=y-assemblyPitch/2,
            highY=y+assemblyPitch/2,
            lowZ=-coreHeight/2,
            highZ=coreHeight/2,
            nx=nxyMesh,
            ny=nxyMesh,
            nz=coreNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="square",
        nx=nXY, ny=nXY,
        pitch=assemblyPitch,
        elementsToPlace=['F']
    )

    # Fill the gap
    if (isAddRing):
        nMesh.fill_lattice_ring_gap(
            name='gap',
            ringRadius=1.725,
            nxLat=nXY,
            nyLat=nXY,
            pitch=assemblyPitch,
            latticeType='square',
            nxBlock=nxyMesh,
            nyBlock=nxyMesh,
            nzBlock=coreNodes,
            zmin=-coreHeight/2,
            zmax=coreHeight/2,
            isAddAllBC=True
        )


    # Merge patches with same name pattern
    nMesh.add_merge_patch_pairs(
        includeFacename=['Top', 'Bottom'],
        excludeFacename=['Wall']
    )
    nMesh.add_merge_patch_pairs(
        includeFacename=['Wall'],
        excludeFacename=['Top', 'Bottom'],
    )
    externalWalls = nMesh.get_standalone_faces(
        includeFacename=['Wall'],
        excludeFacename=['Top', 'Bottom'],
    )

    # Rename patches
    nMesh.merge_patches_with_name(
        name="wall",
        includeFacename=[face.name for face in externalWalls],
        patchType="wall"
    )
    nMesh.merge_patches_with_name(name="bottom", includeFacename=["Bottom_"])
    nMesh.merge_patches_with_name(name="top", includeFacename=["Top_"])

    return(nMesh)


nMesh = createMesh(region="neutroRegion", isAddRing=False)
thMesh = createMesh(region="fluidRegion")


#==============================================================================*
# Fields

inletTemperature = 292.7 + 273.15

timeFolder0 = ffn.TimeFolder(time=0)

# Neutronics
defaultFlux = ffn.fields.Field("defaultFlux", region=nMesh.region)
defaultFlux.dimensions = ffn.fields.Dimension(default='flux')
defaultFlux.internalField = 1e21
defaultFlux.set_boundary_condition('wall', bc.FixedValue(0))
defaultFlux.set_boundary_condition('top', bc.FixedValue(0))
defaultFlux.set_boundary_condition('bottom', bc.FixedValue(0))

# Fluid
T = ffn.fields.Field("T", region=thMesh.region)
T.dimensions = ffn.fields.Dimension(default='T')
T.internalField = inletTemperature
T.set_boundary_condition('wall', bc.ZeroGradient())
T.set_boundary_condition('top', bc.ZeroGradient())
T.set_boundary_condition('bottom', bc.FixedValue(inletTemperature))

nominalMassFlowRate = 3*25000e3/3600 # (25000 t/h / loop)

U = ffn.fields.Field("U", region=thMesh.region)
U.dimensions = ffn.fields.Dimension(default='U')
U.internalField = ffn.common.Vector(0, 0, 1)
U.set_boundary_condition('wall', bc.Slip())
U.set_boundary_condition('top', bc.ZeroGradient())
U.set_boundary_condition('bottom', bc.FlowRateInletVelocity(
    value=U.internalField,
    massFlowRate=nominalMassFlowRate
))

p = ffn.fields.Field("p", region=thMesh.region)
p.dimensions = ffn.fields.Dimension(default='p')
p.internalField = 155e5
p.set_boundary_condition('wall', bc.ZeroGradient())
p.set_boundary_condition('top', bc.FixedValue(p.internalField))
p.set_boundary_condition('bottom', bc.ZeroGradient())

p_rgh = ffn.fields.Field("p_rgh", region=thMesh.region)
p_rgh.dimensions = ffn.fields.Dimension(default='p')
p_rgh.internalField = p.internalField
p_rgh.set_boundary_condition('wall', bc.ZeroGradient())
p_rgh.set_boundary_condition('top', bc.FixedValue(p_rgh.internalField))
p_rgh.set_boundary_condition('bottom', bc.FixedFluxExtrapolatedPressure(value=p_rgh.internalField, gradient=0))


timeFolder0.append(defaultFlux)
timeFolder0.append(T)
timeFolder0.append(U)
timeFolder0.append(p)
timeFolder0.append(p_rgh)


#==============================================================================*
# Neutronics solver

neutronicsSolver = ffn.solvers.NeutronicsSolver(
    solver="diffusionNeutronics",
    mesh=nMesh,
    region=nMesh.region,
    power=3e9,
    keff=0.9388902,
)

neutronicsSolver.nuclearData.import_from_openfoam("XS/nuclearData")


#==============================================================================*
# Thermal-hydraulics solver

thSolver = ffn.solvers.thermal_hydraulics.OnePhase(
    region="fluidRegion",
    removeBaffles=True,
    mesh=thMesh,
    isSetFvSolutionToDefault=False
)

thSolver.fluid.thermophysicalProperties = ffn.thermo.WaterConst(T=inletTemperature)

thSolver.fluid.turbulenceProperties.simulationType = 'laminar'

core = ffn.porous_medium.Structure(
    zones=["core"],
    pitch=fuelPinPitch,
    elementDiameter=cladOD,
    latticeType='square'
)
core.powerModel = ffn.porous_medium.power_models.NuclearFuelPin(
    fuelInnerRadius=0,
    fuelOuterRadius=fuelOD/2,
    cladInnerRadius=cladOD/2 - 0.57e-3,
    cladOuterRadius=cladOD/2,
    fuelMeshSize=30,
    cladMeshSize=5,
    fuelRho=10480,
    fuelCp=250,
    fuelK=3,
    cladRho=7500,
    cladCp=500,
    cladK=20,
    gapH=3000,
    fuelT=inletTemperature,
    cladT=inletTemperature
)
thSolver.structures.append(core)

thSolver.fluid_structure.dragModels.append(
    ffn.porous_medium.drag.ReynoldsPower(
        coeff=0.316, exp=-0.25,
        zones=["core"]
    )
)
thSolver.fluid_structure.heatTransferModels.append(
    ffn.porous_medium.heat_transfer.NusseltReynoldsPrandtlPower(
        const=0, coeff=0.023, expRe=0.8, expPr=0.4,
        zones=['core']
    )
)


thSolver.fvSchemes.divSchemes['default'] = 'none'
thSolver.fvSchemes.divSchemes['div(phi,alpha)'] = 'Gauss vanLeer'
thSolver.fvSchemes.divSchemes['div(phir,alpha)'] = 'Gauss vanLeer'
thSolver.fvSchemes.divSchemes['div\(phi.*,U.*\)'] = 'Gauss upwind'
thSolver.fvSchemes.divSchemes['div(alphaRhoPhi,U)'] = 'Gauss upwind'
thSolver.fvSchemes.divSchemes['div(alphaRhoPhiNu,U)'] = 'Gauss linear'
thSolver.fvSchemes.divSchemes['div(alphaRhoPhi,K)'] = 'Gauss upwind'
thSolver.fvSchemes.divSchemes['div\(alphaRhoPhi.*,k.*\)'] = 'Gauss upwind'
thSolver.fvSchemes.divSchemes['div\(alphaRhoPhi.*,epsilon.*\)'] = 'Gauss upwind'
thSolver.fvSchemes.divSchemes['div\(alphaRhoPhi.*,(h|e).*\)'] = 'Gauss upwind'


fvSolution = ffn.numerics.fvSolution()

fvSolution.append("p_rgh", ffn.numerics.fvSolutionSolver(
    solver='GAMG', smoother='DIC', tolerance=1e-7, relTol=0
))
fvSolution.append("p_rghFinal", ffn.numerics.fvSolutionSolver(
    solver='GAMG', smoother='DIC', tolerance=1e-7, relTol=0
))
fvSolution.append("e", ffn.numerics.fvSolutionSolver(
    solver='smoothSolver', smoother='symGaussSeidel', tolerance=1e-7, relTol=0, minIter=0
))
fvSolution.append("h", ffn.numerics.fvSolutionSolver(
    solver='smoothSolver', smoother='symGaussSeidel', tolerance=1e-7, relTol=0, minIter=0
))
fvSolution.append('".*"', ffn.numerics.fvSolutionSolver(
    solver='PBiCGStab', preconditioner='diagonal', tolerance=1e-7, relTol=0.001
))
thSolver.fvSolution = fvSolution

thSolver.pimpleOptions.nCorrectors = 2
thSolver.pimpleOptions.nOuterCorrectors = 6
thSolver.pimpleOptions.nNonOrthogonalCorrectors = 0
thSolver.pimpleOptions.momentumMode = "faceCentered"
thSolver.pimpleOptions.solveEnergy = True
thSolver.pimpleOptions.solveFluidMechanics = True


# thSolver.add_relaxation_on_equation('h', 0.7)
# thSolver.add_relaxation_on_field('T', 0.7)


#==============================================================================*
# Solvers

solvers = ffn.solvers.Solvers([neutronicsSolver, thSolver])

# for solver in solvers:
#     solver.decomposeParDict.numberOfSubdomains = 8
#     solver.decomposeParDict.method = "scotch"


#==============================================================================*
# Coupling

coupling = ffn.coupling.Coupling(solvers=solvers)

coupling.add_field_transfer(neutronicsSolver, thSolver, "powerDensity", "powerDensityStructure")
coupling.add_field_transfer(neutronicsSolver, thSolver, "secondaryPowerDensity", "powerDensityLiquid")

coupling.add_field_transfer(thSolver, neutronicsSolver, "T", "TCool")
coupling.add_field_transfer(thSolver, neutronicsSolver, "thermo:rho", "rhoCool")
coupling.add_field_transfer(thSolver, neutronicsSolver, "T.fuelAvForNeutronics", "TFuel")
coupling.add_field_transfer(thSolver, neutronicsSolver, "T.cladAvForNeutronics", "TClad")
coupling.add_field_transfer(thSolver, neutronicsSolver, "T.passiveStructure", "TStructMech")


#==============================================================================*
# Settings


model = ffn.Case(
    solvers=solvers,
    coupling=coupling,
    timeFolders=[timeFolder0]
)

settings = model.settings

settings.application = 'GeN-Foam'
settings.endTime = 400
settings.deltaT = 1
settings.writeControl = 'runTime'
settings.writeInterval = 100
settings.writePrecision = 7
settings.runTimeModifiable = True
settings.adjustTimeStep = False

# Add idxField to visualize cellZones, not used in calculations
idxField = neutronicsSolver.create_zone_field(nMesh.cellZones)
timeFolder0.append(idxField)

TOutletFuncObj = ffn.functions.SurfaceFieldValue(
    name=f'surfFieldValue_top',
    fields=['T'],
    operation="areaAverage",
    log=True,
    writeFields=False,
    writeControl=settings.writeControl,
    writeInterval=settings.writeInterval,
    region=thMesh.region,
    regionType='patch',
    regionName="top"
)

TInletFunObj = ffn.functions.SurfaceFieldValue(
    name=f'surfFieldValue_bottom',
    fields=['T'],
    operation="areaAverage",
    log=True,
    writeFields=False,
    writeControl=settings.writeControl,
    writeInterval=settings.writeInterval,
    region=thMesh.region,
    regionType='patch',
    regionName="bottom"
)

model.add_function_object(TOutletFuncObj)
model.add_function_object(TInletFunObj)

print(model)


if (True):
    # Export to OpenFOAM
    ffn.allclean()
    model.export_to_openfoam()

    coupling.plot_coupling_graph()
    coupling.plot_solving_graph()
    coupling.plot_solving_flowchart()


    #==============================================================================*
    # Duplicate and overwrite files

    # for filename in ['nuclearData']:
    #     ffn.copyFolder(f"./XS/{filename}", f"constant/{nMesh.region}")


    #==============================================================================*
    # Preprocessing

    ffn.run_preprocessing(model)

    model.plot_mesh(region=nMesh, show_edges=True)
    model.plot_mesh(region=thMesh, show_edges=True)
    model.plot_mesh(region=nMesh, fieldName=idxField.name, show_edges=True, cmap="tab10", limits=[0.5, 10.5])
    model.plot_slice(region=nMesh, fieldName=idxField.name, show_edges=True, cmap="tab10", normal="x", limits=[0.5, 10.5])
    model.plot_slice(region=nMesh, fieldName=idxField.name, show_edges=True, cmap="tab10", normal="y", limits=[0.5, 10.5])
    model.plot_slice(region=nMesh, fieldName=idxField.name, show_edges=True, cmap="tab10", normal="z", limits=[0.5, 10.5])
    for boundaryName in ['top', 'bottom', 'wall']:
        model.plot_boundary(region=nMesh, boundaryName=boundaryName)
        model.plot_boundary(region=thMesh, boundaryName=boundaryName)


    #==============================================================================*
    # Run

    ffn.run(model)


#==============================================================================*
# Post-processing

# Plot residuals
#---------------

model.plot_residuals(
    parameters=['fluxStar0'],
    title="Neutronics"
)
model.plot_residuals(
    parameters=['p_rgh', 'h'],
    title="Thermal-hydraulics"
)

lastTime = model.get_time_steps()[-1]


# Plot keff
#----------

print(f"keff = {model.keff(time=lastTime)}")

res = model.get_keff_from_log()

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(res['time'], res['keff'])
ax.set_xlabel("Time [s]")
ax.set_ylabel(r"k$_\text{eff}$")
fig.tight_layout()
fig.savefig("fig_results_keff.png")


# Temperature difference
#-----------------------

dataOutlet = TOutletFuncObj.read_from_case(startTime=int(settings.startTime))
dataInlet = TInletFunObj.read_from_case(startTime=int(settings.startTime))

Tin = np.array(dataInlet.tail(1))[0][1]
Tout = np.array(dataOutlet.tail(1))[0][1]

print(f"Core average deltaT = {Tout - Tin:g} K")


# Plot slices
#------------

model.plot_slice(
    region=nMesh,
    time=lastTime,
    fieldName='flux0',
    normal='z',
    cmap='Blues',
    unit='neutron/m2/s',
    show_edges=True,
    limits=[0, None]
)
model.plot_slice(
    region=nMesh,
    time=lastTime,
    fieldName='powerDensity',
    normal='z',
    cmap='inferno',
    unit='W/m3',
    show_edges=True
)

for field in ['TCool', 'TClad', 'TFuel']:
    minT, maxT = model.get_data_range(nMesh, lastTime, field, removeZeros=True)

    model.plot_slice(
        region=nMesh,
        time=lastTime,
        fieldName=field,
        normal='x',
        cmap='RdBu_r',
        unit='K',
        show_edges=True,
        limits=[minT, maxT]
    )
    model.plot_slice(
        region=nMesh,
        time=lastTime,
        fieldName=field,
        normal='z',
        cmap='RdBu_r',
        unit='K',
        show_edges=True,
        limits=[minT, maxT]
    )

model.plot_mesh(
    region=thMesh,
    time=lastTime,
    fieldName='T',
    cmap='RdBu_r',
    unit='K'
)


#==============================================================================*
