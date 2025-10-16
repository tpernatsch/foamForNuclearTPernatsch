#==============================================================================*
# Imports

import foamForNuclear as ffn
import foamForNuclear.mesh as mesh
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.thermomechanicalMaterial as mat

import matplotlib.pyplot as plt


#==============================================================================*

fuelLength = 3000
wedgeAngle = 0.25

fuelRadius = 4.5

fuelMesh = mesh.BlockMesh(scale=0.001)

fuelBlock = fuelMesh.createWedge(
    'fuel', 0.0, fuelRadius, 0, fuelLength, wedgeAngle=wedgeAngle, nr=30, nz=10
)

cladBlock = fuelMesh.createWedge(
    'cladding', 4.565, 5.315, 0, fuelLength, wedgeAngle=wedgeAngle, nr=10, nz=10
)
cladBlockTop = fuelMesh.extrudeTop([cladBlock], 'cladding', dz=200, nz=1)


bcFuelFront = mesh.Face(name="fuelFront", boundaryType="wedge")
bcFuelFront.addSubFace(fuelBlock.frontFace())

bcFuelBack = mesh.Face(name="fuelBack", boundaryType="wedge")
bcFuelBack.addSubFace(fuelBlock.backFace())

bcCladFront = mesh.Face(name="cladFront", boundaryType="wedge")
bcCladFront.addSubFace(cladBlock.frontFace())
bcCladFront.addSubFace(cladBlockTop.frontFace())

bcCladBack = mesh.Face(name="cladBack", boundaryType="wedge")
bcCladBack.addSubFace(cladBlock.backFace())
bcCladBack.addSubFace(cladBlockTop.backFace())


bcFuelOuter = mesh.Face(
    name="fuelOuter",
    boundaryType="regionCoupledOFFBEAT",
    extraParameters={
        "neighbourPatch": "cladInner",
        "neighbourRegion": "region0",
        "owner": 'true',
        'updateAMI': 'true'
    }
)
bcFuelOuter.addSubFace(fuelBlock.rightFace())

bcCladInner = mesh.Face(
    name="cladInner",
    boundaryType="regionCoupledOFFBEAT",
    extraParameters={
        "neighbourPatch": "fuelOuter",
        "neighbourRegion": "region0",
        "owner": 'false',
        'updateAMI': 'true'
    }
)
bcCladInner.addSubFace(cladBlock.leftFace())
bcCladInner.addSubFace(cladBlockTop.leftFace())

bcCladOuter = mesh.Face(name="cladOuter")
bcCladOuter.addSubFace(cladBlock.rightFace())
bcCladOuter.addSubFace(cladBlockTop.rightFace())

for bcName, subface in [
    ("fuelBottom", fuelBlock.bottomFace()),
    ("fuelTop", fuelBlock.topFace()),
    ("cladBottom", cladBlock.bottomFace()),
    ("cladTop", cladBlockTop.topFace()),
]:
    bcEmpty = mesh.Face(name=bcName, boundaryType="empty")
    bcEmpty.addSubFace(subface)
    fuelMesh.addBoundary(bcEmpty)

fuelMesh.addBoundary(bcFuelOuter)
fuelMesh.addBoundary(bcCladInner)
fuelMesh.addBoundary(bcCladOuter)
fuelMesh.addBoundary(bcFuelFront)
fuelMesh.addBoundary(bcFuelBack)
fuelMesh.addBoundary(bcCladFront)
fuelMesh.addBoundary(bcCladBack)


#==============================================================================*
# Time solver

timeFolder0 = ffn.TimeFolder(0)

T = ffn.Field("T")
T.dimensions = ffn.Dimension(default='T')
T.internalField = 300
T.set_boundary_condition('.*Front|.*Back', bc.Wedge())
T.set_boundary_condition(bcFuelOuter, bc.FuelRodGap(value=T.internalField, roughness=2.2e-6))
T.set_boundary_condition(bcCladInner, bc.FuelRodGap(value=T.internalField, roughness=0.5e-6))
T.set_boundary_condition(bcCladOuter, bc.FixedValue(value=600))
T.set_boundary_condition(".*Top|.*Bottom|defaultFaces", bc.Empty())

D = ffn.Field("D")
D.dimensions = ffn.Dimension(length=1)
D.internalField = ffn.Vector(0, 0, 0)
D.set_boundary_condition('.*Front|.*Back', bc.Wedge())
D.set_boundary_condition(f'{bcFuelOuter.name}|{bcCladInner.name}', bc.GapContact(
    value=D.internalField,
    penaltyFactor=0.1,
    frictionCoefficient=0,
    relax=1.0,
    relaxInterfacePressure=0.1
))
coolantPressureList = ffn.Table([
    (0, 1e5),
    (3600, 1e7),
    (1.26e8, 1e7),
])
D.set_boundary_condition(bcCladOuter, bc.CoolantPressure(
    value=D.internalField,
    relax=1.0,
    coolantPressureList=coolantPressureList,
    outOfBounds='clamp'
))
D.set_boundary_condition(".*Top|.*Bottom|defaultFaces", bc.Empty())

neutronFlux0 = ffn.Field("neutronFlux0")
neutronFlux0.dimensions = ffn.Dimension(default='flux')
neutronFlux0.internalField = 0
neutronFlux0.set_boundary_condition('.*Front|.*Back', bc.Wedge())
neutronFlux0.set_boundary_condition('cladInner|cladOuter', bc.FixedValue(0))
neutronFlux0.set_boundary_condition(bcFuelOuter, bc.FixedValue(1))
neutronFlux0.set_boundary_condition(".*Top|.*Bottom|defaultFaces", bc.Empty())

gapGas = ffn.GapGas(gasPressureType="fromModel", gasPressure=2.25e6)

timeFolder0.append(T)
timeFolder0.append(D)
timeFolder0.append(neutronFlux0)
timeFolder0.append(gapGas)


#==============================================================================*
# Solver

offbeatSolver = ffn.OffbeatSolver(
    solver='offbeat',
    mesh=fuelMesh,
    thermalSolverOptions=ffn.SolidConductionThermalSolverOptions(),
    mechanicsSolverOptions=ffn.SmallStrainMechanicsSolverOptions(cylindricalStress=True),
    neutronicsSolverOptions=ffn.DiffusionNeutronicsSolverOptions(),
    materialProperties="byZone",
    burnupOptions=ffn.LassmannBurnupOptions(),
    fgrOptions=ffn.SciantixFgrOptions(nFrequency=1, relax=1)
)

offbeatSolver.globalOptions.pinDirection = ffn.Vector(0, 0, 1)
offbeatSolver.globalOptions.reactorType = "LWR"

# offbeatSolver.thermalSolverOptions.heatFluxSummary = False
# offbeatSolver.thermalSolverOptions.patchForAverage = "fuelOuter"

offbeatSolver.rheologyOptions.thermalExpansion = True
offbeatSolver.rheologyOptions.modifiedPlaneStrain = True
offbeatSolver.rheologyOptions.springModulus = 3500
offbeatSolver.rheologyOptions.planeStress = False
offbeatSolver.rheologyOptions.coolantPressureList = coolantPressureList
offbeatSolver.rheologyOptions.outOfBounds = "clamp"

# offbeatSolver.mechanicsSolverOptions.forceSummary = False
# offbeatSolver.mechanicsSolverOptions.cylindricalStress = True


offbeatSolver.gapGasOptions = ffn.FrapconGapGasOptions(
    gapVolumeOffset=0,
    gasReserveVolume=0,
    gasReserveTemperature=290,
    gapPatches=[bcFuelOuter.name, bcCladInner.name],
    topFuelPatches=["fuelTop"],
    bottomFuelPatches=["fuelBottom"]
)

offbeatSolver.heatSourceOptions = ffn.TimeDependentLhgrHeatSourceOptions(
    timePoints=[0, 3600, 1.26E+08],
    lhgr=[0, 200e2, 100e2],
    materials=[fuelBlock.name],
    axialProfile=ffn.TimeDependentTabulatedAxialProfile(
        timePoints=[
            0.00E+00, 1.58E+07, 3.15E+07, 3.16E+07, 4.73E+07, 6.31E+07, 6.32E+07,
            7.88E+07, 9.46E+07, 9.47E+07, 1.10E+08, 1.26E+08
        ],
        axialLocations=[
            0, 0.090909090909091, 0.181818181818182, 0.272727272727273,
            0.363636363636364, 0.454545454545455, 0.545454545454546,
            0.636363636363636, 0.727272727272727, 0.818181818181818,
            0.909090909090909, 1
        ],
        data=[
            ( 0.469787352107652, 0.913850849204304, 1.19921124080937, 1.32281654412493, 1.3457064151093, 1.29840068174161, 1.21141917200103, 1.08628787728651, 0.942844685784499, 0.767355674904378, 0.541508947858657, 0.271408470243168 ),
            ( 0.625964954658454, 1.03050615457391, 1.13900379439798, 1.13280392926518, 1.10335456988436, 1.08630494076915, 1.07235524422034, 1.06305544652113, 1.04135591855632, 0.977807300945077, 0.811960908642572, 0.457018629789547 ),
            ( 0.807036867064403, 1.11329355973688, 1.10701137116924, 1.05518331548621, 1.01749018408036, 1.00178471266126, 0.997073071235531, 1.00335525980317, 1.02063127836418, 1.023772372648, 0.942103921268676, 0.629565040028555 ),
            ( 0.469787352107652, 0.913850849204304, 1.19921124080937, 1.32281654412493, 1.3457064151093, 1.29840068174161, 1.21141917200103, 1.08628787728651, 0.942844685784499, 0.767355674904378, 0.541508947858657, 0.271408470243168 ),
            ( 0.625964954658454, 1.03050615457391, 1.13900379439798, 1.13280392926518, 1.10335456988436, 1.08630494076915, 1.07235524422034, 1.06305544652113, 1.04135591855632, 0.977807300945077, 0.811960908642572, 0.457018629789547 ),
            ( 0.807036867064403, 1.11329355973688, 1.10701137116924, 1.05518331548621, 1.01749018408036, 1.00178471266126, 0.997073071235531, 1.00335525980317, 1.02063127836418, 1.023772372648, 0.942103921268676, 0.629565040028555 ),
            ( 0.469787352107652, 0.913850849204304, 1.19921124080937, 1.32281654412493, 1.3457064151093, 1.29840068174161, 1.21141917200103, 1.08628787728651, 0.942844685784499, 0.767355674904378, 0.541508947858657, 0.271408470243168 ),
            ( 0.625964954658454, 1.03050615457391, 1.13900379439798, 1.13280392926518, 1.10335456988436, 1.08630494076915, 1.07235524422034, 1.06305544652113, 1.04135591855632, 0.977807300945077, 0.811960908642572, 0.457018629789547 ),
            ( 0.807036867064403, 1.11329355973688, 1.10701137116924, 1.05518331548621, 1.01749018408036, 1.00178471266126, 0.997073071235531, 1.00335525980317, 1.02063127836418, 1.023772372648, 0.942103921268676, 0.629565040028555 ),
            ( 0.469787352107652, 0.913850849204304, 1.19921124080937, 1.32281654412493, 1.3457064151093, 1.29840068174161, 1.21141917200103, 1.08628787728651, 0.942844685784499, 0.767355674904378, 0.541508947858657, 0.271408470243168 ),
            ( 0.625964954658454, 1.03050615457391, 1.13900379439798, 1.13280392926518, 1.10335456988436, 1.08630494076915, 1.07235524422034, 1.06305544652113, 1.04135591855632, 0.977807300945077, 0.811960908642572, 0.457018629789547 ),
            ( 0.807036867064403, 1.11329355973688, 1.10701137116924, 1.05518331548621, 1.01749018408036, 1.00178471266126, 0.997073071235531, 1.00335525980317, 1.02063127836418, 1.023772372648, 0.942103921268676, 0.629565040028555 ),
        ],
        axialInterpolationMethod="linear",
        burnupInterpolationMethod="linear",
    ),
    radialProfile=ffn.FromBurnupRadialProfile()
)

offbeatSolver.fastFluxOptions = ffn.TimeDependentFastFluxOptions(
    timePoints=[0, 3600, 1.26e8],
    fastFlux=[0, 2e13, 1e13],
    timeInterpolationMethod="linear",
    materials=[fuelBlock.name, cladBlock.name],
    axialProfile=ffn.FlatAxialProfile()
)

offbeatSolver.add_material(mat.UO2(
    name=fuelBlock.name,
    Tref=293,
    enrichment=0.045,
    rGrain=2.8e-5,
    GdContent=0,
    theoreticalDensity=10960,
    densityFraction=0.95,
    dishFraction=0,
    resinteringDensityChange=0.3,
    GapCold=0.13e-3,
    DiamCold=fuelRadius * fuelMesh.scale,
    recoveryFraction=0.5,
    outerPatch="fuelOuter",
    isotropicCracking=True,
    nCracksMax=12,
    rheologyModel=mat.ElasticityRheologyModel()
))

offbeatSolver.add_material(mat.Zircalloy(
    name=cladBlock.name,
    Tref=293,
    rheologyModel=mat.MisesPlasticCreepRheologyModel(
        yieldStressModel=mat.HardningYieldStressModel(ffn.Table([
            (0, 250e6)
        ])),
        creepModel=mat.LimbackCreepModel()
    )
))

offbeatSolver.fvSchemes.d2dt2Schemes['default'] = 'steadyState'

offbeatSolver.stressAnalysis.nCorrectors = 3
offbeatSolver.stressAnalysis.maxOuterIter = 1000
offbeatSolver.stressAnalysis.D = ffn.Vector(1e-5, 1, 1e-5)
offbeatSolver.stressAnalysis.T = 1e-5
offbeatSolver.stressAnalysis.neutronFlux = 1e-5
offbeatSolver.stressAnalysis.relD = 1e-6
offbeatSolver.stressAnalysis.relT = 1e-6
offbeatSolver.stressAnalysis.relneutronFlux0 = 1e-6

offbeatSolver.fvSolution.relaxationFactors['fields']['D'] = 0.8


#==============================================================================*
# Settings

model = ffn.Model()
model.add_solver(offbeatSolver)
model.add_time_folder(timeFolder0)

settings: ffn.ControlDict = model.settings

settings.application = 'offbeat'
settings.endTime = 3.15E+07
settings.deltaT = 1
settings.writeControl = 'timeStep'
settings.writeInterval = 100
settings.writePrecision = 6
settings.timePrecision = 20
settings.runTimeModifiable = True
settings.adjustTimeStep = True
settings.maxDeltaT = 6.048e5
settings.minDeltaT = 0.001
settings.maxRelativeDeltaTIncrease = 1e9
settings.minRelativeDeltaTDecrease = 1e9
settings.maxRelativePowerIncrease = 1e9
settings.maxRelativePowerDecrease = 1e9
settings.maxBurnupIncrease = 0.1
settings.maxAverageCreep = 1e-4
settings.maxMaximumCreep = 1e-4
settings.maxFGR = 1e-8

probes = ffn.Probes(
    name="probes",
    fields=['T', 'Bu'],
    probeLocations=[ffn.Vector(0, 0, 1.5)],
    enabled=True
)

settings.add_function_object(probes)


print(model)


#==============================================================================*
# Export to OpenFOAM

ffn.allclean()
model.export_to_openfoam()


#==============================================================================*
# Pre-processing

scalingVector = [100, 1, 1]

ffn.run_preprocessing(model=model)
model.plot_mesh(region=fuelMesh, scalingVector=scalingVector, show_edges=True)
model.plot_mesh(region=fuelMesh, scalingVector=scalingVector, show_edges=True, normal='y')


#==============================================================================*
# Simulation

ffn.run(model=model)


#==============================================================================*
# Post-processing

lastTime = model.get_time_steps()[-1]

for fieldName in ['T', 'Bu']:
    model.plot_mesh(
        region=fuelMesh,
        time=lastTime,
        fieldName=fieldName,
        cmap="RdBu_r",
        normal='y',
        scalingVector=scalingVector
    )


# Extract data from probe

T, points = probes.read_from_case(
    startTime=0,
    fieldName='T'
)
Bu, points = probes.read_from_case(
    startTime=0,
    fieldName='Bu'
)

time = [t/24/3600 for t in list(T.keys())]
Tlist = [val[0] for val in T.values()]
Bulist = [val[0]/1000 for val in Bu.values()]


fig, axT = plt.subplots(figsize=(5, 4), dpi=200)
axBu = axT.twinx()

color = 'tab:red'
axT.plot(time, Tlist, label= "Temperature", color=color)
axT.tick_params(axis='y', labelcolor=color)
axT.set_xlabel('Time, days')
axT.set_ylabel('Temperature, K', color=color)

color = 'tab:blue'
axBu.plot(time, Bulist, '-.', label= "Burnup", color=color)
axBu.tick_params(axis='y', labelcolor=color)
axBu.set_ylabel('Burnup, MWd/kg', color=color)

fig.legend()

fig.tight_layout()
fig.savefig("fig_results_T_Bu.png")


#==============================================================================*
