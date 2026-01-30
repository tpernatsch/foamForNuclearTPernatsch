"""

"""

import numpy as np
import foamForNuclear as ffn
import foamForNuclear.mesh as mesh
import foamForNuclear.boundaryConditions as bc


isSimpleMesh = True


# --- Basis
kg, s, m, J, mol, K, rad = 1, 1, 1, 1, 1, 1, 1

cm = 1e-2 * m
inch = 2.54 * cm



fuelElementPitch    =  0.754 * inch
coreHeight          = 52.000 * inch
coreChannelPitch    =  0.174 * inch
coreChannelDiameter =  0.100 * inch

volumeFractionStructure = 1 - 0.299544

lengthUpperPlenum   = 20 * cm
lengthNozzleChamber = (62.694 - 53.5) * inch

nTotElement = 1500
nElementPerAssembly = 6

nzFuelElement   = 15
nzNozzleChamber = 5
nzUpperPlenum   = 10


lattice = """0 F F
F C F
F F 0"""


#==============================================================================*
# Mesh

def createFuelAssemblyMesh():
    nMesh = mesh.BlockMesh(region='neutroRegion')

    nMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: nMesh.createHexagonPrismAlongZ(
            name="fuelElement",
            zmin=-coreHeight/2,
            zmax=coreHeight/2,
            pitch=fuelElementPitch,
            x=x, y=y,
            nr=1,
            nt=1,
            nz=nzFuelElement,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=3, ny=3, pitch=fuelElementPitch,
        elementsToPlace=['F'],
        # isMergePatches=True
    )
    nMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: nMesh.createHexagonPrismAlongZ(
            name="centralUnloadedFuelElement",
            zmin=-coreHeight/2,
            zmax=coreHeight/2,
            pitch=fuelElementPitch,
            x=x, y=y,
            nr=1,
            nt=1,
            nz=nzFuelElement,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=3, ny=3, pitch=fuelElementPitch,
        elementsToPlace=['C'],
        # isMergePatches=True
    )

    nMesh.addMergePatchPairs()

    nMesh.mergePatchesWithName(name="top", regex=".*Top_.*")
    nMesh.mergePatchesWithName(name="bottom", regex=".*Bottom_.*")
    nMesh.mergePatchesWithName(name="wall", regex="fuelElementWall.*")



    thMesh = mesh.BlockMesh(region='fluidRegion')

    thMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: thMesh.createHexagonPrismAlongZ(
            name="fuelElement",
            zmin=-coreHeight/2,
            zmax=coreHeight/2,
            pitch=fuelElementPitch,
            x=x, y=y,
            nr=1,
            nt=1,
            nz=nzFuelElement,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=3, ny=3, pitch=fuelElementPitch,
        elementsToPlace=['F'],
        # isMergePatches=True
    )
    thMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: thMesh.createHexagonPrismAlongZ(
            name="centralUnloadedFuelElement",
            zmin=-coreHeight/2,
            zmax=coreHeight/2,
            pitch=fuelElementPitch,
            x=x, y=y,
            nr=1,
            nt=1,
            nz=nzFuelElement,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=3, ny=3, pitch=fuelElementPitch,
        elementsToPlace=['C'],
        # isMergePatches=True
    )
    thMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: thMesh.createHexagonPrismAlongZ(
            name="upperPlenum",
            zmin=-coreHeight/2 - lengthNozzleChamber,
            zmax=-coreHeight/2,
            pitch=fuelElementPitch,
            x=x, y=y,
            nr=1,
            nt=1,
            nz=nzUpperPlenum,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=3, ny=3, pitch=fuelElementPitch,
        elementsToPlace=['F', 'C'],
    )
    thMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: thMesh.createHexagonPrismAlongZ(
            name="nozzleChamber",
            zmin=coreHeight/2,
            zmax=coreHeight/2 + lengthNozzleChamber,
            pitch=fuelElementPitch,
            x=x, y=y,
            nr=1,
            nt=1,
            nz=nzNozzleChamber,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=3, ny=3, pitch=fuelElementPitch,
        elementsToPlace=['F', 'C'],
    )

    thMesh.addMergePatchPairs(
        includeFacename=['Top', 'Bottom'],
        excludeFacename=['Wall']
    )
    thMesh.addMergePatchPairs(
        includeFacename=['nozzleChamberWall', 'upperPlenumWall'],
        excludeFacename=['Top', 'Bottom']
    )
    bafflesFaces = thMesh.addBaffles(
        includeFacename=['fuelElementWall', 'centralUnloadedFuelElementWall'],
        excludeFacename=['Top', 'Bottom']
    )

    externalWalls = thMesh.getStandaloneFaces(includeFacename=['Wall'])

    thMesh.mergePatchesWithName(name="wall", includeFacename=[face.name for face in externalWalls], patchType="wall")
    thMesh.mergePatchesWithName(name="outlet", includeFacename=["Top_.*"])
    thMesh.mergePatchesWithName(name="inlet", includeFacename=["Bottom_"])

    return(nMesh, thMesh)


def create1DMesh():
    sqrSide = np.sqrt(np.sqrt(3)/2) * fuelElementPitch

    nMesh = mesh.BlockMesh(region='neutroRegion')

    nMesh.createCube(
        'fuelElement',
        lowX=-sqrSide/2, highX=sqrSide/2,
        lowY=-sqrSide/2, highY=sqrSide/2,
        lowZ=-coreHeight/2, highZ=coreHeight/2,
        nx=1, ny=1, nz=nzFuelElement,
        isAddBoundaryConditions=True
    )

    nMesh.mergePatchesWithName(name="wall", includeFacename=["Wall"], patchType="wall")
    nMesh.mergePatchesWithName(name="top", includeFacename=["Top_"])
    nMesh.mergePatchesWithName(name="bottom", includeFacename=["Bottom_"])

    thMesh = mesh.BlockMesh(region='fluidRegion')

    thMesh.createCube(
        'fuelElement',
        lowX=-sqrSide/2, highX=sqrSide/2,
        lowY=-sqrSide/2, highY=sqrSide/2,
        lowZ=-coreHeight/2, highZ=coreHeight/2,
        nx=1, ny=1, nz=nzFuelElement,
        isAddBoundaryConditions=True
    )

    thMesh.mergePatchesWithName(name="wall", includeFacename=["Wall"], patchType="wall")
    thMesh.mergePatchesWithName(name="outlet", includeFacename=["Top_"])
    thMesh.mergePatchesWithName(name="inlet", includeFacename=["Bottom_"])

    return(nMesh, thMesh)


if isSimpleMesh:
    nMesh, thMesh = create1DMesh()
else:
    nMesh, thMesh = createFuelAssemblyMesh()


#==============================================================================*
# Time folder

inletT = 300

timeFolder0 = ffn.TimeFolder(0)

defaultFlux = ffn.Field("defaultFlux", region=nMesh.region)
defaultFlux.dimensions = ffn.Dimension(default='neutronFlux')
defaultFlux.internalField = 1
defaultFlux.set_boundary_condition("top", bc.FixedValue(0))
defaultFlux.set_boundary_condition("bottom", bc.FixedValue(0))
defaultFlux.set_boundary_condition("wall", bc.ZeroGradient())

TFuel = ffn.Field("TFuel", region=nMesh.region)
TFuel.dimensions = ffn.Dimension(default='T')
TFuel.internalField = inletT
TFuel.set_boundary_condition("top", bc.ZeroGradient())
TFuel.set_boundary_condition("bottom", bc.ZeroGradient())
TFuel.set_boundary_condition("wall", bc.ZeroGradient())

rhoCool = ffn.Field("rhoCool", region=nMesh.region)
rhoCool.dimensions = ffn.Dimension(mass=1, length=-3)
rhoCool.internalField = 11.1
rhoCool.set_boundary_condition("top", bc.ZeroGradient())
rhoCool.set_boundary_condition("bottom", bc.ZeroGradient())
rhoCool.set_boundary_condition("wall", bc.ZeroGradient())

alphat = ffn.Field("alphat", region=thMesh.region)
alphat.dimensions = ffn.Dimension(mass=1, length=-1, time=-1)
alphat.internalField = 0
alphat.set_boundary_condition("inlet", bc.FixedValue(0))
alphat.set_boundary_condition("outlet", bc.FixedValue(0))
alphat.set_boundary_condition("wall", bc.ZeroGradient())
if not isSimpleMesh:
    alphat.set_boundary_condition("baffle0", bc.ZeroGradient())
    alphat.set_boundary_condition("baffle1", bc.ZeroGradient())

epsilon = ffn.Field("epsilon", region=thMesh.region)
epsilon.dimensions = ffn.Dimension(length=2, time=-3)
epsilon.internalField = 0.01155
epsilon.set_boundary_condition("inlet", bc.FixedValue(epsilon.internalField))
epsilon.set_boundary_condition("outlet", bc.ZeroGradient())
epsilon.set_boundary_condition("wall", bc.ZeroGradient())
if not isSimpleMesh:
    epsilon.set_boundary_condition("baffle0", bc.ZeroGradient())
    epsilon.set_boundary_condition("baffle1", bc.ZeroGradient())

k = ffn.Field("k", region=thMesh.region)
k.dimensions = ffn.Dimension(length=2, time=-2)
k.internalField = 0.02675
k.set_boundary_condition("inlet", bc.FixedValue(k.internalField))
k.set_boundary_condition("outlet", bc.ZeroGradient())
k.set_boundary_condition("wall", bc.ZeroGradient())
if not isSimpleMesh:
    k.set_boundary_condition("baffle0", bc.ZeroGradient())
    k.set_boundary_condition("baffle1", bc.ZeroGradient())

nut = ffn.Field("nut", region=thMesh.region)
nut.dimensions = ffn.Dimension(length=2, time=-1)
nut.internalField = 0
nut.set_boundary_condition("inlet", bc.FixedValue(0))
nut.set_boundary_condition("outlet", bc.FixedValue(0))
nut.set_boundary_condition("wall", bc.ZeroGradient())
if not isSimpleMesh:
    nut.set_boundary_condition("baffle0", bc.ZeroGradient())
    nut.set_boundary_condition("baffle1", bc.ZeroGradient())

p = ffn.Field("p", region=thMesh.region)
p.dimensions = ffn.Dimension(default='p')
p.internalField = 3.7e6
p.set_boundary_condition("inlet", bc.ZeroGradient())
# p.set_boundary_condition("outlet", bc.FixedValue(3.427e6))
p.set_boundary_condition("outlet", bc.CoupledUniformExternalValue(3.427e6, inputName="gfpNozzleChamber_in", initValue=3.427e6))
p.set_boundary_condition("wall", bc.ZeroGradient())
if not isSimpleMesh:
    p.set_boundary_condition("baffle0", bc.ZeroGradient())
    p.set_boundary_condition("baffle1", bc.ZeroGradient())

p_rgh = ffn.Field("p_rgh", region=thMesh.region)
p_rgh.dimensions = ffn.Dimension(default='p')
p_rgh.internalField = 3.7e6
p_rgh.set_boundary_condition("inlet", bc.ZeroGradient())
# p_rgh.set_boundary_condition("outlet", bc.FixedValue(3.427e6))
p_rgh.set_boundary_condition("outlet", bc.CoupledUniformExternalValue(3.427e6, inputName="gfpNozzleChamber_in", initValue=3.427e6))
p_rgh.set_boundary_condition("wall", bc.ZeroGradient())
if not isSimpleMesh:
    p_rgh.set_boundary_condition("baffle0", bc.ZeroGradient())
    p_rgh.set_boundary_condition("baffle1", bc.ZeroGradient())

T = ffn.Field("T", region=thMesh.region)
T.dimensions = ffn.Dimension(default='T')
T.internalField = inletT
T.set_boundary_condition("inlet", bc.CoupledUniformExternalValue(T.internalField, inputName="gfT_in", initValue=T.internalField))
# T.set_boundary_condition("inlet", bc.FixedValue(T.internalField))
T.set_boundary_condition("outlet", bc.ZeroGradient())
T.set_boundary_condition("wall", bc.ZeroGradient())
if not isSimpleMesh:
    T.set_boundary_condition("baffle0", bc.ZeroGradient())
    T.set_boundary_condition("baffle1", bc.ZeroGradient())

U = ffn.Field("U", region=thMesh.region)
U.dimensions = ffn.Dimension(default='U')
U.internalField = ffn.Vector(0, 0, 4.74413)
# U.set_boundary_condition('inlet', bc.FixedValue(U.internalField))
mFlowInit = 30
if isSimpleMesh:
    mFlowInit *= 1/nTotElement
else:
    mFlowInit *= nElementPerAssembly/nTotElement
U.set_boundary_condition("inlet", bc.CoupledFlowRateInletVelocity(
    value=U.internalField,
    massFlowRate='gfMassFlowInlet_in',
    mFlowInit=mFlowInit,
    # rho="rho",
    # rho="alphaRhoPhi",
    rho="thermo:rho",
    # rhoInlet=1,
    extrapolateProfile=False
))
U.set_boundary_condition("outlet", bc.ZeroGradient())
U.set_boundary_condition("wall", bc.Slip())
if not isSimpleMesh:
    U.set_boundary_condition("baffle0", bc.Slip())
    U.set_boundary_condition("baffle1", bc.Slip())

timeFolder0.append(defaultFlux)
timeFolder0.append(TFuel)
timeFolder0.append(rhoCool)
timeFolder0.append(alphat)
timeFolder0.append(epsilon)
timeFolder0.append(k)
timeFolder0.append(nut)
timeFolder0.append(p)
timeFolder0.append(p_rgh)
timeFolder0.append(T)
timeFolder0.append(U)




#==============================================================================*
# Neutronics solver

# 937 MW, 6 fuel elements per assembly
# ~1500 fuel elements in the core
# power = 937e6
power = 10e3
if isSimpleMesh:
    power *= 1.0/nTotElement
else:
    power *= nElementPerAssembly/nTotElement

neutronicsSolver = ffn.NeutronicsSolver(
    region=nMesh.region,
    solver="diffusionNeutronics",
    mesh=nMesh,
    power=power

)

neutronicsSolver.neutronTransportOptions.maxNeutronIterations = 50


#==============================================================================*
# Thermal-hydraulics solver

thSolver = ffn.ThermalHydraulicsSolver(
    region=thMesh.region,
    solver="onePhase",
    mesh=thMesh,
    removeBaffles=True,
    isSetFvSolutionToDefault=True
)
thSolver.thermophysicalProperties = ffn.thermophysicalProperty.Hydrogen()
# thSolver.thermophysicalProperties = ffn.thermophysicalProperty.HydrogenPengRobinsonGas()
# thSolver.thermophysicalProperties = ffn.thermophysicalProperty.HydrogenPerfectGas()

thSolver.phaseProperties.pMin = 20e5

thSolver.fvSchemes.divSchemes['default'] = 'none'
thSolver.fvSchemes.divSchemes["div\(phi.*,U.*\)"] = "Gauss upwind"
thSolver.fvSchemes.divSchemes["div(alphaRhoPhi,U)"] = "Gauss upwind"
thSolver.fvSchemes.divSchemes["div(alphaRhoPhiNu,U)"] = "Gauss linear"
thSolver.fvSchemes.divSchemes["div(alphaRhoPhi,K)"] = "Gauss upwind"
thSolver.fvSchemes.divSchemes["div\(alphaRhoPhi.*,k.*\)"] = "Gauss upwind"
thSolver.fvSchemes.divSchemes["div\(alphaRhoPhi.*,epsilon.*\)"] = "Gauss upwind"
thSolver.fvSchemes.divSchemes["div\(alphaRhoPhi.*,(h|e).*\)"] = "Gauss upwind"

thSolution = ffn.fvSolution()

thSolution.append('p_rgh', ffn.fvSolutionSolver(
    solver='GAMG', smoother='DIC', tolerance=1e-5, relTol=1e-2
))
thSolution.append('p_rghFinal', ffn.fvSolutionSolver(
    solver='GAMG', smoother='DIC', tolerance=1e-5, relTol=1e-2
))
thSolution.append('"e.*"', ffn.fvSolutionSolver(
    solver='smoothSolver',
    smoother='symGaussSeidel',
    tolerance=1e-5, relTol=1e-2, minIter=1
))
thSolution.append('"h.*"', ffn.fvSolutionSolver(
    solver='smoothSolver',
    smoother='symGaussSeidel',
    tolerance=1e-5, relTol=1e-2, minIter=1
))
thSolution.append('".*"', ffn.fvSolutionSolver(
    solver='PBiCGStab',
    preconditioner='diagonal',
    tolerance=1e-5, relTol=1e-2
))

thSolver.fvSolution = thSolution

thSolver.pimpleOptions.nCorrectors = 2
thSolver.pimpleOptions.nOuterCorrectors = 1
thSolver.pimpleOptions.nNonOrthogonalCorrectors = 0
thSolver.pimpleOptions.momentumMode = "faceCentered"
thSolver.pimpleOptions.correctUntilConvergence = None
thSolver.pimpleOptions.porousInterfaceSharpness = None
thSolver.pimpleOptions.minMagU = None
thSolver.pimpleOptions.minNOuterCorrectors = None


thSolver.turbulenceProperties.simulationType = 'RAS'
thSolver.turbulenceProperties.RASoptions.RASModel = 'porousKEpsilon'
thSolver.turbulenceProperties.RASoptions.turbulence = True
thSolver.turbulenceProperties.RASoptions.printCoeffs = True
thSolver.turbulenceProperties.porousKEpsilonProperties.append(ffn.PorousKEpsilonPropertiesPerZone(
    zones=['fuelElement'],
    convergenceLength=0.1,
    turbulenceIntensityCoeff=0.16,
    turbulenceIntensityExp=-0.125,
    turbulenceLengthScaleCoeff=0.07,
    DhStruct=19.15e-3
))

fuelElementProperty = ffn.StructureProperty(
    zones=['fuelElement'],
    volumeFraction=volumeFractionStructure,
    Dh=coreChannelDiameter
)
fuelElementProperty.add_power_model(
    ffn.LumpedNuclearStructure(
        volumetricArea=477.305494,
        powerDensity=0,
        nodeFuel=0,
        nodeClad=0,
        nodeMatrix=0,
        kappaMatrix=50,
        heatConductances=ffn.List([2.451162e+07, 1.149520e+07]),
        rhoCp=ffn.List([3671370]),
        volumeFractions=ffn.List([1]),
        powerFractions=ffn.List([1]),
        T0=T.internalField
    )
)
thSolver.add_structure_property(fuelElementProperty)

thSolver.add_drag_model(ffn.Colebrook(
    coeff=1.2776,
    const=-0.406,
    exp=-2.246,
    zones=['fuelElement']
))
thSolver.add_heat_transfer_model(ffn.NusseltReynoldsPrandtlPower(
    const=0,
    coeff=0.023,
    expRe=0.8,
    expPr=0.4,
    expTc=-0.3,
    zones=['fuelElement']
))

if not isSimpleMesh:
    thSolver.turbulenceProperties.porousKEpsilonProperties.append(ffn.PorousKEpsilonPropertiesPerZone(
        zones=['centralUnloadedFuelElement'],
        convergenceLength=0.1,
        turbulenceIntensityCoeff=0.16,
        turbulenceIntensityExp=-0.125,
        turbulenceLengthScaleCoeff=0.07,
        DhStruct=19.15e-3
    ))
    centralUnloadedStructureProp = ffn.StructureProperty(
        zones=['centralUnloadedFuelElement'],
        volumeFraction=0.983257,
        Dh=0.000589/10,
    )
    centralUnloadedStructureProp.add_passive_structure(
        volumetricArea=113.648950,
        rho=8190,
        Cp=500,
        T=T.internalField
    )
    thSolver.add_structure_property(centralUnloadedStructureProp)
    thSolver.add_drag_model(ffn.ReynoldsPower(
        coeff=0.0625,
        exp=-0.32,
        zones=['centralUnloadedFuelElement']
    ))
    thSolver.add_heat_transfer_model(ffn.NusseltReynoldsPrandtlPower(
        const=0,
        coeff=0.023,
        expRe=0.8,
        expPr=0.4,
        expTc=-0.3,
        zones=['centralUnloadedFuelElement']
    ))


#==============================================================================*
# Solvers

solvers = ffn.Solvers([neutronicsSolver, thSolver])

#==============================================================================*
# Coupling

coupling = ffn.Coupling(solvers)

coupling.add_field_transfer(neutronicsSolver, thSolver, "powerDensity", "powerDensityStructure")
coupling.add_field_transfer(neutronicsSolver, thSolver, "secondaryPowerDensity", "powerDensityLiquid")

coupling.add_field_transfer(thSolver, neutronicsSolver, "T", "TCool")
coupling.add_field_transfer(thSolver, neutronicsSolver, "thermo:rho", "rhoCool")
coupling.add_field_transfer(thSolver, neutronicsSolver, "T.fuelAvForNeutronics", "TFuel")
coupling.add_field_transfer(thSolver, neutronicsSolver, "T.cladAvForNeutronics", "TClad")
coupling.add_field_transfer(thSolver, neutronicsSolver, "T.passiveStructure", "TStructMech")


#==============================================================================*
# FMU
# The FMU needs to be linked even during steady-state, this improves stability
# of the results especially if PID controllers are present

massFlowScalingFactor = (nTotElement if isSimpleMesh else nTotElement/nElementPerAssembly) / (1 - volumeFractionStructure)

# Coupling interface
externalCouplingDict = ffn.ExternalCouplingDict()

TnozzleChamberFMU = ffn.ExtPatch(
    name="TnozzleChamberFMU",
    outputName="gfTnozzleChamber_out",
    fieldName="T",
    region="fluidRegion",
    patchName="outlet",
    initValue=inletT
)
pInletFMU = ffn.ExtPatch(
    name="pInletFMU",
    outputName="gfpInlet_out",
    fieldName="p",
    region="fluidRegion",
    patchName="inlet",
    initValue=12e5
)
integratePowerFMU = ffn.FieldIntegralToFMU(
    name="fieldIntegralToFMU",
    nameFMU="gfPower_out",
    fieldName="powerDensity",
    region="neutroRegion",
    cellZone="fuelElement",
    initValue=power
)
massFlowOutletFMU = ffn.MassFlowToFMU(
    name='massFlowOutletFMU',
    nameFMU="gfMassFlowOutlet_out",
    regionType="patch",
    regionName="outlet",
    region="fluidRegion",
    alphaRhoPhiName="alphaRhoPhi",
    scaleFactor=massFlowScalingFactor,
    initValue=30
)
externalCouplingDict.append(TnozzleChamberFMU)
externalCouplingDict.append(pInletFMU)
externalCouplingDict.append(integratePowerFMU)
externalCouplingDict.append(massFlowOutletFMU)


#==============================================================================*
# Function Objects

massFlowInletFO = ffn.MassFlow(
    "mFlowInlet",
    log=True,
    writeFields=False,
    writeControl="adjustableRunTime",
    writeInterval=0.01,
    region="fluidRegion",
    regionType="patch",
    regionName="inlet",
    scaleFactor=massFlowScalingFactor
)
massFlowOutletFO = ffn.MassFlow(
    "mFlowOutlet",
    log=True,
    writeFields=False,
    writeControl="adjustableRunTime",
    writeInterval=0.01,
    region="fluidRegion",
    regionType="patch",
    regionName="outlet",
    scaleFactor=massFlowScalingFactor
)
TBulkInletFO = ffn.TBulk(
    "TInlet",
    log=True,
    writeFields=False,
    writeControl="adjustableRunTime",
    writeInterval=0.01,
    region="fluidRegion",
    regionType="patch",
    regionName="inlet"
)
TBulkOutletFO = ffn.TBulk(
    "TOutlet",
    log=True,
    writeFields=False,
    writeControl="adjustableRunTime",
    writeInterval=0.01,
    region="fluidRegion",
    regionType="patch",
    regionName="outlet"
)


#==============================================================================*
# Model

model = ffn.Model(
    solvers=solvers,
    coupling=coupling,
    timeFolders=[timeFolder0],
    externalCouplingDict=externalCouplingDict
)

model.add_function_object(massFlowInletFO)
model.add_function_object(massFlowOutletFO)
model.add_function_object(TBulkInletFO)
model.add_function_object(TBulkOutletFO)

settings = model.settings

settings.application = "GeN-Foam"
settings.startTime = 0
settings.endTime = 10
settings.deltaT = 1e-6
settings.startFrom = 'latestTime'
settings.writeControl = 'adjustableRunTime'
settings.writeInterval = settings.endTime
settings.adjustTimeStep = True # Not possible to do adaptative time step with FMU !?
settings.runTimeModifiable = True
settings.maxDeltaT = 0.1
settings.maxCo = 1
settings.maxPowerVariation = 0.0001
