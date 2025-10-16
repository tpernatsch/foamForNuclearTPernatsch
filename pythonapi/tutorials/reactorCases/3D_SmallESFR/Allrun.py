"""

"""
#==============================================================================*
# Imports

from copy import copy
import numpy as np
import time
import matplotlib.pyplot as plt

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh


#==============================================================================*
# nMesh

# Use https://foam-for-nuclear.gitlab.io/honeycomb/

lattice = """0 0 0 0 0 0 0 R R R R R R R R
 0 0 0 0 0 0 R O C O O O C O R
  0 0 0 0 0 R C O O O O O O C R
   0 0 0 0 R O O I I I I I O O R
    0 0 0 R O O I C I I C I O O R
     0 0 R O O I I I I I I I O O R
      0 R C O I I I I I I I I O C R
       R O O I C I I I I I C I O O R
        R C O I I I I I I I I O C R 0
         R O O I I I I I I I O O R 0 0
          R O O I C I I C I O O R 0 0 0
           R O O I I I I I O O R 0 0 0 0
            R C O O O O O O C R 0 0 0 0 0
             R O C O O O C O R 0 0 0 0 0 0
              R R R R R R R R 0 0 0 0 0 0 0"""



def createMesh(region: str, isFluidMesh: bool, gap: float=0):
    lc = 0.5
    edge = 0.12171 # 0.11391
    # gap = 0 # 0.009

    coreHeight = 1
    lowReflHeight = 0.3
    lowFGPHeight = 0.36
    diagridHeight = 0.3
    upperFGPHeight = 0.11
    upperReflHeight = 0.38

    coreNodes = 10
    lowReflNodes = 3
    lowFGPNodes = 3
    diagridNodes = 3
    upperFGPNodes = 2
    upperReflNodes = 3

    totalHeight =  lowFGPHeight + lowReflHeight + coreHeight + upperFGPHeight + upperReflHeight

    halfEdge = edge/2
    flatToFlat = edge * np.sqrt(3)
    halfFlatToFlat = flatToFlat/2
    pitch = flatToFlat + gap

    edgeExt = edge + gap/np.sqrt(3)
    halfEdgeExt = edgeExt/2
    flatToFlatExt = edgeExt * np.sqrt(3)
    halfFlatToFlatExt = flatToFlatExt/2


    latticeNX = len(lattice.split("\n"))
    latticeNY = latticeNX
    nr = 1
    nt = 1

    blockMesh = mesh.BlockMesh(region=region)

    blockMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: blockMesh.createHexagonPrismAlongZ(
            name="innerCore",
            zmin=-coreHeight/2,
            zmax=coreHeight/2,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=coreNodes,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['I'],
        # isMergePatches=isMergePatches
    )
    blockMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: blockMesh.createHexagonPrismAlongZ(
            name="outerCore",
            zmin=-coreHeight/2,
            zmax=coreHeight/2,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=coreNodes,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['O'],
        # isMergePatches=isMergePatches
    )
    blockMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: blockMesh.createHexagonPrismAlongZ(
            name="upperGasPlenum",
            zmin=coreHeight/2,
            zmax=coreHeight/2+upperFGPHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=upperFGPNodes,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['I', 'O'],
        # isMergePatches=isMergePatches
    )

    # """
    blockMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: blockMesh.createHexagonPrismAlongZ(
            name="upperReflector",
            zmin=coreHeight/2+upperFGPHeight,
            zmax=coreHeight/2+upperFGPHeight+upperReflHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=upperReflNodes,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['I', 'O'],
        # isMergePatches=isMergePatches
    )
    blockMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: blockMesh.createHexagonPrismAlongZ(
            name="lowerGasPlenum",
            zmin=-coreHeight/2-lowFGPHeight,
            zmax=-coreHeight/2,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=lowFGPNodes,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['I', 'O'],
        # isMergePatches=isMergePatches
    )

    blockMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: blockMesh.createHexagonPrismAlongZ(
            name="lowerReflector",
            zmin=-coreHeight/2-lowFGPHeight-lowReflHeight,
            zmax=-coreHeight/2-lowFGPHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=lowReflNodes,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['I', 'O'],
        # isMergePatches=isMergePatches
    )

    # Radial reflector
    blockMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: blockMesh.createHexagonPrismAlongZ(
            name="radialReflector",
            zmin=-coreHeight/2,
            zmax=coreHeight/2,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=coreNodes,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['R'],
        # isMergePatches=isMergePatches
    )
    blockMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: blockMesh.createHexagonPrismAlongZ(
            name="radialReflector",
            zmin=coreHeight/2,
            zmax=coreHeight/2+upperFGPHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=upperFGPNodes,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['R'],
        # isMergePatches=isMergePatches
    )
    blockMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: blockMesh.createHexagonPrismAlongZ(
            name="radialReflector",
            zmin=coreHeight/2+upperFGPHeight,
            zmax=coreHeight/2+upperFGPHeight+upperReflHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=upperReflNodes,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['R'],
        # isMergePatches=isMergePatches
    )
    blockMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: blockMesh.createHexagonPrismAlongZ(
            name="radialReflector",
            zmin=-coreHeight/2-lowFGPHeight,
            zmax=-coreHeight/2,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=lowFGPNodes,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['R'],
        # isMergePatches=isMergePatches
    )
    blockMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: blockMesh.createHexagonPrismAlongZ(
            name="radialReflector",
            zmin=-coreHeight/2-lowFGPHeight-lowReflHeight,
            zmax=-coreHeight/2-lowFGPHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=lowReflNodes,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['R'],
        # isMergePatches=isMergePatches
    )

    blockMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: blockMesh.createHexagonPrismAlongZ(
            name="diagrid",
            zmin=-coreHeight/2-lowFGPHeight-lowReflHeight-diagridHeight,
            zmax=-coreHeight/2-lowFGPHeight-lowReflHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=diagridNodes,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['I', 'O', 'R', 'C'],
        # isMergePatches=isMergePatches
    )


    # Follower
    blockMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: blockMesh.createHexagonPrismAlongZ(
            name="follower",
            zmin=-coreHeight/2-lowFGPHeight-lowReflHeight,
            zmax=-coreHeight/2-lowFGPHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=lowReflNodes,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['C'],
        # isMergePatches=isMergePatches
    )
    blockMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: blockMesh.createHexagonPrismAlongZ(
            name="follower",
            zmin=-coreHeight/2-lowFGPHeight,
            zmax=-coreHeight/2,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=lowFGPNodes,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['C'],
        # isMergePatches=isMergePatches
    )
    blockMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: blockMesh.createHexagonPrismAlongZ(
            name="follower",
            zmin=-coreHeight/2,
            zmax=coreHeight/2,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=coreNodes,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['C'],
        # isMergePatches=isMergePatches
    )

    # Control rod
    blockMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: blockMesh.createHexagonPrismAlongZ(
            name="controlRod",
            zmin=coreHeight/2,
            zmax=coreHeight/2+upperFGPHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=upperFGPNodes,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['C'],
        # isMergePatches=isMergePatches
    )
    blockMesh.latticePlacement(
        funcElementGenerator=lambda name, x, y: blockMesh.createHexagonPrismAlongZ(
            name="controlRod",
            zmin=coreHeight/2+upperFGPHeight,
            zmax=coreHeight/2+upperFGPHeight+upperReflHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=upperReflNodes,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['C'],
        # isMergePatches=isMergePatches
    )
    # """

    blockMesh.addMergePatchPairs(
        includeFacename=['Top', 'Bottom'],
        excludeFacename=['Wall']
    )
    if (isFluidMesh):
        bafflesFaces = blockMesh.addBaffles(
            includeFacename=['Wall'],
            excludeFacename=['Top', 'Bottom']
        )
    else:
        blockMesh.addMergePatchPairs(
            includeFacename=['Wall'],
            excludeFacename=['Top', 'Bottom']
        )


    externalWalls = blockMesh.getStandaloneFaces(
        includeFacename=['Wall'],
            excludeFacename=['Top', 'Bottom']
    )

    blockMesh.mergePatchesWithName(
        name="wall",
        includeFacename=[face.name for face in externalWalls],
        patchType="wall"
    )
    blockMesh.mergePatchesWithName(
        name="bottom",
        includeFacename=["diagridBottom_"]
    )
    blockMesh.mergePatchesWithName(
        name="top",
        includeFacename=["upperReflectorTop_", "controlRodTop_", "radialReflectorTop_"]
    )

    return(blockMesh)


# nMesh = createMesh(region="neutroRegion", isFluidMesh=False)
# thMesh = createMesh(region="fluidRegion", isFluidMesh=True) #, gap=0.009)
nMesh = mesh.PolyMesh(region="neutroRegion", srcpath="./meshes/polyMeshNeutro")
thMesh = mesh.PolyMesh(region="fluidRegion", srcpath="./meshes/polyMeshFluid")
tmMesh = mesh.PolyMesh(region="thermoMechanicalRegion", srcpath="./meshes/polyMeshMeca")



#==============================================================================*
# Fields

inletTemperature = 668

timeFolder0 = ffn.TimeFolder(time=0)

# Neutronics
defaultFlux = ffn.Field("defaultFlux", region=nMesh.region)
defaultFlux.dimensions = ffn.Dimension(default='flux')
defaultFlux.internalField = 1e21
defaultFlux.set_boundary_condition('wall', bc.FixedValue(0))
defaultFlux.set_boundary_condition('top', bc.FixedValue(0))
defaultFlux.set_boundary_condition('bottom', bc.FixedValue(0))

# Fluid
T = ffn.Field("T", region=thMesh.region)
T.dimensions = ffn.Dimension(default='T')
T.internalField = inletTemperature
T.set_boundary_condition('wall', bc.ZeroGradient())
T.set_boundary_condition('baffle0', bc.ZeroGradient())
T.set_boundary_condition('baffle1', bc.ZeroGradient())
T.set_boundary_condition('top', bc.ZeroGradient())
T.set_boundary_condition('bottom', bc.FixedValue(inletTemperature))

U = ffn.Field("U", region=thMesh.region)
U.dimensions = ffn.Dimension(default='U')
U.internalField = ffn.Vector(0, 0, 4.559114016)
U.set_boundary_condition('wall', bc.Slip())
U.set_boundary_condition('baffle0', bc.Slip())
U.set_boundary_condition('baffle1', bc.Slip())
U.set_boundary_condition('top', bc.ZeroGradient())
U.set_boundary_condition('bottom', bc.FixedValue(U.internalField))

p = ffn.Field("p", region=thMesh.region)
p.dimensions = ffn.Dimension(default='p')
p.internalField = 100000
p.set_boundary_condition('wall', bc.ZeroGradient())
p.set_boundary_condition('baffle0', bc.Calculated(0))
p.set_boundary_condition('baffle1', bc.Calculated(0))
p.set_boundary_condition('top', bc.FixedValue(p.internalField))
p.set_boundary_condition('bottom', bc.ZeroGradient())

p_rgh = ffn.Field("p_rgh", region=thMesh.region)
p_rgh.dimensions = ffn.Dimension(default='p')
p_rgh.internalField = p.internalField
p_rgh.set_boundary_condition('wall', bc.ZeroGradient())
p_rgh.set_boundary_condition('baffle0', bc.ZeroGradient())
p_rgh.set_boundary_condition('baffle1', bc.ZeroGradient())
p_rgh.set_boundary_condition('top', bc.FixedValue(p_rgh.internalField))
p_rgh.set_boundary_condition('bottom', bc.FixedFluxExtrapolatedPressure(value=p_rgh.internalField, gradient=0))

# Thermomechanics
Tmech = ffn.Field("T", region=tmMesh.region)
Tmech.dimensions = ffn.Dimension(default="T")
Tmech.internalField = inletTemperature
Tmech.set_boundary_condition('zeroDispl', bc.ZeroGradient())
Tmech.set_boundary_condition('defaultFaces', bc.ZeroGradient())
Tmech.set_boundary_condition('bottom', bc.FixedValue(inletTemperature))
Tmech.set_boundary_condition('top', bc.ZeroGradient())
Tmech.set_boundary_condition('topCR', bc.ZeroGradient())
Tmech.set_boundary_condition('bottomCR', bc.FixedValue(inletTemperature))

D = ffn.Field("D", region=tmMesh.region)
D.dimensions = ffn.Dimension(length=1)
D.internalField = ffn.Vector(0, 0, 0)
zeroTraction = bc.TractionDisplacement(
    traction=D.internalField, pressure=0, value=D.internalField
)
D.set_boundary_condition('zeroDispl', bc.FixedValue(D.internalField))
D.set_boundary_condition('defaultFaces', zeroTraction)
D.set_boundary_condition('bottom', zeroTraction)
D.set_boundary_condition('top', zeroTraction)
D.set_boundary_condition('topCR', zeroTraction)
D.set_boundary_condition('bottomCR', zeroTraction)

fuelDisp = ffn.Field("fuelDisp", region=tmMesh.region)
fuelDisp.dimensions = ffn.Dimension(length=1)
fuelDisp.internalField = 0
fuelDisp.set_boundary_condition('zeroDispl', bc.FixedValue(0))
fuelDisp.set_boundary_condition('defaultFaces', bc.ZeroGradient())
fuelDisp.set_boundary_condition('bottom', bc.FixedValue(0))
fuelDisp.set_boundary_condition('top', bc.ZeroGradient())
fuelDisp.set_boundary_condition('topCR', bc.ZeroGradient())
fuelDisp.set_boundary_condition('bottomCR', bc.FixedValue(0))

CRDisp = ffn.Field("CRDisp", region=tmMesh.region)
CRDisp.dimensions = ffn.Dimension(length=1)
CRDisp.internalField = 0
CRDisp.set_boundary_condition('zeroDispl', bc.ZeroGradient())
CRDisp.set_boundary_condition('defaultFaces', bc.ZeroGradient())
CRDisp.set_boundary_condition('bottom', bc.ZeroGradient())
CRDisp.set_boundary_condition('top', bc.FixedValue(0))
CRDisp.set_boundary_condition('topCR', bc.FixedValue(0))
CRDisp.set_boundary_condition('bottomCR', bc.ZeroGradient())


timeFolder0.append(defaultFlux)
timeFolder0.append(T)
timeFolder0.append(U)
timeFolder0.append(p)
timeFolder0.append(p_rgh)
timeFolder0.append(Tmech)
timeFolder0.append(D)
timeFolder0.append(fuelDisp)
timeFolder0.append(CRDisp)


#==============================================================================*
# Neutronics solver

neutronicsSolver = ffn.NeutronicsSolver(
    solver="diffusionNeutronics",
    mesh=nMesh,
    region=nMesh.region,
    isMeshDeformation=True,
    displacementFieldName="disp",
    power=8e+08,
    keff=0.9388902,
)
# neutronicsSolver.timeFolder = timeFolder0

# neutronicsSolver.fvSchemes.ddtSchemes['default'] = 'steadyState'
# neutronicsSolver.fvSchemes.laplacianSchemes["div(facePhi_,angularFlux_)"] = "Gauss upwind"

# neutronicsSolver.fvSolution.solvers["\"flux.*\""] = ffn.fvSolutionSolver(
#     solver="PCG",
#     preconditioner="DIC",
#     tolerance=1e-6,
#     relTol=1e-3,
# )
# neutronicsSolver.fvSolution.solvers["\"precStar.*\""] = neutronicsSolver.fvSolution.solvers["\"prec.*\""]

# print(neutronicsSolver.fvSolution.solvers["\"flux.*\""])

# neutronicsSolver.neutronTransportOptions.implicitPredictor = True

# readZone = ffn.NuclearDataZone(f'bu{0}', fuelFraction=0.3)
# readZone.read_from_serpent("../../main_res.m")

# refState = ffn.NuclearDataState('reference')
# refState.add_zone(readZone)

# neutronicsSolver.nuclearData.add_state(refState)

# nuclearData = neutronicsSolver.nuclearData
# nuclearData.axialOrientation = ffn.Vector(1, 1, 1)

# refState = ffn.NuclearDataState(
#     name="reference",
#     zones=[
#         ffn.NuclearDataZone(
#             "zone0",
#             fuelFraction=0.4,
#             removalXS=[5.0082],
#             nuFissionXS=[5.5675],
#             powerXS=[1],
#             scatteringMatrixP0=[[0.1]],
#             discFactor=[1],
#             chiPrompt=[1],
#             chiDelayed=[1],
#             inverseVelocity=[7.1e-08],
#             diffusionCoefficient=[0.1275],
#             integralFlux=[1],
#             decayConstant=[0.0125371, 0.0300828, 0.109879, 0.325484, 1.3036, 9.51817],
#             delayedFraction=[7.2315e-05, 0.000609661, 0.000471181, 0.00118907, 0.000445487, 9.58515e-05],
#         )
#     ]
# )

# nuclearData.add_state(refState)

#==============================================================================*
# Thermal-hydraulics solver

thSolver = ffn.ThermalHydraulicsSolver(
    region="fluidRegion",
    solver="onePhase",
    removeBaffles=True,
    mesh=thMesh,
    isSetFvSolutionToDefault=False
)

thSolver.thermophysicalProperties = ffn.thermophysicalProperty.SodiumConst()
thSolver.thermophysicalProperties.rho = 860
thSolver.thermophysicalProperties.Cp = 1646.97 - 0.831583*inletTemperature + 4.31182e-4*inletTemperature**2
# thSolver.thermophysicalProperties = ffn.thermophysicalProperty.SodiumPolynomial()

thSolver.turbulenceProperties.simulationType = 'laminar'

structureZones = [
    "diagrid", "upperGasPlenum", "upperReflector", "lowerGasPlenum",
    "lowerReflector", "radialReflector", "follower", "controlRod"
]

passiveStructures = ffn.StructureProperty(
    structureZones,
    volumeFraction=0.718520968,
    Dh=0.00365
)
passiveStructures.add_passive_structure(volumetricArea=5, rhoCp=4.8e6, T=inletTemperature)
thSolver.add_structure_property(passiveStructures)


core = ffn.StructureProperty(zones=["innerCore", "outerCore"], volumeFraction=0.718520968, Dh=0.00365)
core.add_power_model(
    ffn.NuclearFuelPin(
        fuelInnerRadius=0.0012,
        fuelOuterRadius=0.004715,
        cladInnerRadius=0.004865,
        cladOuterRadius=0.005365,
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
)
core.add_passive_structure(volumetricArea=5, rhoCp=4.8e6, T=inletTemperature)
thSolver.add_structure_property(core)

thSolver.add_drag_model(
    ffn.ReynoldsPower(
        coeff=0.687, exp=-0.25,
        zones=["innerCore", "outerCore"] + structureZones
    )
)

heatTransferByRegime = ffn.HeatTransferByRegime(
    "lamTurb",
    zones=["innerCore", "outerCore"] + structureZones
)
heatTransferByRegime.add_regime(
    ffn.NusseltReynoldsPrandtlPower(
        const=4, coeff=0, expRe=0, expPr=0, zones=['laminar']
    )
)
heatTransferByRegime.add_regime(
    ffn.NusseltReynoldsPrandtlPower(
        const=4.82, coeff=0.0185, expRe=0.827, expPr=0.827, zones=['turbulent']
    )
)
thSolver.add_heat_transfer_model(heatTransferByRegime)

lamTurb = ffn.RegimeMapOneParameter(name="lamTurb", parameter="Re")
lamTurb.add_regime("laminar", 0, 1000)
lamTurb.add_regime("turbulent", 2300, 2301)
thSolver.add_regime_map_model(lamTurb)


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


fvSolution = ffn.fvSolution()

fvSolution.append("p_rgh", ffn.fvSolutionSolver(
    solver='GAMG', smoother='DIC', tolerance=1e-7, relTol=0
))
fvSolution.append("p_rghFinal", ffn.fvSolutionSolver(
    solver='GAMG', smoother='DIC', tolerance=1e-7, relTol=0
))
fvSolution.append("e", ffn.fvSolutionSolver(
    solver='smoothSolver', smoother='symGaussSeidel', tolerance=1e-7, relTol=0, minIter=0
))
fvSolution.append("h", ffn.fvSolutionSolver(
    solver='smoothSolver', smoother='symGaussSeidel', tolerance=1e-7, relTol=0, minIter=0
))
fvSolution.append('".*"', ffn.fvSolutionSolver(
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
# Thermomechanics solver

tmSolver = ffn.OffbeatSolver(
    region="thermoMechanicalRegion",
    solver="extendedThermoMechanics",
    mesh=tmMesh,
    thermalSolverOptions=ffn.SolidConductionThermalSolverOptions(),
    mechanicsSolverOptions=ffn.SmallStrainMechanicsSolverOptions(),
    materialProperties="byZone",
    # rheology="byMaterial",
    # heatSource="fromLatestTime",
    # sliceMapper="autoAxialSlices"
)

baseMat = ffn.thermomechanicalMaterial.ConstantMaterial(
    name="base",
    rho=1000,
    Cp=200,
    k=5,
    emissivity=0,
    E=1e9,
    nu=0.3,
    G=0,
    alpha=1.8e-5,
    Tref=inletTemperature
)

innerCoreMat = copy(baseMat)
innerCoreMat.name = "innerCore"
innerCoreMat.alphaFuel = 1.1e-5
innerCoreMat.TFuelRef = inletTemperature

outerCoreMat = copy(innerCoreMat)
outerCoreMat.name = "outerCore"

followerMat = copy(baseMat)
followerMat.name = "follower"

controlRodMat = copy(baseMat)
controlRodMat.name = "controlRod"
controlRodMat.alphaCR = 5.4e-5
controlRodMat.TCRRef = inletTemperature

diagridMat = copy(baseMat)
diagridMat.name = "diagrid"

radialReflectorMat = copy(baseMat)
radialReflectorMat.name = "radialReflector"

restMat = copy(baseMat)
restMat.name = "rest"

softStructureMat = copy(baseMat)
softStructureMat.name = "softStructure"

tmSolver.add_material(innerCoreMat)
tmSolver.add_material(outerCoreMat)
tmSolver.add_material(followerMat)
tmSolver.add_material(controlRodMat)
tmSolver.add_material(diagridMat)
tmSolver.add_material(radialReflectorMat)
tmSolver.add_material(restMat)
tmSolver.add_material(softStructureMat)

# tmSolver.add_relaxation_on_field('T', 0.8)


#==============================================================================*
# Solvers

solvers = ffn.Solvers([neutronicsSolver, thSolver, tmSolver])

# for solver in solvers:
#     solver.decomposeParDict.numberOfSubdomains = 8
#     solver.decomposeParDict.method = "scotch"


#==============================================================================*
# Coupling

coupling = ffn.Coupling(solvers=solvers)

coupling.add_field_transfer(neutronicsSolver, thSolver, "powerDensity", "powerDensityNeutronics")
coupling.add_field_transfer(neutronicsSolver, thSolver, "secondaryPowerDensity", "powerDensityNeutronicsToLiquid")

coupling.add_field_transfer(thSolver, neutronicsSolver, "T", "TCool")
coupling.add_field_transfer(thSolver, neutronicsSolver, "thermo:rho", "rhoCool")
coupling.add_field_transfer(thSolver, neutronicsSolver, "T.fuelAvForNeutronics", "TFuel")
coupling.add_field_transfer(thSolver, neutronicsSolver, "T.cladAvForNeutronics", "TClad")
coupling.add_field_transfer(thSolver, neutronicsSolver, "T.passiveStructure", "TStructMech")

coupling.add_field_transfer(tmSolver, neutronicsSolver, "meshDisp", "disp")

coupling.add_field_transfer(neutronicsSolver, tmSolver, "powerDensity", "Q")

coupling.add_field_transfer(thSolver, tmSolver, "T.passiveStructure", "TStructFromTH")
coupling.add_field_transfer(thSolver, tmSolver, "T.fuelAvForNeutronics", "TFuel")


#==============================================================================*
# Settings


model = ffn.Model(solvers=solvers, coupling=coupling, timeFolders=[timeFolder0])

settings: ffn.ControlDict = model.settings

settings.application = 'GeN-Foam'
settings.endTime = 200
settings.deltaT = 1
settings.writeControl = 'runTime'
settings.writeInterval = 10
settings.writePrecision = 7
settings.runTimeModifiable = True
settings.adjustTimeStep = False

# Add idxField to visualize cellZones, not used in calculations
idxField = neutronicsSolver.create_zone_field(nMesh.cellZones)
timeFolder0.append(idxField)

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

    for filename in [
        'nuclearData', 'XSaxialExpansion', 'XSradialExpansion', 'XSref',
        'XSrhoCool500kgm3', 'XSTClad1950K', 'XSTFuel1200K', 'CRmove'
    ]:
        ffn.copyFolder(f"./XS/{filename}", f"constant/{nMesh.region}")


    #==============================================================================*
    # Preprocessing

    ffn.run_preprocessing(model=model)

    model.plot_mesh(region=nMesh, show_edges=True)
    model.plot_mesh(region=thMesh, show_edges=True)
    model.plot_mesh(region=nMesh, fieldName=idxField.name, show_edges=True, cmap="tab10", limits=[0.5, 10.5])
    model.plot_slice(region=nMesh, fieldName=idxField.name, show_edges=True, cmap="tab10", normal="x", limits=[0.5, 10.5])
    model.plot_slice(region=nMesh, fieldName=idxField.name, show_edges=True, cmap="tab10", normal="y", limits=[0.5, 10.5])
    model.plot_slice(region=nMesh, fieldName=idxField.name, show_edges=True, cmap="tab10", normal="z", limits=[0.5, 10.5])
    model.plot_boundary(region=nMesh, boundaryName='top')
    model.plot_boundary(region=thMesh, boundaryName='baffle0')

    #==============================================================================*
    # Run

    ffn.run(model=model)


#==============================================================================*
# Post-processing

model.plot_residuals(
    parameters=['fluxStar0'],
    title="Neutronics"
)
model.plot_residuals(
    parameters=['p_rgh', 'h'],
    title="Thermal-hydraulics"
)
model.plot_residuals(
    parameters=['fuelDisp', 'CRDisp', 'Dx', 'Dy', 'Dz'],
    title="Thermomechanics"
)

lastTime = model.get_time_steps()[-1]

print(f"keff = {model.keff(time=lastTime)}")

res = model.get_keff_from_log()

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(res['time'], res['keff'])
ax.set_xlabel("Time [s]")
ax.set_ylabel(r"k$_\text{eff}$")
fig.tight_layout()
fig.savefig("fig_results_keff.png")

model.plot_mesh(
    region=thMesh,
    time=lastTime,
    fieldName='T',
    cmap='RdBu_r',
    unit='K'
)

model.plot_slice(
    region=nMesh,
    time=lastTime,
    fieldName='flux0',
    normal='z',
    cmap='Blues',
    unit='neutron/m2/s',
    show_edges=True
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

for field in ['TCool', 'TClad', 'TFuel', 'TStructMech']:
    model.plot_slice(
        region=nMesh,
        time=lastTime,
        fieldName=field,
        normal='x',
        cmap='RdBu_r',
        unit='K',
        show_edges=True
    )
    model.plot_slice(
        region=nMesh,
        time=lastTime,
        fieldName=field,
        normal='z',
        cmap='RdBu_r',
        unit='K',
        show_edges=True
    )

    minT, maxT = model.get_data_range(nMesh, lastTime, field, removeZeros=True)

    model.plot_animation(
        region=nMesh,
        fieldName=field,
        cmap="RdBu_r",
        normal="y",
        unit="K",
        fps=2,
        limits=[minT, maxT]
    )


model.plot_animation(
    region=nMesh,
    fieldName='TCool',
    cmap="RdBu_r",
    unit="K",
    fps=1
)


#==============================================================================*
