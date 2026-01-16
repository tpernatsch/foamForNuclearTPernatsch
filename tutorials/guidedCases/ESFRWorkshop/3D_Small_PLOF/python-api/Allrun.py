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
    lc = 0.5;
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

    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            name="innerCore",
            zmin=-coreHeight/2,
            zmax=coreHeight/2,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=coreNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['I'],
        # isMergePatches=isMergePatches
    )
    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            name="outerCore",
            zmin=-coreHeight/2,
            zmax=coreHeight/2,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=coreNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['O'],
        # isMergePatches=isMergePatches
    )
    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            name="upperGasPlenum",
            zmin=coreHeight/2,
            zmax=coreHeight/2+upperFGPHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=upperFGPNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['I', 'O'],
        # isMergePatches=isMergePatches
    )

    # """
    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            name="upperReflector",
            zmin=coreHeight/2+upperFGPHeight,
            zmax=coreHeight/2+upperFGPHeight+upperReflHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=upperReflNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['I', 'O'],
        # isMergePatches=isMergePatches
    )
    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            name="lowerGasPlenum",
            zmin=-coreHeight/2-lowFGPHeight,
            zmax=-coreHeight/2,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=lowFGPNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['I', 'O'],
        # isMergePatches=isMergePatches
    )

    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            name="lowerReflector",
            zmin=-coreHeight/2-lowFGPHeight-lowReflHeight,
            zmax=-coreHeight/2-lowFGPHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=lowReflNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['I', 'O'],
        # isMergePatches=isMergePatches
    )

    # Radial reflector
    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            name="radialReflector",
            zmin=-coreHeight/2,
            zmax=coreHeight/2,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=coreNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['R'],
        # isMergePatches=isMergePatches
    )
    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            name="radialReflector",
            zmin=coreHeight/2,
            zmax=coreHeight/2+upperFGPHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=upperFGPNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['R'],
        # isMergePatches=isMergePatches
    )
    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            name="radialReflector",
            zmin=coreHeight/2+upperFGPHeight,
            zmax=coreHeight/2+upperFGPHeight+upperReflHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=upperReflNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['R'],
        # isMergePatches=isMergePatches
    )
    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            name="radialReflector",
            zmin=-coreHeight/2-lowFGPHeight,
            zmax=-coreHeight/2,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=lowFGPNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['R'],
        # isMergePatches=isMergePatches
    )
    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            name="radialReflector",
            zmin=-coreHeight/2-lowFGPHeight-lowReflHeight,
            zmax=-coreHeight/2-lowFGPHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=lowReflNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['R'],
        # isMergePatches=isMergePatches
    )

    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            name="diagrid",
            zmin=-coreHeight/2-lowFGPHeight-lowReflHeight-diagridHeight,
            zmax=-coreHeight/2-lowFGPHeight-lowReflHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=diagridNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['I', 'O', 'R', 'C'],
        # isMergePatches=isMergePatches
    )


    # Follower
    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            name="follower",
            zmin=-coreHeight/2-lowFGPHeight-lowReflHeight,
            zmax=-coreHeight/2-lowFGPHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=lowReflNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['C'],
        # isMergePatches=isMergePatches
    )
    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            name="follower",
            zmin=-coreHeight/2-lowFGPHeight,
            zmax=-coreHeight/2,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=lowFGPNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['C'],
        # isMergePatches=isMergePatches
    )
    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            name="follower",
            zmin=-coreHeight/2,
            zmax=coreHeight/2,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=coreNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['C'],
        # isMergePatches=isMergePatches
    )

    # Control rod
    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            name="controlRod",
            zmin=coreHeight/2,
            zmax=coreHeight/2+upperFGPHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=upperFGPNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['C'],
        # isMergePatches=isMergePatches
    )
    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            name="controlRod",
            zmin=coreHeight/2+upperFGPHeight,
            zmax=coreHeight/2+upperFGPHeight+upperReflHeight,
            pitch=flatToFlat,
            x=x, y=y,
            nr=nr,
            nt=nt,
            nz=upperReflNodes,
            isAddAllBC=True
        ),
        lattice=lattice,
        latticeType="hexagon",
        nx=latticeNY, ny=latticeNY, pitch=pitch,
        elementsToPlace=['C'],
        # isMergePatches=isMergePatches
    )
    # """

    blockMesh.add_merge_patch_pairs(
        includeFacename=['Top', 'Bottom'],
        excludeFacename=['Wall']
    )
    if (isFluidMesh):
        bafflesFaces = blockMesh.add_baffles(
            includeFacename=['Wall'],
            excludeFacename=['Top', 'Bottom']
        )
    else:
        blockMesh.add_merge_patch_pairs(
            includeFacename=['Wall'],
            excludeFacename=['Top', 'Bottom'],
        )


    externalWalls = blockMesh.get_standalone_faces(includeFacename=['Wall'])

    blockMesh.merge_patches_with_name(name="wall", regex='|'.join([face.name for face in externalWalls]), patchType="wall")
    blockMesh.merge_patches_with_name(name="bottom", regex="diagridBottom_.*")
    blockMesh.merge_patches_with_name(name="top", regex="upperReflectorTop_.*")
    blockMesh.merge_patches_with_name(name="top", regex="controlRodTop_.*")
    blockMesh.merge_patches_with_name(name="top", regex="radialReflectorTop_.*")

    return(blockMesh)


# nMesh = createMesh(region="neutroRegion", isFluidMesh=False)
# thMesh = createMesh(region="fluidRegion", isFluidMesh=True) #, gap=0.009)
nMesh = mesh.PolyMesh(
    region="neutroRegion",
    srcpath="./meshes/polyMeshNeutro"
)
thMesh = mesh.PolyMesh(region="fluidRegion", srcpath="./meshes/polyMeshFluid")
tmMesh = mesh.PolyMesh(region="thermoMechanicalRegion", srcpath="./meshes/polyMeshMeca")
# nMesh = mesh.PolyMesh(region="neutroRegion", srcpath="../steadyState_THNeutronicsAndTM/constant/neutroRegion/polyMesh")
# thMesh = mesh.PolyMesh(region="fluidRegion", srcpath="../steadyState_THNeutronicsAndTM/constant/fluidRegion/polyMesh")
# tmMesh = mesh.PolyMesh(region="thermoMechanicalRegion", srcpath="../steadyState_THNeutronicsAndTM/constant/thermoMechanicalRegion/polyMesh")



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
U.set_boundary_condition('bottom', bc.UniformFixedValue(
    uniformValue=ffn.Table([
        (200, U.internalField),
        (300, U.internalField),
        (310, ffn.Vector(0, 0, 2.279557)),
        (320, ffn.Vector(0, 0, 1.1397785)),
        (3000, ffn.Vector(0, 0, 1.1397785)),
    ])
))



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
neutronicsSolver.neutronTransportOptions.maxNeutronIterations = 50

neutronicsSolver.controlRodMove.add_control_rod_movement(
    cellZoneName="controlRod",
    startTime=310,
    endTime=310.25,
    speed=9,
    initialDistanceFromMeshCR=0,
    followerName="follower"
)

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
# thSolver.thermophysicalProperties.rho = 860
# thSolver.thermophysicalProperties.Cp = 1646.97 - 0.831583*inletTemperature + 4.31182e-4*inletTemperature**2
thSolver.thermophysicalProperties = ffn.thermophysicalProperty.SodiumPolynomial()
thSolver.thermophysicalProperties.Cp = ffn.Polynome(1646.97, -0.831587, 4.31182e-04)
thSolver.thermophysicalProperties.rho = ffn.Polynome(860)
thSolver.thermophysicalProperties.mu = ffn.Polynome(2e-4)
thSolver.thermophysicalProperties.kappa = ffn.Polynome(60)

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

core = ffn.StructureProperty(
    zones=["innerCore", "outerCore"],
    volumeFraction=0.718520968,
    Dh=0.00365
)
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
core.add_passive_structure(
    volumetricArea=5, rhoCp=4.8e6, T=inletTemperature
)
thSolver.add_structure_property(core)

thSolver.add_drag_model(
    ffn.ReynoldsPower(
        coeff=0.687, exp=-0.25,
        zones=["innerCore", "outerCore"] + structureZones
    )
)

heatTransferByRegime = ffn.ByRegime(
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

thSolver.pimpleOptions = ffn.PimpleOptions(
    nCorrectors=2,
    nNonOrthogonalCorrectors=0,
    nOuterCorrectors=6,
    solveEnergy=True,
    solveFluidMechanics=True,
    momentumMode="faceCentered",
    correctUntilConvergence=None,
    porousInterfaceSharpness=None,
    minMagU=None,
    minNOuterCorrectors=None,
)

thSolver.add_relaxation_on_equation('h', 0.7)
thSolver.add_relaxation_on_field('T', 0.7)


#==============================================================================*
# Thermomechanics solver

tmSolver = ffn.OffbeatSolver(
    region="thermoMechanicalRegion",
    solver="extendedThermoMechanics",
    mesh=tmMesh,
    thermalSolver="solidConduction",
    mechanicsSolver="smallStrain",
    materialProperties="byZone",
    rheology="byMaterial",
    heatSource="fromLatestTime",
    sliceMapper="autoAxialSlices"
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

tmSolver.add_relaxation_on_field('D', 0.7)


#==============================================================================*
# Solvers

solvers = ffn.Solvers([thSolver, neutronicsSolver, tmSolver])

for solver in solvers:
    solver.decomposeParDict.numberOfSubdomains = 4
    solver.decomposeParDict.method = 'simple'
    solver.decomposeParDict.simpleCoeffs['n'] = ffn.Vector(1, 1, 4)

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

coupling.plot_coupling_graph()
coupling.plot_solving_graph()


#==============================================================================*
# Settings


model = ffn.Model(
    solvers=solvers,
    coupling=coupling,
    timeFolders=[timeFolder0]
)

settings: ffn.ControlDict = model.settings

settings.application = 'GeN-Foam'
settings.endTime = 200
settings.deltaT = 1
settings.writeControl = 'runTime'
settings.writeInterval = 10
settings.writePrecision = 7
settings.runTimeModifiable = True
settings.adjustTimeStep = False


#==============================================================================*
# Main

if __name__ == "__main__":
    # Export to OpenFOAM
    ffn.allclean()
    model.export_to_openfoam()

    print(model)

    #==============================================================================*
    # Duplicate and overwrite files

    for filename in [
        'nuclearData', 'XSaxialExpansion', 'XSradialExpansion', 'XSref',
        'XSrhoCool500kgm3', 'XSTClad1950K', 'XSTFuel1200K'
    ]:
        ffn.copyFolder(f"../steadyState_THNeutronicsAndTM/constant/neutroRegion/{filename}", f"constant/{nMesh.region}")


    #==============================================================================*
    # Preprocessing

    ffn.run_preprocessing(model=model)

    model.plot_mesh(region=nMesh, show_edges=True)
    model.plot_mesh(region=thMesh, show_edges=True)
    model.plot_boundary(region=nMesh, boundaryName='top')
    model.plot_boundary(region=thMesh, boundaryName='baffle0')

    #==============================================================================*
    # Run

    ffn.run(model=model)


#==============================================================================*
