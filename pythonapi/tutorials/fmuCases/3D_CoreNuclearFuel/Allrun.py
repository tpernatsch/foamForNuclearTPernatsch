"""

"""
#==============================================================================*
# Imports

from copy import copy, deepcopy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh
import foamForNuclear.thermophysicalProperty as thermo



def latticeReplace(lattice, elementsToInclude: list[str]=None):
    flatLattice = ' '.join([line for line in lattice.split('\n') if line != ""])
    elements = {e for e in flatLattice.split()}

    newLattice = lattice
    for element in elementsToInclude:
        newLattice = newLattice.replace(element, 'F')

    for element in elements:
        if (element not in elementsToInclude):
            newLattice = newLattice.replace(element, 'e')

    return(newLattice)


#==============================================================================*
# Mesh

pitchFuelAssembly = 0.167

activeCoreHeight = 1.02 # m
upperFuelPlenumHeight = 0.945 # m
lowerFuelPlenumHeight = 0.925 # m
diagridHeight = 0.610 # m

activeCoreSegments = 15
upperFuelPlenumSegments = 5
lowerFuelPlenumSegments = 5
diagridSegments = 5

lattice = """
ce ce ce ce ce ce ce ce ce ce ce ce ce ce ce ce ce ce ce ce ce
 ce ce ce ce ce ce ce ce ce ce ce ce RA RA RA RA RA RA ce ce ce
  ce ce ce ce ce ce ce ce ce ce RA RA RA RA RA RA RA RA RA ce ce
   ce ce ce ce ce ce ce ce RA RA RA RA OA OA OA OA RA RA RA RA ce
    ce ce ce ce ce ce ce RA RA RA OA OA OA OA OA OA OA RA RA RA ce
     ce ce ce ce ce ce RA RA OA OA OA OA OA OA OA OA OA OA RA RA ce
      ce ce ce ce ce RA RA OA OA OA IA IA IA IA IA OA OA OA RA RA ce
       ce ce ce ce RA RA OA OA OA IA IA IA IA IA IA OA OA OA RA RA ce
        ce ce ce RA RA OA OA OA IA IA IA IA IA IA IA OA OA OA RA RA ce
         ce ce ce RA RA OA OA IA IA IA EA EA IA IA IA OA OA RA RA ce ce
          ce ce RA RA OA OA IA IA IA EA EA EA IA IA IA OA OA RA RA ce ce
           ce ce RA RA OA OA IA IA IA EA EA IA IA IA OA OA RA RA ce ce ce
            ce RA RA OA OA OA IA IA IA IA IA IA IA OA OA OA RA RA ce ce ce
             ce RA RA OA OA OA IA IA IA IA IA IA OA OA OA RA RA ce ce ce ce
              ce RA RA OA OA OA IA IA IA IA IA OA OA OA RA RA ce ce ce ce ce
               ce RA RA OA OA OA OA OA OA OA OA OA OA RA RA ce ce ce ce ce ce
                ce RA RA RA OA OA OA OA OA OA OA RA RA RA ce ce ce ce ce ce ce
                 ce RA RA RA RA OA OA OA OA RA RA RA RA ce ce ce ce ce ce ce ce
                  ce ce RA RA RA RA RA RA RA RA RA ce ce ce ce ce ce ce ce ce ce
                   ce ce ce RA RA RA RA RA RA ce ce ce ce ce ce ce ce ce ce ce ce
                    ce ce ce ce ce ce ce ce ce ce ce ce ce ce ce ce ce ce ce ce ce
"""


innerCoreLattice = latticeReplace(lattice, ['IA'])
outerCoreLattice = latticeReplace(lattice, ['OA'])
nXY = len([line for line in innerCoreLattice.split('\n') if line != ''])

nPinsInnerCore = len([element for element in ' '.join(innerCoreLattice.split("\n")).split() if element == 'F'])
nPinsOuterCore = len([element for element in ' '.join(outerCoreLattice.split("\n")).split() if element == 'F'])


def generateNeutronicMesh(region: str):
    blockMesh = mesh.BlockMesh(region=region)

    # Upper section
    for zoneName, elements in [
        ('target', ['EA']),
        ('reflector', ['RA']),
        ('upperFuelPlenum', ['IA', 'OA'])
    ]:
        blockMesh.lattice_placement(
            funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
                zoneName,
                zmin=activeCoreHeight/2,
                zmax=activeCoreHeight/2+upperFuelPlenumHeight,
                pitch=pitchFuelAssembly,
                x=x, y=y,
                nr=1, nt=1, nz=upperFuelPlenumSegments,
                isAddBoundaryConditions=True
            ),
            lattice=lattice,
            nx=nXY, ny=nXY,
            pitch=pitchFuelAssembly,
            latticeType='hexagon',
            elementsToPlace=elements
        )

    # Core section
    for zoneName, elements in [
        ('innerCore', ['IA']),
        ('outerCore', ['OA']),
        ('target',    ['EA']),
        ('reflector', ['RA'])
    ]:
        blockMesh.lattice_placement(
            funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
                zoneName,
                zmin=-activeCoreHeight/2,
                zmax=activeCoreHeight/2,
                pitch=pitchFuelAssembly,
                x=x, y=y,
                nr=1, nt=1, nz=activeCoreSegments,
                isAddBoundaryConditions=True
            ),
            lattice=lattice,
            nx=nXY, ny=nXY,
            pitch=pitchFuelAssembly,
            latticeType='hexagon',
            elementsToPlace=elements
        )

    # Lower section
    for zoneName, elements in [
        ('target', ['EA']),
        ('reflector', ['RA']),
        ('lowerFuelPlenum', ['IA', 'OA'])
    ]:
        blockMesh.lattice_placement(
            funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
                zoneName,
                zmin=-activeCoreHeight/2-lowerFuelPlenumHeight,
                zmax=-activeCoreHeight/2,
                pitch=pitchFuelAssembly,
                x=x, y=y,
                nr=1, nt=1, nz=lowerFuelPlenumSegments,
                isAddBoundaryConditions=True
            ),
            lattice=lattice,
            nx=nXY, ny=nXY,
            pitch=pitchFuelAssembly,
            latticeType='hexagon',
            elementsToPlace=elements
        )

    # Diagrid
    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            'diagrid',
            zmin=-activeCoreHeight/2-lowerFuelPlenumHeight-diagridHeight,
            zmax=-activeCoreHeight/2-lowerFuelPlenumHeight,
            pitch=pitchFuelAssembly,
            x=x, y=y,
            nr=1, nt=1, nz=diagridSegments,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        nx=nXY, ny=nXY,
        pitch=pitchFuelAssembly,
        latticeType='hexagon',
        elementsToPlace=['IA', 'OA', 'EA', 'RA']
    )

    # Merge patches
    blockMesh.add_merge_patch_pairs(
        includeFacename=['Top', 'Bottom'],
        excludeFacename=['Wall']
    )
    blockMesh.add_merge_patch_pairs(
        includeFacename=['Wall'],
        excludeFacename=['Top', 'Bottom'],
    )

    externalWalls = blockMesh.get_standalone_faces(
        includeFacename=['Wall'],
        excludeFacename=['Top', 'Bottom']
    )

    # Merge into unique patches
    blockMesh.merge_patches_with_name(
        name="wall",
        includeFacename=[face.name for face in externalWalls],
        patchType="wall"
    )
    blockMesh.merge_patches_with_name(
        name="bottom",
        includeFacename=["diagridBottom_"]
    )
    blockMesh.merge_patches_with_name(
        name="top",
        includeFacename=["reflectorTop_", "targetTop_", "upperFuelPlenumTop_"]
    )

    return(blockMesh)


def generateThermalhydraulicMesh(region: str):
    blockMesh = mesh.BlockMesh(region=region)

    # # Upper section
    for zoneName, elements in [
        ('target', ['EA']),
        ('reflector', ['RA']),
        ('upperPlenum', ['IA', 'OA'])
    ]:
        blockMesh.lattice_placement(
            funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
                zoneName,
                zmin=activeCoreHeight/2+0.12,
                zmax=activeCoreHeight/2+upperFuelPlenumHeight,
                pitch=pitchFuelAssembly,
                x=x, y=y,
                nr=1, nt=1, nz=upperFuelPlenumSegments,
                isAddBoundaryConditions=True
            ),
            lattice=lattice,
            nx=nXY, ny=nXY,
            pitch=pitchFuelAssembly,
            latticeType='hexagon',
            elementsToPlace=elements
        )

    # Core section
    for zoneName, elements in [
        ('innerCore', ['IA']),
        ('outerCore', ['OA']),
        ('target',    ['EA']),
        ('reflector', ['RA'])
    ]:
        blockMesh.lattice_placement(
            funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
                zoneName,
                zmin=-activeCoreHeight/2,
                zmax=activeCoreHeight/2,
                pitch=pitchFuelAssembly,
                x=x, y=y,
                nr=1, nt=1, nz=activeCoreSegments,
                isAddBoundaryConditions=True
            ),
            lattice=lattice,
            nx=nXY, ny=nXY,
            pitch=pitchFuelAssembly,
            latticeType='hexagon',
            elementsToPlace=elements
        )
        blockMesh.lattice_placement(
            funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
                zoneName,
                zmin=-activeCoreHeight/2-0.55,
                zmax=-activeCoreHeight/2,
                pitch=pitchFuelAssembly,
                x=x, y=y,
                nr=1, nt=1, nz=5,
                isAddBoundaryConditions=True
            ),
            lattice=lattice,
            nx=nXY, ny=nXY,
            pitch=pitchFuelAssembly,
            latticeType='hexagon',
            elementsToPlace=elements
        )
        blockMesh.lattice_placement(
            funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
                zoneName,
                zmin=activeCoreHeight/2,
                zmax=activeCoreHeight/2+0.12,
                pitch=pitchFuelAssembly,
                x=x, y=y,
                nr=1, nt=1, nz=2,
                isAddBoundaryConditions=True
            ),
            lattice=lattice,
            nx=nXY, ny=nXY,
            pitch=pitchFuelAssembly,
            latticeType='hexagon',
            elementsToPlace=elements
        )

    # Lower section
    for zoneName, elements in [
        ('target', ['EA']),
        ('reflector', ['RA']),
        ('lowerPlenum', ['IA', 'OA'])
    ]:
        blockMesh.lattice_placement(
            funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
                zoneName,
                zmin=-activeCoreHeight/2-lowerFuelPlenumHeight,
                zmax=-activeCoreHeight/2-0.55,
                pitch=pitchFuelAssembly,
                x=x, y=y,
                nr=1, nt=1, nz=lowerFuelPlenumSegments,
                isAddBoundaryConditions=True
            ),
            lattice=lattice,
            nx=nXY, ny=nXY,
            pitch=pitchFuelAssembly,
            latticeType='hexagon',
            elementsToPlace=elements
        )

    # Diagrid
    blockMesh.lattice_placement(
        funcElementGenerator=lambda name, x, y: blockMesh.create_hexagon_prism_along_z(
            'diagrid',
            zmin=-activeCoreHeight/2-lowerFuelPlenumHeight-diagridHeight,
            zmax=-activeCoreHeight/2-lowerFuelPlenumHeight,
            pitch=pitchFuelAssembly,
            x=x, y=y,
            nr=1, nt=1, nz=diagridSegments,
            isAddBoundaryConditions=True
        ),
        lattice=lattice,
        nx=nXY, ny=nXY,
        pitch=pitchFuelAssembly,
        latticeType='hexagon',
        elementsToPlace=['IA', 'OA', 'EA', 'RA']
    )

    # Merge patches first
    blockMesh.add_merge_patch_pairs(
        includeFacename=['Top', 'Bottom'],
        excludeFacename=['Wall']
    )
    # Create baffles
    bafflesFaces = blockMesh.add_baffles(
        includeFacename=['Wall'],
        excludeFacename=['Top', 'Bottom']
    )
    # Extract lateral external faces
    externalWalls = blockMesh.get_standalone_faces(
        includeFacename=['Wall'],
        excludeFacename=['Top', 'Bottom']
    )

    # Merge into unique patches
    blockMesh.merge_patches_with_name(
        name="wall",
        includeFacename=[face.name for face in externalWalls],
        patchType="wall"
    )
    blockMesh.merge_patches_with_name(
        name="inlet",
        includeFacename=["diagridBottom_"]
    )
    blockMesh.merge_patches_with_name(
        name="outlet",
        includeFacename=["reflectorTop_", "targetTop_", "upperPlenumTop_"]
    )

    return(blockMesh)


if False:
    nMesh = generateNeutronicMesh(region="neutroRegion")
else:
    nMesh = mesh.PolyMesh(region='neutroRegion', srcpath="../storedData/neutroPolyMesh")

if False:
    thMesh = generateThermalhydraulicMesh(region="fluidRegion")
else:
    thMesh = mesh.PolyMesh(region='fluidRegion', srcpath="../storedData/fluidPolyMesh")

#==============================================================================*
# Time folder

timeFolder0 = ffn.TimeFolder(0)

defaultFlux = ffn.Field("defaultFlux", region=nMesh.region)
defaultFlux.dimensions = ffn.Dimension(default='neutronFlux')
defaultFlux.internalField = 1
defaultFlux.set_boundary_condition('wall', bc.FixedValue(0))
defaultFlux.set_boundary_condition('bottom', bc.FixedValue(0))
defaultFlux.set_boundary_condition('top', bc.FixedValue(0))

T = ffn.Field("T", region=thMesh.region)
T.dimensions = ffn.Dimension(default='T')
T.internalField = 673.15
T.set_boundary_condition('inlet', bc.FixedValue(T.internalField))
T.set_boundary_condition('outlet', bc.ZeroGradient())
T.set_boundary_condition('wall', bc.ZeroGradient())
T.set_boundary_condition('baffle0', bc.ZeroGradient())
T.set_boundary_condition('baffle1', bc.ZeroGradient())

U = ffn.Field("U", region=thMesh.region)
U.dimensions = ffn.Dimension(default='U')
U.internalField = ffn.Vector(0, 0, 0.767)
U.set_boundary_condition('inlet', bc.FixedValue(U.internalField))
U.set_boundary_condition('outlet', bc.ZeroGradient())
U.set_boundary_condition('wall', bc.Slip())
U.set_boundary_condition('baffle0', bc.Slip())
U.set_boundary_condition('baffle1', bc.Slip())

p = ffn.Field("p", region=thMesh.region)
p.dimensions = ffn.Dimension(default='p')
p.internalField = 1e5
p.set_boundary_condition('inlet', bc.Calculated(p.internalField))
p.set_boundary_condition('outlet', bc.Calculated(p.internalField))
p.set_boundary_condition('wall', bc.Calculated(p.internalField))
p.set_boundary_condition('baffle0', bc.ZeroGradient())
p.set_boundary_condition('baffle1', bc.ZeroGradient())

p_rgh = ffn.Field("p_rgh", region=thMesh.region)
p_rgh.dimensions = ffn.Dimension(default='p')
p_rgh.internalField = 1e5
p_rgh.set_boundary_condition('inlet', bc.FixedFluxPressure(p_rgh.internalField))
p_rgh.set_boundary_condition('outlet', bc.FixedValue(p_rgh.internalField))
p_rgh.set_boundary_condition('wall', bc.FixedFluxPressure(p_rgh.internalField))
p_rgh.set_boundary_condition('baffle0', bc.ZeroGradient())
p_rgh.set_boundary_condition('baffle1', bc.ZeroGradient())

timeFolder0.append(defaultFlux)
timeFolder0.append(T)
timeFolder0.append(U)
timeFolder0.append(p)
timeFolder0.append(p_rgh)


#==============================================================================*
# Thermal-hydraulics solver

thSolver = ffn.ThermalHydraulicsSolver(
    "fluidRegion", "onePhase",
    mesh=thMesh,
    removeBaffles=True,
    isSetFvSolutionToDefault=False
)


thSolver.thermophysicalProperties = thermo.LeadBoussinesq()

thSolver.turbulenceProperties.simulationType = 'laminar'

rCladOuter = 0.00525  # Outer cladding radius
pitch = 0.01386     # Fuel pin pitch

innerCore = ffn.StructureProperty(
    zones=['innerCore'],
    pitch=pitch,
    elementDiameter=2*rCladOuter,
    latticeType='hexagon'
)

innerCore.powerModel = ffn.NuclearFuelFMU(
    volumetricArea=2.0*innerCore.volumeFraction / rCladOuter,
    xyPosLattice=ffn.XYPosLattice(
        'hexagonal',
        nElements=[nXY, nXY],
        pitch=pitchFuelAssembly,
        lattice=innerCoreLattice
    ),
    axialPowerInterpolationMethod="linear",
    axialLocationsNameToFMU="gfAxialLocations",
    axialLocations=[
        -1.06000, -0.93929, -0.81857, -0.69786, -0.57714, -0.45643, -0.33571,
        -0.21500, -0.09429, 0.02643, 0.14714, 0.26786, 0.38857, 0.50929, 0.63000
    ],
    fuelLength=activeCoreHeight,
    fuelFraction=3.635160e-01,
    radialBasisFunctionMethod='polyharmonicSpline',
    initPowerDensity=0,
    initHeatFlux=0,
    initTstruct=T.internalField,
    initRhoCpdTdt=0,
    avgPowerDensityNameToFMU=[f'gfAvgPowerDensity{i}' for i in range(nPinsInnerCore)],
    axialProfilePowerDensityNameToFMU=[f'gfAxialProfilePowerDensity{i}' for i in range(nPinsInnerCore)],
    TfluidNameToFMU=[f'gfTfluid{i}' for i in range(nPinsInnerCore)],
    htcNameToFMU=[f'gfHtc{i}' for i in range(nPinsInnerCore)],
    heatFluxNameFromFMU=[f'gfSurfacePower{i}' for i in range(nPinsInnerCore)],
    rhoCpdTdtNameFromFMU=[f'gfEnthalpyPower{i}' for i in range(nPinsInnerCore)],
    TFuelNameFromFMU=[f'gfTFuel{i}' for i in range(nPinsInnerCore)],
    TCladNameFromFMU=[f'gfTClad{i}' for i in range(nPinsInnerCore)],
    relaxationFactorTfluid=0.8,
    relaxationFactorTstruct=0.8,
    relaxationFactorPowerDensity=0.8,
    relaxationFactorHtc=0.8
)

outerCore = deepcopy(innerCore)
outerCore.zones = ["outerCore"]
outerCore.powerModel.xyPosLattice.lattice = outerCoreLattice
outerCore.powerModel.avgPowerDensityNameToFMU = [f'gfAvgPowerDensity{i}' for i in range(nPinsInnerCore, nPinsInnerCore+nPinsOuterCore)]
outerCore.powerModel.axialProfilePowerDensityNameToFMU = [f'gfAxialProfilePowerDensity{i}' for i in range(nPinsInnerCore, nPinsInnerCore+nPinsOuterCore)]
outerCore.powerModel.TfluidNameToFMU = [f'gfTfluid{i}' for i in range(nPinsInnerCore, nPinsInnerCore+nPinsOuterCore)]
outerCore.powerModel.htcNameToFMU = [f'gfHtc{i}' for i in range(nPinsInnerCore, nPinsInnerCore+nPinsOuterCore)]
outerCore.powerModel.heatFluxNameFromFMU = [f'gfSurfacePower{i}' for i in range(nPinsInnerCore, nPinsInnerCore+nPinsOuterCore)]
outerCore.powerModel.rhoCpdTdtNameFromFMU = [f'gfEnthalpyPower{i}' for i in range(nPinsInnerCore, nPinsInnerCore+nPinsOuterCore)]
outerCore.powerModel.TFuelNameFromFMU = [f'gfTFuel{i}' for i in range(nPinsInnerCore, nPinsInnerCore+nPinsOuterCore)]
outerCore.powerModel.TCladNameFromFMU = [f'gfTClad{i}' for i in range(nPinsInnerCore, nPinsInnerCore+nPinsOuterCore)]

innerCore.powerModel.xyPosLattice.plotLattice()
outerCore.powerModel.xyPosLattice.plotLattice()

fig, ax = plt.subplots()
innerCore.powerModel.xyPosLattice.plotLattice(ax=ax)
outerCore.powerModel.xyPosLattice.plotLattice(ax=ax)
ax.legend(
    bbox_to_anchor=(0., 1.02, 1., .102),
    loc='lower left',
    ncols=2,
    mode="expand",
    borderaxespad=0.
)
fig.savefig("fig_lattice_placement_core.png")
plt.close()

thSolver.add_structure_property(innerCore)
thSolver.add_structure_property(outerCore)

structureNames = ["reflector", "lowerPlenum", "diagrid", "upperPlenum", "target"]
structures = ffn.StructureProperty(
    zones=structureNames,
    volumeFraction=innerCore.volumeFraction,
    Dh=innerCore.Dh
)
structures.add_passive_structure(volumetricArea=0.91993, rhoCp=7582180, T=T.internalField)
thSolver.add_structure_property(structures)


# Regime map
regimeMap = ffn.RegimeMapOneParameter("lamTurb", "Re")
regimeMap.add_regime("laminar", 0, 1000)
regimeMap.add_regime("turbulent", 2300, 2301)
thSolver.add_regime_map_model(regimeMap)


# Drag
thSolver.add_drag_model(
    ffn.ReynoldsPower(
        coeff=0.687,
        exp=-0.25,
        zones=['innerCore', 'outerCore'] + structureNames
    )
)


# Heat transfer
heatTransfer = ffn.HeatTransferByRegime(
    regimeMap="lamTurb",
    zones=['innerCore', 'outerCore'] + structureNames
)
heatTransfer.add_regime(
    ffn.NusseltReynoldsPrandtlPower(4, 0, 0, 0, zones=['laminar'])
)
heatTransfer.add_regime(
    ffn.NusseltReynoldsPrandtlPower(7.923508, 0.031694, 0.770000, 0.770000, zones=['turbulent'])
)
thSolver.add_heat_transfer_model(heatTransfer)


#  fvSolution
thSolution = ffn.fvSolution()

thSolution.append('"p_rgh.*"', ffn.fvSolutionSolver(
    solver='GAMG', smoother='DIC', tolerance=1e-5, relTol=0
))
thSolution.append('"e.*"', ffn.fvSolutionSolver(
    solver='smoothSolver',
    smoother='symGaussSeidel',
    tolerance=1e-5, relTol=0, minIter=0
))
thSolution.append('"h.*"', ffn.fvSolutionSolver(
    solver='smoothSolver',
    smoother='symGaussSeidel',
    tolerance=1e-5, relTol=0, minIter=0
))
thSolution.append('".*"', ffn.fvSolutionSolver(
    solver='PBiCGStab',
    preconditioner='diagonal',
    tolerance=1e-5, relTol=0.001
))

thSolver.fvSolution = thSolution

thSolver.pimpleOptions.nCorrectors = 2
thSolver.pimpleOptions.nOuterCorrectors = 1
thSolver.pimpleOptions.nNonOrthogonalCorrectors = 0
thSolver.pimpleOptions.momentumMode = 'faceCentered'


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

thSolver.fvSchemes.laplacianSchemes['default'] = 'Gauss linear uncorrected'

thSolver.fvSchemes.snGradSchemes['default'] = 'uncorrected'


#==============================================================================*
# Neutronics solver

neutronicsSolver = ffn.NeutronicsSolver(
    "neutroRegion",
    "diffusionNeutronics",
    mesh=nMesh,
    power=300e6,
    eigenvalueNeutronics=False,
    externalSourceNeutronics=True
)

# Beam trip shape
tBT = 10 # s
tRamp = 1e-6

neutronicsSolver.externalSource.isExternalSource = True
neutronicsSolver.externalSource.externalSourceMode = "transient"
neutronicsSolver.externalSource.beamEnergy = 1.28160e-10
neutronicsSolver.externalSource.nuSource = 17.3
neutronicsSolver.externalSource.externalSourceModulationTimeProfile = ffn.TimeProfile(
    type='table',
    startTime=100,
    table=[
        (0,             1),
        (tRamp,         0),
        (tRamp+tBT,     0),
        (2*tRamp+tBT,   1)
    ]
)

neutronicsSolver.fvSchemes.divSchemes['div(facePhi_,angularFlux_)'] = "Gauss upwind"

neutronicsSolver.neutronTransportOptions.maxNeutronIterations = 30

#==============================================================================*
# Solvers

solvers = ffn.Solvers([neutronicsSolver, thSolver])

#==============================================================================*
# Coupling

coupling = ffn.Coupling(solvers)

coupling.add_field_transfer(neutronicsSolver, thSolver, 'powerDensity', 'powerDensityStructure')
coupling.add_field_transfer(neutronicsSolver, thSolver, 'secondaryPowerDensity', 'powerDensityLiquid')
coupling.add_field_transfer(thSolver, neutronicsSolver, 'T', 'TCool')
coupling.add_field_transfer(thSolver, neutronicsSolver, 'thermo:rho', 'rhoCool')
coupling.add_field_transfer(thSolver, neutronicsSolver, 'T.fuelAvForNeutronics', 'TFuel')
coupling.add_field_transfer(thSolver, neutronicsSolver, 'T.cladAvForNeutronics', 'TClad')
coupling.add_field_transfer(thSolver, neutronicsSolver, 'T.passiveStructure', 'TStructMech')



#==============================================================================*
# Model

endTime = 200

model = ffn.Model(
    solvers=solvers,
    coupling=coupling,
    timeFolders=[timeFolder0],
    caseFolder='caseGenfoam'
)

settings = model.settings

settings.application = "GeN-Foam"
settings.startTime = 0
settings.endTime = endTime
settings.deltaT = 0.5
settings.startFrom = 'latestTime'
settings.writeControl = 'adjustableRunTime'
settings.writeInterval = 2
settings.adjustTimeStep = True
settings.runTimeModifiable = True
settings.maxDeltaT = 10
settings.maxCo = 20
settings.maxPowerVariation = 0.025
settings.solveFMI = True

coupling.plot_coupling_graph()
coupling.plot_solving_graph()
coupling.plot_solving_flowchart()

print(model)

model.export_to_openfoam()

# Copy files only when the case is exported
for i in range(4):
    ffn.copyFolder(
        src=f'./storedData/externalSourceFlux{i}',
        dst=f"{model.caseFolder}/{settings.startTime}/{nMesh.region}"
    )
ffn.copyFolder('./storedData/nuclearData', f"{model.caseFolder}/constant/{nMesh.region}")

ffn.run_preprocessing(model=model, verbose=True)

model.plot_mesh(region=thMesh, show_edges=True)


# Export GeN-Foam as an FMU
ffn.generateCaseAsFMU(model.caseFolder, "CoreFMU")
ffn.generateCaseAsFMU("caseOffbeat", "FuelPinFMU")



connections = []
fmuNames = {"core": "CoreFMU.fmu"}
initialParametersBeforeFMUInit = {
    "core": {"port": 8001, 'outputPath': "CoreFMU"}
}

for i in range(nPinsInnerCore+nPinsOuterCore):
    fmuNames[f"fuelPin{i}"] = "FuelPinFMU.fmu"

    initialParametersBeforeFMUInit[f"fuelPin{i}"] = {
        'port': 8002+i,
        'outputPath': f'FuelPinFMU{i}',
        'timeStep': model.settings.deltaT,
        'writeInterval': model.settings.writeInterval if i in [39, 128] else model.settings.endTime
    }

    # GeN-Foam to OFFBEAT
    connections.append(
        ("core", "gfAxialLocations", f"fuelPin{i}", "offbeatAxialLocations")
    )
    connections.append(
        ("core", f"gfAvgPowerDensity{i}", f"fuelPin{i}", "offbeatQavg")
    )
    connections.append(
        ("core", f"gfAxialProfilePowerDensity{i}", f"fuelPin{i}", "offbeatHeatSrcAxial")
    )
    # connections.append(
    #     ("core", f"gfTstruct{i}", f"fuelPin{i}", "offbeatTstructAxial")
    # )
    connections.append(
        ("core", f"gfTfluid{i}", f"fuelPin{i}", "offbeatTfluidAxial")
    )
    connections.append(
        ("core", f"gfHtc{i}", f"fuelPin{i}", "offbeatHtcAxial")
    )
    # OFFBEAT to GeN-Foam
    connections.append(
        (f"fuelPin{i}", "offbeatTFuelAxial", "core", f"gfTFuel{i}")
    )
    connections.append(
        (f"fuelPin{i}", "offbeatTCladAxial", "core", f"gfTClad{i}")
    )
    connections.append(
        (f"fuelPin{i}", "offbeatSurfacePowerAxial", "core", f"gfSurfacePower{i}")
    )
    connections.append(
        (f"fuelPin{i}", "offbeatRhoCpdTdtAxial", "core", f"gfEnthalpyPower{i}")
    )

MasterRunner = ffn.FMPyMasterRunnerParallel(
# MasterRunner = ffn.FMPyMasterRunner(
    fmuNames=fmuNames,
    connections=connections,
    startTime=0,
    endTime=endTime,
    stepSize=settings.deltaT,
    minStepSize=0.05,
    maxStepSize=5,
    minNImplicitSteps=5,
    maxNImplicitSteps=101,
    adjustStepSize=False,
    # writeInterval=1,
    parametersToRecord={
        "core": [f"gfAvgPowerDensity0"]
    },
    initialParametersBeforeFMUInit=initialParametersBeforeFMUInit
)
MasterRunner.simulate()



#==============================================================================*
# Main


# if __name__ == "__main__":

#     # Steady-state calculations
#     model, neutronicsSolver = generateSteadyStateModel()

#     ffn.run(model=model, is_preprocessing=True)

#     # Transient
#     model = generateTransientModel(model=model, neutronicsSolver=neutronicsSolver)

#     # Export GeN-Foam as an FMU
#     ffn.generateCaseAsFMU(model.caseFolder, "CoreFMU")

#     runCoupledFMU()



#     # Post processing
#     model.caseFolder = 'CoreFMU1'

#     model.plot_residuals(
#         parameters=['fluxStar0'],
#         title="Neutronics"
#     )
#     model.plot_residuals(
#         parameters=['p_rgh', 'h'],
#         title="Thermal-hydraulics"
#     )

#     dataCsvCtrl = pd.read_csv("ExternalReactivityController_0.csv")

#     res = model.get_parameters_from_point_kinetics()


#     fig, ax = plt.subplots()
#     # ax.plot(res['time'], [e * 1e-6 for e in res['totalPower']], label="GeN-Foam PK")
#     ax.plot(dataCsvCtrl['time'], [e * 1e-6 for e in dataCsvCtrl['pid.u_s']], label="Command")
#     for i in range(8):
#         dataCsvCore = pd.read_csv(f"CoreFMU_{i+1}.csv")
#         ax.plot(dataCsvCore['time'], [e * 1e-6 for e in dataCsvCore['gfPower']], label=f"GeN-Foam integralFMI {i+1}")
#     ax.set_xlabel("Time [s]")
#     ax.set_ylabel("Power [MW]")
#     ax.legend()
#     fig.tight_layout()
#     fig.savefig("fig_results_power.png")
#     plt.close()


#     fig, ax = plt.subplots()
#     ax.plot(res['time'], res["totalReactivity"], label="Total")
#     ax.plot(res['time'], res["extReactivity"], label="External")
#     ax.plot(res['time'], res["dopplerReactivity"], label="Doppler")
#     ax.plot(res['time'], res["TFuelReactivity"], label="TFuel")
#     ax.set_xlabel("Time [s]")
#     ax.set_ylabel("Reactivity [pcm]")
#     ax.legend()
#     fig.tight_layout()
#     fig.savefig("fig_results_reactivity.png")
#     plt.close()

#==============================================================================*
