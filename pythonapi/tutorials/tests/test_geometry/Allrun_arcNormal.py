import foamForNuclear as ffn
import foamForNuclear.mesh as mesh
import numpy as np


thMesh = mesh.BlockMesh(region='fluidRegion')

if (True):
    cube = thMesh.createCube(
        "cube",
        lowX=0, lowY=0, lowZ=0,
        highX=1, highY=1, highZ=1,
        nx=2, ny=3, nz=4,
        isAddBoundaryConditions=False
    )
    arc1 = thMesh.extrudeNormalArc(
        cube, 'right', 'arc1', arcCenter=ffn.Vector(1, 2, 0.5), angleSpan=45, nt=30
    )
    arc2 = thMesh.extrudeNormalArc(
        cube, 'left', 'arc2', arcCenter=ffn.Vector(0, 2, 0.5), angleSpan=-90, nt=30
    )
    arc3 = thMesh.extrudeNormalArc(
        cube, 'front', 'arc3', arcCenter=ffn.Vector(2, 0, 0.5), angleSpan=120, nt=30
    )
    arc4 = thMesh.extrudeNormalArc(
        cube, 'back', 'arc4', arcCenter=ffn.Vector(2, 1, 0.5), angleSpan=-90, nt=30
    )
    arc5 = thMesh.extrudeNormalArc(
        cube, 'top', 'arc5', arcCenter=ffn.Vector(0, 2, 1), angleSpan=-90, nt=10, rotationAxis='x'
    )
    arc6 = thMesh.extrudeNormalArc(
        arc5, 'top', 'arc6', arcCenter=ffn.Vector(2, 2, 1), angleSpan=-90, nt=10, rotationAxis='z'
    )
    arc7 = thMesh.extrudeNormalArc(
        cube, 'bottom', 'arc7', arcCenter=ffn.Vector(-2, 0, 0), angleSpan=-90, nt=10, rotationAxis='y'
    )

else:
    equivalentHydraulicDiameter = 0.3

    # Manifold
    manifoldBlocks = thMesh.createPipeCylindricalManifoldAlongZ(
        name='manifold',
        nEntries=3,
        innerRadius=1, outerRadius=2,
        equivalentHydraulicDiameter=equivalentHydraulicDiameter,
        lowZ=2,
        nr=4,
        nt=10
    )
    pipeEntries = manifoldBlocks[::2]

    for i, pipeEntry in enumerate(pipeEntries):
        pipe0 = thMesh.extrudeNormal(pipeEntry, 'front', f'hotLeg{i}_pipe0', 0.5, 2)

        direction = pipe0.getFaceNormal('front')
        direction.rotateZ(45*np.pi/180)

        pipe1 = thMesh.addPipe1DFromDirection(
            name=f'hotLeg{i}_pipe1',
            originPosition=pipe0,
            direction=direction,
            length=2,
            equivalentHydraulicDiameter=equivalentHydraulicDiameter,
            elbowRadius=0.5,
            n=4,
            isAddBoundaryConditions=True,
            originPositionOutletFaceName='front'
        )
        # pipe2 = thMesh.addPipe1DFromDirection(
        #     name=f'hotLeg{i}_pipe2',
        #     originPosition=pipe1,
        #     direction=ffn.Vector(0, 0, 1),
        #     length=3,
        #     equivalentHydraulicDiameter=equivalentHydraulicDiameter,
        #     elbowRadius=1,
        #     n=6,
        #     isAddBoundaryConditions=True
        # )

# thMesh.addMergePatchPairs()
# thMesh.mergePatchesWithName(name='outerClad', includeFacename=['OuterWall'])
# thMesh.mergePatchesWithName(name='topFuel', includeFacename=['cylinderTop_'])
# thMesh.mergePatchesWithName(name='bottomFuel', includeFacename=['cylinderBottom_'])

solver = ffn.ThermalHydraulicsSolver(region=thMesh.region, mesh=thMesh, solver='onePhase')


model = ffn.Model()
model.settings.application = 'dummy'

model.solvers.append(solver)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=thMesh, show_edges=True)
for normal in ['x', 'y', 'z']:
    model.plot_mesh(region=thMesh, show_edges=True, normal=normal)

# model.plot_boundary(region=thMesh, boundaryName='defaultFaces', show_edges=True)
# model.plot_boundary(region=thMesh, boundaryName='topFuel', show_edges=True)
# model.plot_boundary(region=thMesh, boundaryName='bottomFuel', show_edges=True)
# model.plot_boundary(region=thMesh, boundaryName='outerClad', show_edges=True)
