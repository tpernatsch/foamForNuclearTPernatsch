import numpy as np
import foamForNuclear as ffn
import foamForNuclear.mesh as mesh
import foamForNuclear.boundaryConditions as bc


isGenerateFromRandom: bool = False
isGenerateFromLattice: bool = True
isGenerateFromCoordinates: bool = False


nMesh = mesh.PolyMesh(region='neutroMesh')

if (isGenerateFromRandom):

    rng = np.random.default_rng(11)
    # create a set of points in 3D
    # nPoints = 20
    # points = [
    #     [p[0], p[1], p[2]] #, f'zone{i}']
    #     for i, p in enumerate(rng.uniform(size=(nPoints, 3)))
    # ]
    # points = rng.uniform(size=(nPoints, 3))

    points = [
        [0, 0, -1],
        [-1, 0, 0],
        [0, 0, 0],
        [1, 0, 0],
        [0, -1, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]

    nMesh.add_voronoi_points_from_coordinates(coords=points)

    boundaryBox = [[-1, 1], [-1, 1], [-1, 1]]


elif (isGenerateFromLattice):
    lattice = """
e e e e e e e e B B B B B B B B B
 e e e e e e e B R R R R R R R R B
  e e e e e e B R O C O O O C O R B
   e e e e e B R C O O O O O O C R B
    e e e e B R O O I I I I I O O R B
     e e e B R O O I C I I C I O O R B
      e e B R O O I I I I I I I O O R B
       e B R C O I I I I I I I I O C R B
        B R O O I C I I I I I C I O O R B
         B R C O I I I I I I I I O C R B e
          B R O O I I I I I I I O O R B e e
           B R O O I C I I C I O O R B e e e
            B R O O I I I I I O O R B e e e e
             B R C O O O O O O C R B e e e e e
              B R O C O O O C O R B e e e e e e
               B R R R R R R R R B e e e e e e e
                B B B B B B B B B e e e e e e e e
    """

    # lattice = """
    # e e R R R
    #  e R I I R
    #   R I I I R
    #    R I I R e
    #     R R R e e
    # """

    nXY = len([line for line in lattice.split('\n') if line.strip() != ""])

    edge = 0.12171
    pitch = edge * np.sqrt(3)

    fuelLength = 2
    diagridThickness = 0.2

    nMesh.add_voronoi_points_from_lattice(
        lattice=lattice,
        latticeType='hexagon',
        nx=nXY, ny=nXY,
        pitch=pitch,
        elementsToPlace=['I', 'O', 'C', 'R', 'B'] # B for border
    )

    latticeDiagrid = lattice
    for elementName in ['I', 'O', 'C', 'R', 'B']:
        latticeDiagrid = latticeDiagrid.replace(elementName, 'diagrid')

    nMesh.add_voronoi_points_from_lattice(
        lattice=latticeDiagrid,
        latticeType='hexagon',
        nx=nXY, ny=nXY,
        pitch=pitch,
        z=-fuelLength,
        elementsToPlace=['diagrid']
    )

    nMesh.add_voronoi_points_from_lattice(
        lattice=latticeDiagrid,
        latticeType='hexagon',
        nx=nXY, ny=nXY,
        pitch=pitch,
        z=fuelLength,
        elementsToPlace=['diagrid']
    )


elif (isGenerateFromCoordinates):
    nMesh.add_voronoi_points_from_coordinates(
        coords="""
    -17.320508075688775;-30;e
    -17.320508075688775;-20;e
    -17.320508075688775;-10;e
    -17.320508075688775;0;e
    -17.320508075688775;10;e
    -8.660254037844387;-25;e
    -8.660254037844387;-15;e
    -8.660254037844387;-5;F
    -8.660254037844387;5;F
    -8.660254037844387;15;e
    0;-20;e
    0;-10;F
    0;0;F
    0;10;F
    0;20;e
    8.660254037844387;-15;e
    8.660254037844387;-5;F
    8.660254037844387;5;F
    8.660254037844387;15;e
    8.660254037844387;25;e
    17.320508075688775;-10;e
    17.320508075688775;0;e
    17.320508075688775;10;e
    17.320508075688775;20;e
    17.320508075688775;30;e
    """,
        zLevel=0
    )

    boundaryBox = [[-40, 40], [-40, 40], [-40, 40]]


# Generate the mesh
nMesh.generate_mesh_from_voronoi_points(
    cellZonesToStrip=['diagrid', 'B'],
    # isAddInletOuletBC=True,
    # isAddBaffles=True
)

# for i, vec in enumerate(nMesh.pointsList):
#     if (abs(vec.z + 0.55) <= 1e-6):
#         nMesh.pointsList[i].z = -fuelLength/2

# Plot the Voronoi graph
# import matplotlib.pyplot as plt
# fig, ax = nMesh.plot_voronoi_cells()
# plt.show()


timeFolder0 = ffn.TimeFolder(time=0)

defaultFlux = ffn.Field(
    name="defaultFlux",
    dimensions=ffn.Dimension(default='flux'),
    region=nMesh.region
)
defaultFlux.internalField = 1
defaultFlux.set_boundary_condition('inlet', bc.FixedValue(0))
defaultFlux.set_boundary_condition('outlet', bc.FixedValue(0))
defaultFlux.set_boundary_condition('defaultFaces', bc.FixedValue(0))

timeFolder0.append(defaultFlux)

solver = ffn.NeutronicsSolver(
    region=nMesh.region,
    mesh=nMesh,
    solver='diffusionNeutronics'
)


model = ffn.Model(timeFolders=[timeFolder0])
model.settings.application = 'dummy'

model.solvers.append(solver)

# Add idxField to visualize cellZones, not used in calculations
idxField = solver.create_zone_field(nMesh.cellZones)
idxField.set_boundary_condition('inlet', bc.ZeroGradient())
idxField.set_boundary_condition('outlet', bc.ZeroGradient())
timeFolder0.append(idxField)

print(model)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=nMesh, show_edges=True)

model.plot_mesh(region=nMesh, show_edges=True, fieldName=idxField.name, cmap="tab10", limits=[0.5, 10.5])
model.plot_mesh(region=nMesh, show_edges=True, fieldName=idxField.name, cmap="tab10", limits=[0.5, 10.5], normal='x')

for boundary in ['baffleZone_master', 'inlet', 'outlet', 'defaultFaces']:
    try:
        model.plot_boundary(region=nMesh, boundaryName=boundary)
    except:
        pass
