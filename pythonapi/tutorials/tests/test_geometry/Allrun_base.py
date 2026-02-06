import numpy as np
import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


def createDummyCase(nMesh: mesh.Mesh):
    solver = ffn.NeutronicsSolver(
        region=nMesh.region,
        mesh=nMesh,
        solver='diffusionNeutronics'
    )

    model = ffn.Model()
    model.settings.application = 'dummy'

    model.solvers.append(solver)

    ffn.cleanCase0()
    # ffn.allclean()
    model.export_to_openfoam()

    ffn.run_preprocessing(model=model)

    model.plot_mesh(region=nMesh, show_edges=True)
    # model.plot_mesh(region=nMesh, show_edges=True, normal='z')
    # try:
    #     model.plot_boundary(region=nMesh, boundaryName="wall", show_edges=True)
    # except:
    #     pass


def create_cube():
    nMesh = mesh.BlockMesh(region="cube")

    block1 = nMesh.create_cube(
        name='block',
        lowX=0, highX=1,
        lowY=0, highY=1,
        lowZ=0, highZ=1,
        nx=2, ny=3, nz=4
    )

    # block2 = nMesh.extrude_front(
    #     targetBlocks=block1,
    #     name="block2",
    #     dy=1,
    #     ny=2
    # )

    # p2 = block2.points[2]

    # p2.x = p2.x + 0.5

    return(nMesh)


def createCubeWithHoleCylz():
    nMesh = mesh.BlockMesh(region="square_hole_cylz")

    nMesh.create_cube_with_hole_along_z(
        "block",
        lowX=0, highX=0.01,
        lowY=0, highY=0.01,
        lowZ=0, highZ=0.01,
        radius=0.003,
        nx=4, ny=4, nz=2, nt=4,
        isAddAllBC=True
    )

    return(nMesh)


def createCubeWithHoleSqr():
    nMesh = mesh.BlockMesh(region="square_hole_sqr")

    nMesh.create_cube_with_hole_along_z(
        "block",
        lowX=0, highX=0.01,
        lowY=0, highY=0.01,
        lowZ=0, highZ=0.01,
        radius=0.003,
        nx=4, ny=4, nz=2, nt=4,
        isAddAllBC=True,
        isHoleCylinder=False
    )

    return(nMesh)


def create_wedge():
    nMesh = mesh.BlockMesh(region="wedge")

    nMesh.create_wedge(
        name='block',
        innerRadius=0,
        outerRadius=0.5,
        lowZ=0, highZ=1,
        wedgeAngle=15,
        nr=4,
        nz=5
    )

    return(nMesh)


def createHexagonPrism():
    nMesh = mesh.BlockMesh(region="hexagon")

    nMesh.create_hexagon_prism_along_z(
        name="hex",
        zmin=0,
        zmax=1,
        pitch=1,
        x=1, y=1,
        nr=2,
        nt=1,
        nz=3,
        isAddAllBC=True
    )
    nMesh.isMergeCoincidentPoints = True

    return(nMesh)


def createHexagonPrismFine():
    nMesh = mesh.BlockMesh(region="hexagon_fine")

    nMesh.create_hexagon_prism_along_z(
        name="hex",
        zmin=0,
        zmax=1,
        pitch=1,
        x=1, y=1,
        nr=3,
        nt=3,
        nz=4,
        isAddAllBC=True
    )
    nMesh.isMergeCoincidentPoints = True

    return(nMesh)


def createHexagonPrismWithHoleCylz():
    nMesh = mesh.BlockMesh(region="hexagon_hole_cylz")

    nMesh.create_hexagon_prism_with_hole_along_z(
        "hexHole",
        zmin=0, zmax=0.01,
        pitch=0.01,
        radius=0.003,
        nr=4, nt=4, nz=2,
        isAddAllBC=True
    )

    return(nMesh)


def createHexagonPrismWithHoleHex():
    nMesh = mesh.BlockMesh(region="hexagon_hole_hex")

    nMesh.create_hexagon_prism_with_hole_along_z(
        "hexHole",
        zmin=0, zmax=0.01,
        pitch=0.01,
        radius=0.003,
        nr=4, nt=4, nz=2,
        isAddAllBC=True,
        isHoleCylinder=False
    )

    return(nMesh)


def createCylinderZ():
    nMesh = mesh.BlockMesh(region="cylz")

    nMesh.create_cylinder_along_z(
        name='block',
        radius=0.5,
        lowZ=0, highZ=1,
        nx=3, ny=3, nz=4,
        x=1, y=-1
    )

    return(nMesh)


def createRingZ():
    nMesh = mesh.BlockMesh(region="ring")

    nMesh.create_ring_along_z(
        name='block',
        innerRadius=0.2, outerRadius=0.5,
        lowZ=0, highZ=1,
        nr=3, nt=10, nz=4,
        x=1, y=-1
    )

    return(nMesh)


def createRingSectorZ():
    nMesh = mesh.BlockMesh(region="ring_sector")

    nMesh.create_ring_sector_along_z(
        name='block',
        innerRadius=0.2, outerRadius=0.5,
        angleStart=0, angleArc=60,
        lowZ=0, highZ=1,
        nr=3, nt=5, nz=4,
        x=1, y=-1
    )

    return(nMesh)


def createTriangularChannel():
    nMesh = mesh.BlockMesh(region="hexagon_channel")

    nMesh.create_triangular_channel(
        name='block',
        lowZ=0, highZ=1,
        pitch=1,
        radius=0.3,
        nx=3, ny=3, nz=4
    )

    return(nMesh)


def create_sphere():
    nMesh = mesh.BlockMesh(region="sphere")

    nMesh.create_sphere(
        name="sphere",
        radius=1,
        nCenter=10,
        nBorder=2,
        isAddAllBC=True
    )
    nMesh.merge_patches_with_name(name="wall", includeFacename=[".*Wall.*"], patchType="wall")

    return(nMesh)


def create_hollow_half_sphere():
    nMesh = mesh.BlockMesh(region="hollow_half_sphere")

    nMesh.create_hollow_half_sphere(
        name="sph1",
        innerRadius=0.5, outerRadius=0.8,
        nr=3,
        nt=10,
        isAddAllBC=True
    )

    return(nMesh)


def create_half_sphere():
    nMesh = mesh.BlockMesh(region="half_sphere")

    nMesh.create_half_sphere(
        name="halfSphere",
        radius=1,
        z=-1,
        nCenter=10,
        nBorder=10,
        isAddAllBC=True
    )

    return(nMesh)


def create_sphere_1D():
    nMesh = mesh.BlockMesh(region="sphere_1D")

    nMesh.create_sphere_1D(
        name='core',
        innerRadius=0, outerRadius=1,
        opening=5,
        nr=10,
        isAddWedgeBC=True
    )
    # nMesh.create_sphere_1D(
    #     name='reflector',
    #     innerRadius=1, outerRadius=1.2,
    #     opening=5,
    #     nr=5,
    #     isAddWedgeBC=True,
    #     isAddOuterBC=True
    # )
    # nMesh.isMergeCoincidentPoints = True
    # nMesh.add_merge_patch_pairs(includeFacename=['Wall'], excludeFacename=['Wedge'])

    return(nMesh)


def createCubeWithCornerHole():
    nMesh = mesh.BlockMesh(region="cube_corner_hole")

    nMesh.create_cube_with_corner_hole_along_z(
        name="block",
        lowX=0, highX=1,
        lowY=0, highY=1,
        lowZ=0, highZ=1,
        radius=0.25,
        nx=3, ny=3, nz=4, nt=3,
        isHoleCylinder=True,
        isAddAllBC=True
    )

    # nMesh.create_ring_sector_along_z(
    #     name="sector",
    #     innerRadius=0,  outerRadius=0.25,
    #     angleStart=180, angleArc=45,
    #     lowZ=0, highZ=1,
    #     x=1, y=1,
    #     nr=2, nt=5, nz=4,
    #     isAddAllBC=True
    # )
    # nMesh.create_ring_sector_along_z(
    #     name="sector",
    #     innerRadius=0,      outerRadius=0.25,
    #     angleStart=180+45,  angleArc=45,
    #     lowZ=0, highZ=1,
    #     x=1, y=1,
    #     nr=2, nt=5, nz=4,
    #     isAddAllBC=True
    # )
    # nMesh.isMergeCoincidentPoints = True

    # externalFaces = nMesh.get_standalone_faces(excludeFacename=["Top", "Bottom"])
    # nMesh.merge_patches_with_name(name='top', includeFacename=['Top'])
    # nMesh.merge_patches_with_name(name='bottom', includeFacename=['Bottom'])
    # nMesh.merge_patches_with_name(name='wall', includeFacename=[f.name for f in externalFaces], isStrict=True)

    # ---

    # lxy = 1

    # xCenter, yCenter = 1, 1

    # r = lxy/2 * np.sqrt(2)
    # deg = np.pi/180

    # for theta in np.arange(45, 360+45, 90):
    #     rcost = r * np.cos(theta*deg)
    #     rsint = r * np.sin(theta*deg)
    #     nMesh.create_cube_with_corner_hole_along_z(
    #         name="block",
    #         lowX=xCenter+rcost-lxy/2, highX=xCenter+rcost+lxy/2,
    #         lowY=yCenter+rsint-lxy/2, highY=yCenter+rsint+lxy/2,
    #         lowZ=0, highZ=1,
    #         radius=0.25,
    #         nx=3, ny=3, nz=4, nt=3,
    #         isHoleCylinder=True,
    #         isAddAllBC=True,
    #         edgeFaceOrientation=theta-45
    #     )

    return(nMesh)


def createQuarterCylinderZ():
    nMesh = mesh.BlockMesh(region="quarter_cylinder")

    nMesh.create_quarter_cylinder_along_z(
        name="block",
        radius=1,
        lowZ=0, highZ=1,
        nx=3, ny=5, nz=3,
        angleStart=0,
        isAddAllBC=True
    )

    return(nMesh)



for funcMeshGen in [
    create_cube,
    createCubeWithCornerHole,
    createCubeWithHoleCylz,
    createCubeWithHoleSqr,
    create_wedge,
    createHexagonPrism,
    createHexagonPrismFine,
    createHexagonPrismWithHoleCylz,
    createHexagonPrismWithHoleHex,
    createCylinderZ,
    createQuarterCylinderZ,
    createRingZ,
    createRingSectorZ,
    createTriangularChannel,
    create_sphere,
    # create_hollow_half_sphere,
    create_half_sphere,
    create_sphere_1D
]:
    createDummyCase(nMesh=funcMeshGen())
