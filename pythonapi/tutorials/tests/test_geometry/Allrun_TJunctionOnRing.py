import numpy as np
import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


nMesh = mesh.BlockMesh(region='neutroRegion')



def createTJunctionOnRing(
        name: str,
        innerRadius: float,
        outerRadius: float,
        theta: float,
        thetaOpening: float,
        pipeRadius: float,
        pipeLength: float,
        lowZ: float,
        highZ: float,
        x: float=0, y: float=0,
        nr: int=1,
        nt: int=1,
        nz: int=1
    ) -> mesh.Block:

    sqrt2 = np.sqrt(2)

    pipeDiameter = 2*pipeRadius
    pipeDiameter3 = pipeDiameter/3
    dPipeEnd = pipeLength + outerRadius

    pipeApperture3 = np.pi * pipeDiameter3 / (2*outerRadius) /2
    pipeApperture23 = np.pi * (pipeRadius/sqrt2 - pipeDiameter3/2) / (2*outerRadius) / sqrt2

    thetaHalf = theta + thetaOpening/2
    theta15 = thetaHalf - pipeApperture3/2 - pipeApperture23
    theta25 = thetaHalf - pipeApperture3/2
    theta35 = thetaHalf + pipeApperture3/2
    theta45 = theta35 + pipeApperture23


    dz2 = (highZ-lowZ)/2
    z15 = lowZ + dz2 - pipeRadius/sqrt2
    z25 = lowZ + dz2 - pipeDiameter3/2
    z35 = highZ - dz2 + pipeDiameter3/2
    z45 = highZ - dz2 + pipeRadius/sqrt2
    zPipe = (highZ+lowZ)/2

    cost, sint = np.cos(theta), np.sin(theta)
    costhalf, sinthalf = np.cos(thetaHalf), np.sin(thetaHalf)
    cost15, sint15 = np.cos(theta15), np.sin(theta15)
    cost25, sint25 = np.cos(theta25), np.sin(theta25)
    cost35, sint35 = np.cos(theta35), np.sin(theta35)
    cost45, sint45 = np.cos(theta45), np.sin(theta45)
    costt, sintt = np.cos(theta+thetaOpening), np.sin(theta+thetaOpening)

    blockA1 = nMesh.create_block(name, [
        mesh.Point(x+outerRadius*cost, y+outerRadius*sint, lowZ),
        mesh.Point(x+outerRadius*cost15, y+outerRadius*sint15, lowZ),
        mesh.Point(x+innerRadius*cost15, y+innerRadius*sint15, lowZ),
        mesh.Point(x+innerRadius*cost, y+innerRadius*sint, lowZ),
        mesh.Point(x+outerRadius*cost, y+outerRadius*sint, z15),
        mesh.Point(x+outerRadius*cost15, y+outerRadius*sint15, z15),
        mesh.Point(x+innerRadius*cost15, y+innerRadius*sint15, z15),
        mesh.Point(x+innerRadius*cost, y+innerRadius*sint, z15),
    ], nx=nt, ny=nr, nz=nz)

    blockA2 = nMesh.create_block(name, [
        blockA1.points[1],
        mesh.Point(x+outerRadius*cost45, y+outerRadius*sint45, lowZ),
        mesh.Point(x+innerRadius*cost45, y+innerRadius*sint45, lowZ),
        blockA1.points[2],
        blockA1.points[5],
        mesh.Point(x+outerRadius*cost45, y+outerRadius*sint45, z15),
        mesh.Point(x+innerRadius*cost45, y+innerRadius*sint45, z15),
        blockA1.points[6],
    ], nx=nt, ny=nr, nz=nz)

    blockA3 = nMesh.create_block(name, [
        blockA2.points[1],
        mesh.Point(x+outerRadius*costt, y+outerRadius*sintt, lowZ),
        mesh.Point(x+innerRadius*costt, y+innerRadius*sintt, lowZ),
        blockA2.points[2],
        blockA2.points[5],
        mesh.Point(x+outerRadius*costt, y+outerRadius*sintt, z15),
        mesh.Point(x+innerRadius*costt, y+innerRadius*sintt, z15),
        blockA2.points[6],
    ], nx=nt, ny=nr, nz=nz)

    # blockA4 = nMesh.create_block(name, [
    #     blockA3.points[1],
    #     mesh.Point(x+outerRadius*cost45, y+outerRadius*sint45, lowZ),
    #     mesh.Point(x+innerRadius*cost45, y+innerRadius*sint45, lowZ),
    #     blockA3.points[2],
    #     blockA3.points[5],
    #     mesh.Point(x+outerRadius*cost45, y+outerRadius*sint45, z15),
    #     mesh.Point(x+innerRadius*cost45, y+innerRadius*sint45, z15),
    #     blockA3.points[6],
    # ], nx=nt, ny=nr, nz=nz)

    # blockA5 = nMesh.create_block(name, [
    #     blockA4.points[1],
    #     mesh.Point(x+outerRadius*costt, y+outerRadius*sintt, lowZ),
    #     mesh.Point(x+innerRadius*costt, y+innerRadius*sintt, lowZ),
    #     blockA4.points[2],
    #     blockA4.points[5],
    #     mesh.Point(x+outerRadius*costt, y+outerRadius*sintt, z15),
    #     mesh.Point(x+innerRadius*costt, y+innerRadius*sintt, z15),
    #     blockA4.points[6],
    # ], nx=nt, ny=nr, nz=nz)

    (blockB1, blockB3) = nMesh.extrude_top(
        [blockA1, blockA3],
        name,
        dz=z45-z15,
        nz=nz
    )

    (blockC1, blockC3) = nMesh.extrude_top(
        [blockB1, blockB3],
        name,
        dz=highZ-z45,
        nz=nz
    )

    blockC2 = nMesh.create_block(name, [
        blockC1.points[1],
        blockC3.points[0],
        blockC3.points[3],
        blockC1.points[2],
        blockC1.points[5],
        blockC3.points[4],
        blockC3.points[7],
        blockC1.points[6],
    ], nx=blockA1.nx, ny=blockA1.ny, nz=nz)

    # (blockD1, blockD2, blockD3, blockD4, blockD5) = nMesh.extrude_top(
    #     [blockC1, blockC2, blockC3, blockC4, blockC5],
    #     name,
    #     dz=z45-z35,
    #     nz=nz
    # )

    # (blockE1, blockE2, blockE3, blockE4, blockE5) = nMesh.extrude_top(
    #     [blockD1, blockD2, blockD3, blockD4, blockD5],
    #     name,
    #     dz=highZ-z45,
    #     nz=nz
    # )

    blockPC = nMesh.create_block(name, [
        mesh.Point(x+outerRadius*cost25, y+outerRadius*sint25, z25),
        mesh.Point(x+outerRadius*cost35, y+outerRadius*sint35, z25),
        mesh.Point(x+innerRadius*cost35, y+innerRadius*sint35, z25),
        mesh.Point(x+innerRadius*cost25, y+innerRadius*sint25, z25),
        mesh.Point(x+outerRadius*cost25, y+outerRadius*sint25, z35),
        mesh.Point(x+outerRadius*cost35, y+outerRadius*sint35, z35),
        mesh.Point(x+innerRadius*cost35, y+innerRadius*sint35, z35),
        mesh.Point(x+innerRadius*cost25, y+innerRadius*sint25, z35),
    ], nx=blockA1.nx, ny=blockA1.ny, nz=nz)

    blockPBottom = nMesh.add_top(blockA2, name, [
        blockPC.points[0],
        blockPC.points[1],
        blockPC.points[2],
        blockPC.points[3],
    ], nz=5)

    blockPTop = nMesh.add_bottom(blockC2, name, blockPC.topFace(), nz=blockPBottom.nz)
    blockPLeft = nMesh.add_right(blockB1, name, [
        blockPC.points[0],
        blockPC.points[3],
        blockPC.points[4],
        blockPC.points[7],
    ], nx=blockPBottom.nz)
    blockPRight = nMesh.add_left(blockB3, name, [
        blockPC.points[1],
        blockPC.points[2],
        blockPC.points[5],
        blockPC.points[6],
    ], nx=blockPBottom.nz)

    # Add arcs
    for block in [
        blockA1, blockA2, blockA3,
        blockB1, blockB3,
        blockC1, blockC3
    ]:
        block.addEdge('arc', 0, 1, x, y, isOrigin=True)
        block.addEdge('arc', 2, 3, x, y, isOrigin=True)

    for block in [blockC1, blockC2, blockC3]:
        block.addEdge('arc', 4, 5, x, y, isOrigin=True)
        block.addEdge('arc', 6, 7, x, y, isOrigin=True)


    blockPC.addEdge('arc', 2, 3, x, y, isOrigin=True)
    blockPC.addEdge('arc', 6, 7, x, y, isOrigin=True)
    blockA2.addEdge('arc', 6, 7, x, y, isOrigin=True)
    blockC2.addEdge('arc', 2, 3, x, y, isOrigin=True)

    # Pipe
    xPipe, yPipe = x + dPipeEnd*costhalf, y + dPipeEnd*sinthalf
    pipeCenter = nMesh.add_front(blockPC, name, [
        mesh.Point(xPipe + pipeDiameter3/2*sinthalf, yPipe - pipeDiameter3/2*costhalf, z25),
        mesh.Point(xPipe - pipeDiameter3/2*sinthalf, yPipe + pipeDiameter3/2*costhalf, z25),
        mesh.Point(xPipe + pipeDiameter3/2*sinthalf, yPipe - pipeDiameter3/2*costhalf, z35),
        mesh.Point(xPipe - pipeDiameter3/2*sinthalf, yPipe + pipeDiameter3/2*costhalf, z35),
    ], 5)

    pipeBottom = nMesh.add_front(blockPBottom, name, [
        mesh.Point(xPipe + pipeRadius/sqrt2*sinthalf, yPipe - pipeRadius/sqrt2*costhalf, z15),
        mesh.Point(xPipe - pipeRadius/sqrt2*sinthalf, yPipe + pipeRadius/sqrt2*costhalf, z15),
        pipeCenter.points[0],
        pipeCenter.points[1],
    ], ny=pipeCenter.ny)

    pipeTop = nMesh.add_front(blockPTop, name, [
        pipeCenter.points[4],
        pipeCenter.points[5],
        mesh.Point(xPipe + pipeRadius/sqrt2*sinthalf, yPipe - pipeRadius/sqrt2*costhalf, z45),
        mesh.Point(xPipe - pipeRadius/sqrt2*sinthalf, yPipe + pipeRadius/sqrt2*costhalf, z45),
    ], ny=pipeCenter.ny)

    pipeLeft = nMesh.add_front(blockPLeft, name, [
        pipeBottom.points[0],
        pipeBottom.points[4],
        pipeTop.points[4],
        pipeTop.points[0],
    ], ny=pipeCenter.ny)

    pipeRight = nMesh.add_front(blockPRight, name, [
        pipeBottom.points[5],
        pipeBottom.points[1],
        pipeTop.points[1],
        pipeTop.points[5],
    ], ny=pipeCenter.ny)

    # Pipe end
    pipeBottom.addEdge('arc', 0, 1, xPipe, yPipe, zPipe, isOrigin=True)
    pipeTop.addEdge('arc', 4, 5, xPipe, yPipe, zPipe, isOrigin=True)
    pipeLeft.addEdge('arc', 0, 4, xPipe, yPipe, zPipe, isOrigin=True)
    pipeRight.addEdge('arc', 5, 1, xPipe, yPipe, zPipe, isOrigin=True)

    nMesh.add_cylinder(
        name=name+"_outerRadius",
        point1=(x, y, lowZ-1),
        point2=(x, y, highZ+1),
        radius=outerRadius
    )
    nMesh.add_cylinder(
        name=name+"_cylinderPipe",
        point1=(x, y, zPipe),
        point2=(xPipe, yPipe, zPipe),
        radius=pipeRadius
    )
    pipeBottom.add_face_projection("bottom", name+"_cylinderPipe")
    pipeTop.add_face_projection("top", name+"_cylinderPipe")
    pipeLeft.add_face_projection("left", name+"_cylinderPipe")
    pipeRight.add_face_projection("right", name+"_cylinderPipe")

    # blockA1.add_edge_projection(4, 5, ["outerRadius"])
    # blockA3.add_edge_projection(4, 5, ["outerRadius"])
    # blockB1.add_edge_projection(4, 5, ["outerRadius"])
    # blockB3.add_edge_projection(4, 5, ["outerRadius"])

    pipeBottom.add_edge_projection(2, 3, [name+"_cylinderPipe", name+"_outerRadius"])
    pipeTop.add_edge_projection(6, 7, [name+"_cylinderPipe", name+"_outerRadius"])
    pipeLeft.add_edge_projection(3, 7, [name+"_cylinderPipe", name+"_outerRadius"])
    pipeRight.add_edge_projection(2, 6, [name+"_cylinderPipe", name+"_outerRadius"])


    return()


name = "block"

innerRadius = 1
outerRadius = 1.2
pipeRadius = 0.2
pipeLength = 0.2
theta = -np.pi/3
thetaOpening = np.pi/4
lowZ = 0
highZ = 2

createTJunctionOnRing(
    f"{name}_{0}",
    innerRadius, outerRadius,
    np.pi/3, thetaOpening,
    pipeRadius, pipeLength,
    lowZ, highZ,
    x=0, y=0,
    nr=1, nt=10, nz=10
)


# for i, theta in enumerate(np.arange(0, 2*np.pi, thetaOpening)):
#     createTJunctionOnRing(
#         f"{name}_{i}",
#         innerRadius, outerRadius,
#         theta, thetaOpening,
#         pipeRadius, pipeLength,
#         lowZ, highZ,
#         x=0, y=0,
#         nr=1, nt=10, nz=10
#     )



solver = ffn.NeutronicsSolver(
    region=nMesh.region,
    mesh=nMesh,
    solver='diffusionNeutronics'
)


model = ffn.Model()
model.settings.application = 'dummy'

model.solvers.append(solver)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=nMesh, show_edges=True)
try:
    model.plot_boundary(region=nMesh, boundaryName="wall", show_edges=True)
except:
    pass
