"""

"""
#==============================================================================*
# Imports

import numpy as np
import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh

def thermalHydraulicMesh(region):

    crocus = ffn.BlockMesh(region=region)

    crocus.isReducedCells = False

    coreRadius = 1.3/2
    expansionRadius = 0.8
    vesselHalfLength = 0.9
    coreHeight = 1.2
    coreBaseHeight = 0.15
    opening = np.pi/2 - np.pi/2.5

    nzCore = 30
    nzExpansionTanks = 30

    sqrt2 = np.sqrt(2)


    #---------------------------------------------------------------------------
    # Core

    coreCenter = crocus.createCube("core", -coreRadius/6, -coreRadius/6, 0, coreRadius/6, coreRadius/6, coreHeight, 3, 3, nzCore)

    coreCenterLeft = crocus.extrudeLeft([coreCenter], "core", coreRadius/3, nx=coreCenter.nx)
    coreCenterRight = crocus.extrudeRight([coreCenter], "core", coreRadius/3, nx=coreCenter.nx)

    coreFrontLeft, coreCenterFront, coreFrontRight = crocus.extrudeFront(
        [coreCenterLeft, coreCenter, coreCenterRight], "core", coreRadius/3, ny=coreCenter.ny
    )

    coreBackLeft, coreCenterBack, coreBackRight = crocus.extrudeBack(
        [coreCenterLeft, coreCenter, coreCenterRight], "core", coreRadius/3, ny=coreCenter.ny
    )


    coreFront = crocus.addFront(coreCenterFront, "coreFront", [
        mesh.Point(-coreRadius*np.sin(opening), -coreRadius*np.cos(opening), 0),
        mesh.Point(+coreRadius*np.sin(opening), -coreRadius*np.cos(opening), 0),
        mesh.Point(-coreRadius*np.sin(opening), -coreRadius*np.cos(opening), coreHeight),
        mesh.Point(+coreRadius*np.sin(opening), -coreRadius*np.cos(opening), coreHeight),
    ], ny=4)

    coreFront.addEdge("arc", 0, 1, x=0, y=-coreRadius)
    coreFront.addEdge("arc", 4, 5, x=0, y=-coreRadius)

    ffn.PointBottomRight0 = mesh.Point(coreRadius/sqrt2, -coreRadius/sqrt2, 0)
    ffn.PointBottomRightH = mesh.Point(coreRadius/sqrt2, -coreRadius/sqrt2, coreHeight)

    coreFront1 = crocus.addRight(coreFront, "coreFront", [
        ffn.PointBottomRight0,
        coreFrontRight.points[1],
        ffn.PointBottomRightH,
        coreFrontRight.points[5],
    ], nx=coreFrontRight.nx)

    coreFront1.addEdge("arc", 0, 1, x=coreRadius*np.sin(2*opening), y=-coreRadius*np.cos(2*opening))
    coreFront1.addEdge("arc", 4, 5, x=coreRadius*np.sin(2*opening), y=-coreRadius*np.cos(2*opening))

    ffn.PointBottomLeft0 = mesh.Point(-coreRadius/sqrt2, -coreRadius/sqrt2, 0)
    ffn.PointBottomLeftH = mesh.Point(-coreRadius/sqrt2, -coreRadius/sqrt2, coreHeight)

    coreFront2 = crocus.addLeft(coreFront, "coreFront", [
        ffn.PointBottomLeft0,
        coreFrontLeft.points[0],
        ffn.PointBottomLeftH,
        coreFrontLeft.points[4],
    ], nx=coreFrontLeft.nx)

    coreFront2.addEdge("arc", 0, 1, x=-coreRadius*np.sin(2*opening), y=-coreRadius*np.cos(2*opening))
    coreFront2.addEdge("arc", 4, 5, x=-coreRadius*np.sin(2*opening), y=-coreRadius*np.cos(2*opening))

    coreRight = crocus.addRight(coreCenterRight, "coreRight", [
        mesh.Point(+coreRadius*np.cos(opening), -coreRadius*np.sin(opening), 0),
        mesh.Point(+coreRadius*np.cos(opening), +coreRadius*np.sin(opening), 0),
        mesh.Point(+coreRadius*np.cos(opening), -coreRadius*np.sin(opening), coreHeight),
        mesh.Point(+coreRadius*np.cos(opening), +coreRadius*np.sin(opening), coreHeight),
    ], nx=coreFront.ny)

    coreRight.addEdge("arc", 1, 2, x=coreRadius, y=0)
    coreRight.addEdge("arc", 5, 6, x=coreRadius, y=0)

    ffn.PointTopRight0 = mesh.Point(coreRadius/sqrt2, coreRadius/sqrt2, 0)
    ffn.PointTopRightH = mesh.Point(coreRadius/sqrt2, coreRadius/sqrt2, coreHeight)

    coreRight1 = crocus.addBack(coreRight, "coreRight", [
        ffn.PointTopRight0,
        coreBackRight.points[2],
        ffn.PointTopRightH,
        coreBackRight.points[6],
    ], ny=coreBackRight.ny)

    coreRight1.addEdge("arc", 1, 2, x=coreRadius*np.cos(2*opening), y=coreRadius*np.sin(2*opening))
    coreRight1.addEdge("arc", 5, 6, x=coreRadius*np.cos(2*opening), y=coreRadius*np.sin(2*opening))

    coreRight2 = crocus.addFront(coreRight, "coreRight", [
        coreFrontRight.points[1],
        ffn.PointBottomRight0,
        coreFrontRight.points[5],
        ffn.PointBottomRightH,
    ], ny=coreFrontRight.ny)

    coreRight2.addEdge("arc", 1, 2, x=coreRadius*np.cos(2*opening), y=-coreRadius*np.sin(2*opening))
    coreRight2.addEdge("arc", 5, 6, x=coreRadius*np.cos(2*opening), y=-coreRadius*np.sin(2*opening))

    coreLeft = crocus.addLeft(coreCenterLeft, "coreLeft", [
        mesh.Point(-coreRadius*np.cos(opening), -coreRadius*np.sin(opening), 0),
        mesh.Point(-coreRadius*np.cos(opening), +coreRadius*np.sin(opening), 0),
        mesh.Point(-coreRadius*np.cos(opening), -coreRadius*np.sin(opening), coreHeight),
        mesh.Point(-coreRadius*np.cos(opening), +coreRadius*np.sin(opening), coreHeight),
    ], nx=coreFront.ny)

    coreLeft.addEdge("arc", 3, 0, x=-coreRadius, y=0)
    coreLeft.addEdge("arc", 7, 4, x=-coreRadius, y=0)

    coreLeft1 = crocus.addFront(coreLeft, "coreLeft", [
        ffn.PointBottomLeft0,
        coreFrontLeft.points[0],
        ffn.PointBottomLeftH,
        coreFrontLeft.points[4],
    ], ny=coreFrontLeft.ny)

    coreLeft1.addEdge("arc", 3, 0, x=-coreRadius*np.cos(2*opening), y=-coreRadius*np.sin(2*opening))
    coreLeft1.addEdge("arc", 7, 4, x=-coreRadius*np.cos(2*opening), y=-coreRadius*np.sin(2*opening))

    ffn.PointTopLeft0 = mesh.Point(-coreRadius/sqrt2, coreRadius/sqrt2, 0)
    ffn.PointTopLeftH = mesh.Point(-coreRadius/sqrt2, coreRadius/sqrt2, coreHeight)

    coreLeft2 = crocus.addBack(coreLeft, "coreLeft", [
        coreBackLeft.points[3],
        ffn.PointTopLeft0,
        coreBackLeft.points[7],
        ffn.PointTopLeftH,
    ], ny=coreFrontRight.ny)

    coreLeft2.addEdge("arc", 3, 0, x=-coreRadius*np.cos(2*opening), y=coreRadius*np.sin(2*opening))
    coreLeft2.addEdge("arc", 7, 4, x=-coreRadius*np.cos(2*opening), y=coreRadius*np.sin(2*opening))

    coreBack = crocus.addBack(coreCenterBack, "coreBack", [
        mesh.Point(+coreRadius*np.sin(opening), +coreRadius*np.cos(opening), 0),
        mesh.Point(-coreRadius*np.sin(opening), +coreRadius*np.cos(opening), 0),
        mesh.Point(+coreRadius*np.sin(opening), +coreRadius*np.cos(opening), coreHeight),
        mesh.Point(-coreRadius*np.sin(opening), +coreRadius*np.cos(opening), coreHeight),
    ], ny=coreFront.ny)

    coreBack.addEdge("arc", 2, 3, x=0, y=coreRadius)
    coreBack.addEdge("arc", 6, 7, x=0, y=coreRadius)

    coreBack1 = crocus.addLeft(coreBack, "coreBack", [
        coreBackLeft.points[3],
        ffn.PointTopLeft0,
        coreBackLeft.points[7],
        ffn.PointTopLeftH,
    ], nx=coreFrontLeft.nx)

    coreBack1.addEdge("arc", 2, 3, x=-coreRadius*np.sin(2*opening), y=coreRadius*np.cos(2*opening))
    coreBack1.addEdge("arc", 6, 7, x=-coreRadius*np.sin(2*opening), y=coreRadius*np.cos(2*opening))

    coreBack2 = crocus.addRight(coreBack, "coreBack", [
        coreBackRight.points[2],
        ffn.PointTopRight0,
        coreBackRight.points[6],
        ffn.PointTopRightH,
    ], nx=coreBackRight.nx)

    coreBack2.addEdge("arc", 2, 3, x=coreRadius*np.sin(2*opening), y=coreRadius*np.cos(2*opening))
    coreBack2.addEdge("arc", 6, 7, x=coreRadius*np.sin(2*opening), y=coreRadius*np.cos(2*opening))


    (
        baseCenterLeft, baseCenterBack, baseCenterRight,
        baseBackLeft, baseBackRight, baseFrontLeft, baseFrontRight,
        baseFront1, baseFront2,
        baseLeft, baseLeft1, baseLeft2,
        baseRight, baseRight1, baseRight2,
        baseBack, baseBack1, baseBack2
    ) = crocus.extrudeBottom([
        coreCenterLeft, coreCenterBack, coreCenterRight,
        coreBackLeft, coreBackRight, coreFrontLeft, coreFrontRight,
        coreFront1, coreFront2,
        coreLeft, coreLeft1, coreLeft2,
        coreRight, coreRight1, coreRight2,
        coreBack, coreBack1, coreBack2,
    ], "coreBase", coreBaseHeight, nz=5)

    #---------------------------------------------------------------------------
    # Expansion back right

    channelBackRight1 = crocus.addRight(baseRight1, "channelBackRight", [
        mesh.Point(expansionRadius*np.cos(opening), expansionRadius*np.sin(opening), -coreBaseHeight),
        mesh.Point(expansionRadius/sqrt2, expansionRadius/sqrt2, -coreBaseHeight),
        mesh.Point(expansionRadius*np.cos(opening), expansionRadius*np.sin(opening), 0),
        mesh.Point(expansionRadius/sqrt2, expansionRadius/sqrt2, 0),
    ], nx=2)

    expansionBaseBackRight1 = crocus.addRight(channelBackRight1, "expansionBaseBackRight", [
        mesh.Point(vesselHalfLength, expansionRadius*np.sin(opening), -coreBaseHeight),
        mesh.Point(vesselHalfLength, vesselHalfLength, -coreBaseHeight),
        mesh.Point(vesselHalfLength, expansionRadius*np.sin(opening), 0),
        mesh.Point(vesselHalfLength, vesselHalfLength, 0),
    ], nx=3)

    expansionBaseBackRight2 = crocus.addFront(expansionBaseBackRight1, "expansionBaseBackRight", [
        mesh.Point(expansionRadius*np.cos(opening/2), expansionRadius*np.sin(opening/2), -coreBaseHeight),
        mesh.Point(vesselHalfLength, expansionRadius*np.sin(opening/2), -coreBaseHeight),
        mesh.Point(expansionRadius*np.cos(opening/2), expansionRadius*np.sin(opening/2), 0),
        mesh.Point(vesselHalfLength, expansionRadius*np.sin(opening/2), 0),
    ], ny=4)


    channelBackRight2 = crocus.addBack(baseBack2, "channelBackRight", [
        channelBackRight1.points[2],
        mesh.Point(expansionRadius*np.sin(opening), expansionRadius*np.cos(opening), -coreBaseHeight),
        channelBackRight1.points[6],
        mesh.Point(expansionRadius*np.sin(opening), expansionRadius*np.cos(opening), 0),
    ], ny=channelBackRight1.nx)

    expansionBaseBackRight3 = crocus.addBack(channelBackRight2, "expansionBaseBackRight", [
        expansionBaseBackRight1.points[2],
        mesh.Point(expansionRadius*np.sin(opening), vesselHalfLength, -coreBaseHeight),
        expansionBaseBackRight1.points[6],
        mesh.Point(expansionRadius*np.sin(opening), vesselHalfLength, 0)
    ], ny=expansionBaseBackRight1.nx)

    expansionBaseBackRight4 = crocus.addLeft(expansionBaseBackRight3, "expansionBaseBackRight", [
        mesh.Point(expansionRadius*np.sin(opening/2), expansionRadius*np.cos(opening/2), -coreBaseHeight),
        mesh.Point(expansionRadius*np.sin(opening/2), vesselHalfLength, -coreBaseHeight),
        mesh.Point(expansionRadius*np.sin(opening/2), expansionRadius*np.cos(opening/2), 0),
        mesh.Point(expansionRadius*np.sin(opening/2), vesselHalfLength, 0),
    ], nx=expansionBaseBackRight2.ny)


    (
        expansionBackRight1, expansionBackRight2,
        expansionBackRight3, expansionBackRight4
    ) = crocus.extrudeTop([
        expansionBaseBackRight1, expansionBaseBackRight2,
        expansionBaseBackRight3, expansionBaseBackRight4
    ], "expansionBackRight", coreHeight, nz=nzExpansionTanks)

    expansionBackRight1.addEdge("arc", 3, 0, x=expansionRadius*np.cos(opening/2), y=expansionRadius*np.sin(opening/2))
    expansionBackRight1.addEdge("arc", 7, 4, x=expansionRadius*np.cos(opening/2), y=expansionRadius*np.sin(opening/2))
    expansionBackRight3.addEdge("arc", 0, 1, x=expansionRadius*np.cos(opening/2), y=expansionRadius*np.sin(opening/2))
    expansionBackRight3.addEdge("arc", 4, 5, x=expansionRadius*np.cos(opening/2), y=expansionRadius*np.sin(opening/2))


    #---------------------------------------------------------------------------
    # Expansion back left


    channelBackLeft1 = crocus.addLeft(baseLeft2, "channelBackLeft", [
        mesh.Point(-expansionRadius*np.cos(opening), expansionRadius*np.sin(opening), -coreBaseHeight),
        mesh.Point(-expansionRadius/sqrt2, expansionRadius/sqrt2, -coreBaseHeight),
        mesh.Point(-expansionRadius*np.cos(opening), expansionRadius*np.sin(opening), 0),
        mesh.Point(-expansionRadius/sqrt2, expansionRadius/sqrt2, 0),
    ], nx=channelBackRight1.nx)

    expansionBaseBackLeft1 = crocus.addLeft(channelBackLeft1, "expansionBaseBackLeft", [
        mesh.Point(-vesselHalfLength, expansionRadius*np.sin(opening), -coreBaseHeight),
        mesh.Point(-vesselHalfLength, vesselHalfLength, -coreBaseHeight),
        mesh.Point(-vesselHalfLength, expansionRadius*np.sin(opening), 0),
        mesh.Point(-vesselHalfLength, vesselHalfLength, 0),
    ], nx=expansionBaseBackRight1.nx)

    expansionBaseBackLeft2 = crocus.addFront(expansionBaseBackLeft1, "expansionBaseBackLeft", [
        mesh.Point(-vesselHalfLength, expansionRadius*np.sin(opening/2), -coreBaseHeight),
        mesh.Point(-expansionRadius*np.cos(opening/2), expansionRadius*np.sin(opening/2), -coreBaseHeight),
        mesh.Point(-vesselHalfLength, expansionRadius*np.sin(opening/2), 0),
        mesh.Point(-expansionRadius*np.cos(opening/2), expansionRadius*np.sin(opening/2), 0),
    ], ny=expansionBaseBackRight2.ny)


    channelBackLeft2 = crocus.addBack(baseBack1, "channelBackLeft", [
        mesh.Point(-expansionRadius*np.sin(opening), expansionRadius*np.cos(opening), -coreBaseHeight),
        channelBackLeft1.points[3],
        mesh.Point(-expansionRadius*np.sin(opening), expansionRadius*np.cos(opening), 0),
        channelBackLeft1.points[7],
    ], ny=channelBackLeft1.nx)

    expansionBaseBackLeft3 = crocus.addBack(channelBackLeft2, "expansionBaseBackLeft", [
        mesh.Point(-expansionRadius*np.sin(opening), vesselHalfLength, -coreBaseHeight),
        expansionBaseBackLeft1.points[3],
        mesh.Point(-expansionRadius*np.sin(opening), vesselHalfLength, 0),
        expansionBaseBackLeft1.points[7],
    ], ny=expansionBaseBackLeft1.nx)

    expansionBaseBackLeft4 = crocus.addRight(expansionBaseBackLeft3, "expansionBaseBackLeft", [
        mesh.Point(-expansionRadius*np.sin(opening/2), expansionRadius*np.cos(opening/2), -coreBaseHeight),
        mesh.Point(-expansionRadius*np.sin(opening/2), vesselHalfLength, -coreBaseHeight),
        mesh.Point(-expansionRadius*np.sin(opening/2), expansionRadius*np.cos(opening/2), 0),
        mesh.Point(-expansionRadius*np.sin(opening/2), vesselHalfLength, 0),
    ], nx=expansionBaseBackLeft2.ny)



    (
        expansionBackLeft1, expansionBackLeft2,
        expansionBackLeft3, expansionBackLeft4
    ) = crocus.extrudeTop([
        expansionBaseBackLeft1, expansionBaseBackLeft2,
        expansionBaseBackLeft3, expansionBaseBackLeft4
    ], "expansionBackLeft", coreHeight, nz=expansionBackRight1.nz)

    expansionBackLeft1.addEdge("arc", 1, 2, x=-expansionRadius*np.sin(opening/2), y=expansionRadius*np.cos(opening/2))
    expansionBackLeft1.addEdge("arc", 5, 6, x=-expansionRadius*np.sin(opening/2), y=expansionRadius*np.cos(opening/2))
    expansionBackLeft3.addEdge("arc", 0, 1, x=-expansionRadius*np.sin(opening/2), y=expansionRadius*np.cos(opening/2))
    expansionBackLeft3.addEdge("arc", 4, 5, x=-expansionRadius*np.sin(opening/2), y=expansionRadius*np.cos(opening/2))


    #---------------------------------------------------------------------------
    # Expansion front left

    channelFrontLeft1 = crocus.addLeft(baseLeft1, "channelFrontLeft", [
        mesh.Point(-expansionRadius/sqrt2, -expansionRadius/sqrt2, -coreBaseHeight),
        mesh.Point(-expansionRadius*np.cos(opening), -expansionRadius*np.sin(opening), -coreBaseHeight),
        mesh.Point(-expansionRadius/sqrt2, -expansionRadius/sqrt2, 0),
        mesh.Point(-expansionRadius*np.cos(opening), -expansionRadius*np.sin(opening), 0),
    ], nx=channelBackRight1.nx)

    expansionBaseFrontLeft1 = crocus.addLeft(channelFrontLeft1, "expansionBaseFrontLeft", [
        mesh.Point(-vesselHalfLength, -vesselHalfLength, -coreBaseHeight),
        mesh.Point(-vesselHalfLength, -expansionRadius*np.sin(opening), -coreBaseHeight),
        mesh.Point(-vesselHalfLength, -vesselHalfLength, 0),
        mesh.Point(-vesselHalfLength, -expansionRadius*np.sin(opening), 0),
    ], nx=expansionBaseBackRight1.nx)

    expansionBaseFrontLeft2 = crocus.addBack(expansionBaseFrontLeft1, "expansionBaseFrontLeft", [
        mesh.Point(-expansionRadius*np.cos(opening/2), -expansionRadius*np.sin(opening/2), -coreBaseHeight),
        mesh.Point(-vesselHalfLength, -expansionRadius*np.sin(opening/2), -coreBaseHeight),
        mesh.Point(-expansionRadius*np.cos(opening/2), -expansionRadius*np.sin(opening/2), 0),
        mesh.Point(-vesselHalfLength, -expansionRadius*np.sin(opening/2), 0),
    ], ny=expansionBaseBackRight2.ny)


    channelFrontLeft2 = crocus.addFront(baseFront2, "channelFrontLeft", [
        channelFrontLeft1.points[0],
        mesh.Point(-expansionRadius*np.sin(opening), -expansionRadius*np.cos(opening), -coreBaseHeight),
        channelFrontLeft1.points[4],
        mesh.Point(-expansionRadius*np.sin(opening), -expansionRadius*np.cos(opening), 0),
    ], ny=channelFrontLeft1.nx)

    expansionBaseFrontLeft3 = crocus.addFront(channelFrontLeft2, "expansionBaseFrontLeft", [
        expansionBaseFrontLeft1.points[0],
        mesh.Point(-expansionRadius*np.sin(opening), -vesselHalfLength, -coreBaseHeight),
        expansionBaseFrontLeft1.points[4],
        mesh.Point(-expansionRadius*np.sin(opening), -vesselHalfLength, 0)
    ], ny=expansionBaseFrontLeft1.nx)

    expansionBaseFrontLeft4 = crocus.addRight(expansionBaseFrontLeft3, "expansionBaseFrontLeft", [
        mesh.Point(-expansionRadius*np.sin(opening/2), -vesselHalfLength, -coreBaseHeight),
        mesh.Point(-expansionRadius*np.sin(opening/2), -expansionRadius*np.cos(opening/2), -coreBaseHeight),
        mesh.Point(-expansionRadius*np.sin(opening/2), -vesselHalfLength, 0),
        mesh.Point(-expansionRadius*np.sin(opening/2), -expansionRadius*np.cos(opening/2), 0),
    ], nx=expansionBaseFrontLeft2.ny)


    (
        expansionFrontLeft1, expansionFrontLeft2,
        expansionFrontLeft3, expansionFrontLeft4
    ) = crocus.extrudeTop([
        expansionBaseFrontLeft1, expansionBaseFrontLeft2,
        expansionBaseFrontLeft3, expansionBaseFrontLeft4
    ], "expansionFrontLeft", coreHeight, nz=expansionBackRight1.nz)

    expansionFrontLeft1.addEdge("arc", 1, 2, x=-expansionRadius*np.cos(opening/2), y=-expansionRadius*np.sin(opening/2))
    expansionFrontLeft1.addEdge("arc", 5, 6, x=-expansionRadius*np.cos(opening/2), y=-expansionRadius*np.sin(opening/2))
    expansionFrontLeft3.addEdge("arc", 2, 3, x=-expansionRadius*np.cos(opening/2), y=-expansionRadius*np.sin(opening/2))
    expansionFrontLeft3.addEdge("arc", 6, 7, x=-expansionRadius*np.cos(opening/2), y=-expansionRadius*np.sin(opening/2))


    #---------------------------------------------------------------------------
    # Expansion front right

    channelFrontRight1 = crocus.addRight(baseRight2, "channelFrontRight", [
        mesh.Point(expansionRadius/sqrt2, -expansionRadius/sqrt2, -coreBaseHeight),
        mesh.Point(expansionRadius*np.cos(opening), -expansionRadius*np.sin(opening), -coreBaseHeight),
        mesh.Point(expansionRadius/sqrt2, -expansionRadius/sqrt2, 0),
        mesh.Point(expansionRadius*np.cos(opening), -expansionRadius*np.sin(opening), 0),
    ], nx=channelBackRight1.nx)

    expansionBaseFrontRight1 = crocus.addRight(channelFrontRight1, "expansionBaseFrontRight", [
        mesh.Point(vesselHalfLength, -vesselHalfLength, -coreBaseHeight),
        mesh.Point(vesselHalfLength, -expansionRadius*np.sin(opening), -coreBaseHeight),
        mesh.Point(vesselHalfLength, -vesselHalfLength, 0),
        mesh.Point(vesselHalfLength, -expansionRadius*np.sin(opening), 0),
    ], nx=expansionBaseBackRight1.nx)

    expansionBaseFrontRight2 = crocus.addBack(expansionBaseFrontRight1, "expansionBaseFrontRight", [
        mesh.Point(vesselHalfLength, -expansionRadius*np.sin(opening/2), -coreBaseHeight),
        mesh.Point(expansionRadius*np.cos(opening/2), -expansionRadius*np.sin(opening/2), -coreBaseHeight),
        mesh.Point(vesselHalfLength, -expansionRadius*np.sin(opening/2), 0),
        mesh.Point(expansionRadius*np.cos(opening/2), -expansionRadius*np.sin(opening/2), 0),
    ], ny=expansionBaseBackRight2.ny)


    channelFrontRight2 = crocus.addFront(baseFront1, "channelFrontRight", [
        mesh.Point(expansionRadius*np.sin(opening), -expansionRadius*np.cos(opening), -coreBaseHeight),
        channelFrontRight1.points[1],
        mesh.Point(expansionRadius*np.sin(opening), -expansionRadius*np.cos(opening), 0),
        channelFrontRight1.points[5],
    ], ny=channelFrontRight1.nx)

    expansionBaseFrontRight3 = crocus.addFront(channelFrontRight2, "expansionBaseFrontRight", [
        mesh.Point(expansionRadius*np.sin(opening), -vesselHalfLength, -coreBaseHeight),
        expansionBaseFrontRight1.points[1],
        mesh.Point(expansionRadius*np.sin(opening), -vesselHalfLength, 0),
        expansionBaseFrontRight1.points[5],
    ], ny=expansionBaseFrontRight1.nx)

    expansionBaseFrontRight4 = crocus.addLeft(expansionBaseFrontRight3, "expansionBaseFrontRight", [
        mesh.Point(expansionRadius*np.sin(opening/2), -vesselHalfLength, -coreBaseHeight),
        mesh.Point(expansionRadius*np.sin(opening/2), -expansionRadius*np.cos(opening/2), -coreBaseHeight),
        mesh.Point(expansionRadius*np.sin(opening/2), -vesselHalfLength, 0),
        mesh.Point(expansionRadius*np.sin(opening/2), -expansionRadius*np.cos(opening/2), 0),
    ], nx=expansionBaseFrontRight2.ny)


    (
        expansionFrontRight1, expansionFrontRight2,
        expansionFrontRight3, expansionFrontRight4
    ) = crocus.extrudeTop([
        expansionBaseFrontRight1, expansionBaseFrontRight2,
        expansionBaseFrontRight3, expansionBaseFrontRight4
    ], "expansionFrontRight", coreHeight, nz=expansionBackRight1.nz)

    expansionFrontRight1.addEdge("arc", 3, 0, x=expansionRadius*np.sin(opening/2), y=-expansionRadius*np.cos(opening/2))
    expansionFrontRight1.addEdge("arc", 7, 4, x=expansionRadius*np.sin(opening/2), y=-expansionRadius*np.cos(opening/2))
    expansionFrontRight3.addEdge("arc", 2, 3, x=expansionRadius*np.sin(opening/2), y=-expansionRadius*np.cos(opening/2))
    expansionFrontRight3.addEdge("arc", 6, 7, x=expansionRadius*np.sin(opening/2), y=-expansionRadius*np.cos(opening/2))


    #---------------------------------------------------------------------------
    # Patches

    wall = ffn.Face("wall", boundaryType="wall")
    wall.addSubFace(coreFront.frontFace())
    wall.addSubFace(coreFront1.frontFace())
    wall.addSubFace(coreFront2.frontFace())
    wall.addSubFace(coreRight.rightFace())
    wall.addSubFace(coreRight1.rightFace())
    wall.addSubFace(coreRight2.rightFace())
    wall.addSubFace(coreBack.backFace())
    wall.addSubFace(coreBack1.backFace())
    wall.addSubFace(coreBack2.backFace())
    wall.addSubFace(coreLeft.leftFace())
    wall.addSubFace(coreLeft1.leftFace())
    wall.addSubFace(coreLeft2.leftFace())
    wall.addSubFace(channelBackLeft1.topFace())
    wall.addSubFace(channelBackLeft2.topFace())
    wall.addSubFace(channelBackRight1.topFace())
    wall.addSubFace(channelBackRight2.topFace())
    wall.addSubFace(channelFrontLeft1.topFace())
    wall.addSubFace(channelFrontLeft2.topFace())
    wall.addSubFace(channelFrontRight1.topFace())
    wall.addSubFace(channelFrontRight2.topFace())
    for block in [
        coreCenter, coreCenterFront, coreFront,
        baseCenterLeft, baseCenterBack, baseCenterRight,
        baseBackLeft, baseBackRight, baseFrontLeft, baseFrontRight,
        baseFront1, baseFront2,
        baseLeft, baseLeft1, baseLeft2,
        baseRight, baseRight1, baseRight2,
        baseBack, baseBack1, baseBack2,
        channelBackRight1, channelBackRight2,
        channelBackLeft1, channelBackLeft2,
        channelFrontRight1, channelFrontRight2,
        channelFrontLeft1, channelFrontLeft2,
        expansionBaseBackRight1, expansionBaseBackRight2, expansionBaseBackRight3, expansionBaseBackRight4,
        expansionBaseBackLeft1, expansionBaseBackLeft2, expansionBaseBackLeft3, expansionBaseBackLeft4,
        expansionBaseFrontRight1, expansionBaseFrontRight2, expansionBaseFrontRight3, expansionBaseFrontRight4,
        expansionBaseFrontLeft1, expansionBaseFrontLeft2, expansionBaseFrontLeft3, expansionBaseFrontLeft4
    ]:
        wall.addSubFace(block.bottomFace())
    wall.addSubFace(baseCenterBack.frontFace())
    wall.addSubFace(baseCenterLeft.rightFace())
    wall.addSubFace(baseCenterRight.leftFace())
    wall.addSubFace(baseFrontLeft.rightFace())
    wall.addSubFace(baseFrontRight.leftFace())
    wall.addSubFace(baseBack.backFace())
    wall.addSubFace(baseLeft.leftFace())
    wall.addSubFace(baseFront1.leftFace())
    wall.addSubFace(baseFront2.rightFace())
    wall.addSubFace(baseRight.rightFace())
    wall.addSubFace(channelFrontLeft1.backFace())
    wall.addSubFace(channelFrontLeft2.rightFace())
    wall.addSubFace(channelBackLeft1.frontFace())
    wall.addSubFace(channelBackLeft2.rightFace())
    wall.addSubFace(channelBackRight1.frontFace())
    wall.addSubFace(channelBackRight2.leftFace())
    wall.addSubFace(channelFrontRight1.backFace())
    wall.addSubFace(channelFrontRight2.leftFace())
    for block in [
        expansionBaseFrontLeft1,
        expansionBaseFrontLeft2,
        expansionBaseFrontRight2,
        expansionBaseFrontRight4,
        expansionBaseBackLeft1,
        expansionBaseBackLeft2,
        expansionBaseBackRight2,
        expansionBaseBackRight4,
        expansionFrontLeft1,
        expansionFrontLeft2,
        expansionFrontRight1,
        expansionFrontRight2,
        expansionFrontRight4,
        expansionBackLeft1,
        expansionBackLeft2,
        expansionBackRight1,
        expansionBackRight2,
        expansionBackRight4,
    ]:
        wall.addSubFace(block.leftFace())
    for block in [
        expansionBaseFrontLeft2,
        expansionBaseFrontLeft4,
        expansionBaseFrontRight1,
        expansionBaseFrontRight2,
        expansionBaseBackLeft2,
        expansionBaseBackLeft4,
        expansionBaseBackRight1,
        expansionBaseBackRight2,
        expansionFrontRight1,
        expansionFrontRight2,
        expansionFrontLeft1,
        expansionFrontLeft2,
        expansionFrontLeft4,
        expansionBackLeft1,
        expansionBackLeft2,
        expansionBackLeft4,
        expansionBackRight1,
        expansionBackRight2,
    ]:
        wall.addSubFace(block.rightFace())
    for block in [
        expansionBaseFrontLeft2,
        expansionBaseFrontLeft4,
        expansionBaseFrontRight2,
        expansionBaseFrontRight4,
        expansionBaseBackLeft3,
        expansionBaseBackLeft4,
        expansionBaseBackRight3,
        expansionBaseBackRight4,
        expansionFrontLeft2,
        expansionFrontLeft3,
        expansionFrontLeft4,
        expansionFrontRight2,
        expansionFrontRight3,
        expansionFrontRight4,
        expansionBackLeft3,
        expansionBackLeft4,
        expansionBackRight3,
        expansionBackRight4,
    ]:
        wall.addSubFace(block.backFace())
    for block in [
        expansionBaseFrontLeft3,
        expansionBaseFrontLeft4,
        expansionBaseFrontRight3,
        expansionBaseFrontRight4,
        expansionBaseBackLeft2,
        expansionBaseBackLeft4,
        expansionBaseBackRight2,
        expansionBaseBackRight4,
        expansionFrontLeft3,
        expansionFrontLeft4,
        expansionFrontRight3,
        expansionFrontRight4,
        expansionBackLeft2,
        expansionBackLeft3,
        expansionBackLeft4,
        expansionBackRight2,
        expansionBackRight3,
        expansionBackRight4,
    ]:
        wall.addSubFace(block.frontFace())

    atmosphereCore = ffn.Face("atmosphereCore", boundaryType="patch")
    atmosphereCore.addSubFace(coreCenter.topFace())
    atmosphereCore.addSubFace(coreCenterLeft.topFace())
    atmosphereCore.addSubFace(coreCenterFront.topFace())
    atmosphereCore.addSubFace(coreCenterRight.topFace())
    atmosphereCore.addSubFace(coreCenterBack.topFace())
    atmosphereCore.addSubFace(coreLeft.topFace())
    atmosphereCore.addSubFace(coreLeft1.topFace())
    atmosphereCore.addSubFace(coreLeft2.topFace())
    atmosphereCore.addSubFace(coreRight.topFace())
    atmosphereCore.addSubFace(coreRight1.topFace())
    atmosphereCore.addSubFace(coreRight2.topFace())
    atmosphereCore.addSubFace(coreFront.topFace())
    atmosphereCore.addSubFace(coreFront1.topFace())
    atmosphereCore.addSubFace(coreFront2.topFace())
    atmosphereCore.addSubFace(coreFrontLeft.topFace())
    atmosphereCore.addSubFace(coreFrontRight.topFace())
    atmosphereCore.addSubFace(coreBack.topFace())
    atmosphereCore.addSubFace(coreBackLeft.topFace())
    atmosphereCore.addSubFace(coreBackRight.topFace())
    atmosphereCore.addSubFace(coreBack1.topFace())
    atmosphereCore.addSubFace(coreBack2.topFace())


    atmosphereExpansionFrontRight = ffn.Face("atmosphereExpansionFrontRight", boundaryType="patch")
    atmosphereExpansionFrontRight.addSubFace(expansionFrontRight1.topFace())
    atmosphereExpansionFrontRight.addSubFace(expansionFrontRight2.topFace())
    atmosphereExpansionFrontRight.addSubFace(expansionFrontRight3.topFace())
    atmosphereExpansionFrontRight.addSubFace(expansionFrontRight4.topFace())

    atmosphereExpansionBackRight = ffn.Face("atmosphereExpansionBackRight", boundaryType="patch")
    atmosphereExpansionBackRight.addSubFace(expansionBackRight1.topFace())
    atmosphereExpansionBackRight.addSubFace(expansionBackRight2.topFace())
    atmosphereExpansionBackRight.addSubFace(expansionBackRight3.topFace())
    atmosphereExpansionBackRight.addSubFace(expansionBackRight4.topFace())

    atmosphereExpansionBackLeft = ffn.Face("atmosphereExpansionBackLeft", boundaryType="patch")
    atmosphereExpansionBackLeft.addSubFace(expansionBackLeft1.topFace())
    atmosphereExpansionBackLeft.addSubFace(expansionBackLeft2.topFace())
    atmosphereExpansionBackLeft.addSubFace(expansionBackLeft3.topFace())
    atmosphereExpansionBackLeft.addSubFace(expansionBackLeft4.topFace())

    atmosphereExpansionFrontLeft = ffn.Face("atmosphereExpansionFrontLeft", boundaryType="patch")
    atmosphereExpansionFrontLeft.addSubFace(expansionFrontLeft1.topFace())
    atmosphereExpansionFrontLeft.addSubFace(expansionFrontLeft2.topFace())
    atmosphereExpansionFrontLeft.addSubFace(expansionFrontLeft3.topFace())
    atmosphereExpansionFrontLeft.addSubFace(expansionFrontLeft4.topFace())

    crocus.addBoundary(wall)
    crocus.addBoundary(atmosphereCore)
    crocus.addBoundary(atmosphereExpansionFrontRight)
    crocus.addBoundary(atmosphereExpansionBackRight)
    crocus.addBoundary(atmosphereExpansionBackLeft)
    crocus.addBoundary(atmosphereExpansionFrontLeft)

    return(crocus)



thMesh = thermalHydraulicMesh(region='fluidRegion')


#==============================================================================*
# Time Folder

timeFolder0 = ffn.TimeFolder(0)

isOpenValve = True

T = ffn.Field("T", region=thMesh.region)
T.dimensions = ffn.Dimension(default='T')
T.internalField = 300
T.set_boundary_condition("wall", bc.FixedValue(T.internalField))
T.set_boundary_condition('"atmosphere.*"', bc.ZeroGradient())


U = ffn.Field("U", region=thMesh.region)
U.dimensions = ffn.Dimension(default='U')
U.internalField = ffn.Vector(0, 0, 0)
U.set_boundary_condition("wall", bc.NoSlip())
U.set_boundary_condition('atmosphereCore', bc.PressureInletOutletVelocity(ffn.Vector(0, 0, 0)))
if (isOpenValve):
    U.set_boundary_condition('"atmosphereExpansion.*"', bc.PressureInletOutletVelocity(ffn.Vector(0, 0, 0)))
else:
    U.set_boundary_condition('"atmosphereExpansion.*"', bc.NoSlip())

p = ffn.Field("p", region=thMesh.region)
p.dimensions = ffn.Dimension(default='p')
p.internalField = 1e5
p.set_boundary_condition("wall", bc.FixedFluxPressure(p.internalField))
p.set_boundary_condition('atmosphereCore', bc.TotalPressure(p.internalField))
if (isOpenValve):
    p.set_boundary_condition('"atmosphereExpansion.*"', bc.TotalPressure(p.internalField))
else:
    p.set_boundary_condition('"atmosphereExpansion.*"', bc.FixedFluxPressure(p.internalField))

p_rgh = ffn.Field("p_rgh", region=thMesh.region)
p_rgh.dimensions = ffn.Dimension(default='p')
p_rgh.internalField = 1.1e5
p_rgh.set_boundary_condition("wall", bc.FixedFluxPressure(p_rgh.internalField))
p_rgh.set_boundary_condition('atmosphereCore', bc.TotalPressure(p_rgh.internalField))
if (isOpenValve):
    p_rgh.set_boundary_condition('"atmosphereExpansion.*"', bc.TotalPressure(p_rgh.internalField))
else:
    p_rgh.set_boundary_condition('"atmosphereExpansion.*"', bc.FixedFluxPressure(p_rgh.internalField))

alpha_water = ffn.Field("alpha.water", region=thMesh.region)
alpha_water.dimensions = ffn.Dimension()
alpha_water.internalField = 0
alpha_water.set_boundary_condition("wall", bc.ZeroGradient())
alpha_water.set_boundary_condition('atmosphereCore', bc.InletOutlet(inletValue=0, value=0))
if (isOpenValve):
    alpha_water.set_boundary_condition('"atmosphereExpansion.*"', bc.InletOutlet(inletValue=0, value=0))
else:
    alpha_water.set_boundary_condition('"atmosphereExpansion.*"', bc.ZeroGradient())

timeFolder0.append(T)
timeFolder0.append(U)
timeFolder0.append(p)
timeFolder0.append(p_rgh)
timeFolder0.append(alpha_water)


#==============================================================================*
# Solver

thSolver = ffn.CompressibleInterFoam(
    "fluidRegion",
    mesh=thMesh,
    removeBaffles=False,
    isSetFvSolutionToDefault=True,
    isSetFvSchemesToDefault=True
)

# thSolver.decomposeParDict.numberOfSubdomains = 8
# thSolver.decomposeParDict.method = 'scotch'

thSolver.turbulenceProperties.simulationType = 'laminar'

multiPhase = ffn.thermophysicalProperty.MultiPhase()

multiPhase.appendPhase(ffn.thermophysicalProperty.WaterPerfectFluid(ext='water'))
multiPhase.appendPhase(ffn.thermophysicalProperty.AirPerfectGas(ext='air'))

thSolver.thermophysicalProperties = multiPhase

thSolver.setFieldsDict.add_default_values(alpha_water, 0)

thSolver.setFieldsDict.add_cylinder_to_cell(
    fieldValues=[(alpha_water, 1)],
    radius=0.75,
    point2=ffn.Vector(0, 0, 0.98)
)
thSolver.setFieldsDict.add_box_to_cell(
    fieldValues=[(alpha_water, 1)],
    lowCorner=ffn.Vector(-1, -1, -1),
    highCorner=ffn.Vector(1, 1, 0.1)
)


#==============================================================================*
# Model

model = ffn.Model()

model.settings.application = "GeN-Foam"
model.settings.endTime = 10
model.settings.deltaT = 0.0001
model.settings.writeControl = "adjustableRunTime"
model.settings.writeInterval = 0.1
model.settings.adjustTimeStep = True
model.settings.runTimeModifiable = True
model.settings.maxCo = 0.1
model.settings.maxDeltaT = 0.1
model.settings.maxAlphaCo = 0.1

model.solvers.append(thSolver)
model.add_time_folder(timeFolder0)

print(model)

#==============================================================================*
# Run

if (True):
    ffn.allclean()
    model.export_to_openfoam()

    # ffn.run_preprocessing(model=model)
    ffn.run(model=model, is_preprocessing=True)

#==============================================================================*
# Post-processing

model.plot_mesh(region=thMesh.region, show_edges=True)
model.plot_mesh(region=thMesh.region, show_edges=True, normal='z')
model.plot_boundary(region=thMesh.region, boundaryName="wall")
model.plot_boundary(region=thMesh.region, boundaryName="atmosphereCore")
model.plot_boundary(region=thMesh.region, boundaryName='atmosphereExpansionFrontRight')
model.plot_boundary(region=thMesh.region, boundaryName='atmosphereExpansionFrontLeft')
model.plot_boundary(region=thMesh.region, boundaryName='atmosphereExpansionBackRight')
model.plot_boundary(region=thMesh.region, boundaryName='atmosphereExpansionBackLeft')

model.plot_animation(
    region=thMesh.region,
    fieldName="U",
    cmap='Blues',
    threshold=[0.9, 1.1],
    thresholdFieldName="alpha.water",
    unit='m/s'
)

#==============================================================================*
