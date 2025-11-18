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

    coreCenter = crocus.create_cube("core", -coreRadius/6, -coreRadius/6, 0, coreRadius/6, coreRadius/6, coreHeight, 3, 3, nzCore)

    coreCenterLeft = crocus.extrude_left([coreCenter], "core", coreRadius/3, nx=coreCenter.nx)
    coreCenterRight = crocus.extrude_right([coreCenter], "core", coreRadius/3, nx=coreCenter.nx)

    coreFrontLeft, coreCenterFront, coreFrontRight = crocus.extrude_front(
        [coreCenterLeft, coreCenter, coreCenterRight], "core", coreRadius/3, ny=coreCenter.ny
    )

    coreBackLeft, coreCenterBack, coreBackRight = crocus.extrude_back(
        [coreCenterLeft, coreCenter, coreCenterRight], "core", coreRadius/3, ny=coreCenter.ny
    )


    coreFront = crocus.add_front(coreCenterFront, "coreFront", [
        mesh.Point(-coreRadius*np.sin(opening), -coreRadius*np.cos(opening), 0),
        mesh.Point(+coreRadius*np.sin(opening), -coreRadius*np.cos(opening), 0),
        mesh.Point(-coreRadius*np.sin(opening), -coreRadius*np.cos(opening), coreHeight),
        mesh.Point(+coreRadius*np.sin(opening), -coreRadius*np.cos(opening), coreHeight),
    ], ny=4)

    coreFront.addEdge("arc", 0, 1, x=0, y=-coreRadius)
    coreFront.addEdge("arc", 4, 5, x=0, y=-coreRadius)

    ffn.PointBottomRight0 = mesh.Point(coreRadius/sqrt2, -coreRadius/sqrt2, 0)
    ffn.PointBottomRightH = mesh.Point(coreRadius/sqrt2, -coreRadius/sqrt2, coreHeight)

    coreFront1 = crocus.add_right(coreFront, "coreFront", [
        ffn.PointBottomRight0,
        coreFrontRight.points[1],
        ffn.PointBottomRightH,
        coreFrontRight.points[5],
    ], nx=coreFrontRight.nx)

    coreFront1.addEdge("arc", 0, 1, x=coreRadius*np.sin(2*opening), y=-coreRadius*np.cos(2*opening))
    coreFront1.addEdge("arc", 4, 5, x=coreRadius*np.sin(2*opening), y=-coreRadius*np.cos(2*opening))

    ffn.PointBottomLeft0 = mesh.Point(-coreRadius/sqrt2, -coreRadius/sqrt2, 0)
    ffn.PointBottomLeftH = mesh.Point(-coreRadius/sqrt2, -coreRadius/sqrt2, coreHeight)

    coreFront2 = crocus.add_left(coreFront, "coreFront", [
        ffn.PointBottomLeft0,
        coreFrontLeft.points[0],
        ffn.PointBottomLeftH,
        coreFrontLeft.points[4],
    ], nx=coreFrontLeft.nx)

    coreFront2.addEdge("arc", 0, 1, x=-coreRadius*np.sin(2*opening), y=-coreRadius*np.cos(2*opening))
    coreFront2.addEdge("arc", 4, 5, x=-coreRadius*np.sin(2*opening), y=-coreRadius*np.cos(2*opening))

    coreRight = crocus.add_right(coreCenterRight, "coreRight", [
        mesh.Point(+coreRadius*np.cos(opening), -coreRadius*np.sin(opening), 0),
        mesh.Point(+coreRadius*np.cos(opening), +coreRadius*np.sin(opening), 0),
        mesh.Point(+coreRadius*np.cos(opening), -coreRadius*np.sin(opening), coreHeight),
        mesh.Point(+coreRadius*np.cos(opening), +coreRadius*np.sin(opening), coreHeight),
    ], nx=coreFront.ny)

    coreRight.addEdge("arc", 1, 2, x=coreRadius, y=0)
    coreRight.addEdge("arc", 5, 6, x=coreRadius, y=0)

    ffn.PointTopRight0 = mesh.Point(coreRadius/sqrt2, coreRadius/sqrt2, 0)
    ffn.PointTopRightH = mesh.Point(coreRadius/sqrt2, coreRadius/sqrt2, coreHeight)

    coreRight1 = crocus.add_back(coreRight, "coreRight", [
        ffn.PointTopRight0,
        coreBackRight.points[2],
        ffn.PointTopRightH,
        coreBackRight.points[6],
    ], ny=coreBackRight.ny)

    coreRight1.addEdge("arc", 1, 2, x=coreRadius*np.cos(2*opening), y=coreRadius*np.sin(2*opening))
    coreRight1.addEdge("arc", 5, 6, x=coreRadius*np.cos(2*opening), y=coreRadius*np.sin(2*opening))

    coreRight2 = crocus.add_front(coreRight, "coreRight", [
        coreFrontRight.points[1],
        ffn.PointBottomRight0,
        coreFrontRight.points[5],
        ffn.PointBottomRightH,
    ], ny=coreFrontRight.ny)

    coreRight2.addEdge("arc", 1, 2, x=coreRadius*np.cos(2*opening), y=-coreRadius*np.sin(2*opening))
    coreRight2.addEdge("arc", 5, 6, x=coreRadius*np.cos(2*opening), y=-coreRadius*np.sin(2*opening))

    coreLeft = crocus.add_left(coreCenterLeft, "coreLeft", [
        mesh.Point(-coreRadius*np.cos(opening), -coreRadius*np.sin(opening), 0),
        mesh.Point(-coreRadius*np.cos(opening), +coreRadius*np.sin(opening), 0),
        mesh.Point(-coreRadius*np.cos(opening), -coreRadius*np.sin(opening), coreHeight),
        mesh.Point(-coreRadius*np.cos(opening), +coreRadius*np.sin(opening), coreHeight),
    ], nx=coreFront.ny)

    coreLeft.addEdge("arc", 3, 0, x=-coreRadius, y=0)
    coreLeft.addEdge("arc", 7, 4, x=-coreRadius, y=0)

    coreLeft1 = crocus.add_front(coreLeft, "coreLeft", [
        ffn.PointBottomLeft0,
        coreFrontLeft.points[0],
        ffn.PointBottomLeftH,
        coreFrontLeft.points[4],
    ], ny=coreFrontLeft.ny)

    coreLeft1.addEdge("arc", 3, 0, x=-coreRadius*np.cos(2*opening), y=-coreRadius*np.sin(2*opening))
    coreLeft1.addEdge("arc", 7, 4, x=-coreRadius*np.cos(2*opening), y=-coreRadius*np.sin(2*opening))

    ffn.PointTopLeft0 = mesh.Point(-coreRadius/sqrt2, coreRadius/sqrt2, 0)
    ffn.PointTopLeftH = mesh.Point(-coreRadius/sqrt2, coreRadius/sqrt2, coreHeight)

    coreLeft2 = crocus.add_back(coreLeft, "coreLeft", [
        coreBackLeft.points[3],
        ffn.PointTopLeft0,
        coreBackLeft.points[7],
        ffn.PointTopLeftH,
    ], ny=coreFrontRight.ny)

    coreLeft2.addEdge("arc", 3, 0, x=-coreRadius*np.cos(2*opening), y=coreRadius*np.sin(2*opening))
    coreLeft2.addEdge("arc", 7, 4, x=-coreRadius*np.cos(2*opening), y=coreRadius*np.sin(2*opening))

    coreBack = crocus.add_back(coreCenterBack, "coreBack", [
        mesh.Point(+coreRadius*np.sin(opening), +coreRadius*np.cos(opening), 0),
        mesh.Point(-coreRadius*np.sin(opening), +coreRadius*np.cos(opening), 0),
        mesh.Point(+coreRadius*np.sin(opening), +coreRadius*np.cos(opening), coreHeight),
        mesh.Point(-coreRadius*np.sin(opening), +coreRadius*np.cos(opening), coreHeight),
    ], ny=coreFront.ny)

    coreBack.addEdge("arc", 2, 3, x=0, y=coreRadius)
    coreBack.addEdge("arc", 6, 7, x=0, y=coreRadius)

    coreBack1 = crocus.add_left(coreBack, "coreBack", [
        coreBackLeft.points[3],
        ffn.PointTopLeft0,
        coreBackLeft.points[7],
        ffn.PointTopLeftH,
    ], nx=coreFrontLeft.nx)

    coreBack1.addEdge("arc", 2, 3, x=-coreRadius*np.sin(2*opening), y=coreRadius*np.cos(2*opening))
    coreBack1.addEdge("arc", 6, 7, x=-coreRadius*np.sin(2*opening), y=coreRadius*np.cos(2*opening))

    coreBack2 = crocus.add_right(coreBack, "coreBack", [
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
    ) = crocus.extrude_bottom([
        coreCenterLeft, coreCenterBack, coreCenterRight,
        coreBackLeft, coreBackRight, coreFrontLeft, coreFrontRight,
        coreFront1, coreFront2,
        coreLeft, coreLeft1, coreLeft2,
        coreRight, coreRight1, coreRight2,
        coreBack, coreBack1, coreBack2,
    ], "coreBase", coreBaseHeight, nz=5)

    #---------------------------------------------------------------------------
    # Expansion back right

    channelBackRight1 = crocus.add_right(baseRight1, "channelBackRight", [
        mesh.Point(expansionRadius*np.cos(opening), expansionRadius*np.sin(opening), -coreBaseHeight),
        mesh.Point(expansionRadius/sqrt2, expansionRadius/sqrt2, -coreBaseHeight),
        mesh.Point(expansionRadius*np.cos(opening), expansionRadius*np.sin(opening), 0),
        mesh.Point(expansionRadius/sqrt2, expansionRadius/sqrt2, 0),
    ], nx=2)

    expansionBaseBackRight1 = crocus.add_right(channelBackRight1, "expansionBaseBackRight", [
        mesh.Point(vesselHalfLength, expansionRadius*np.sin(opening), -coreBaseHeight),
        mesh.Point(vesselHalfLength, vesselHalfLength, -coreBaseHeight),
        mesh.Point(vesselHalfLength, expansionRadius*np.sin(opening), 0),
        mesh.Point(vesselHalfLength, vesselHalfLength, 0),
    ], nx=3)

    expansionBaseBackRight2 = crocus.add_front(expansionBaseBackRight1, "expansionBaseBackRight", [
        mesh.Point(expansionRadius*np.cos(opening/2), expansionRadius*np.sin(opening/2), -coreBaseHeight),
        mesh.Point(vesselHalfLength, expansionRadius*np.sin(opening/2), -coreBaseHeight),
        mesh.Point(expansionRadius*np.cos(opening/2), expansionRadius*np.sin(opening/2), 0),
        mesh.Point(vesselHalfLength, expansionRadius*np.sin(opening/2), 0),
    ], ny=4)


    channelBackRight2 = crocus.add_back(baseBack2, "channelBackRight", [
        channelBackRight1.points[2],
        mesh.Point(expansionRadius*np.sin(opening), expansionRadius*np.cos(opening), -coreBaseHeight),
        channelBackRight1.points[6],
        mesh.Point(expansionRadius*np.sin(opening), expansionRadius*np.cos(opening), 0),
    ], ny=channelBackRight1.nx)

    expansionBaseBackRight3 = crocus.add_back(channelBackRight2, "expansionBaseBackRight", [
        expansionBaseBackRight1.points[2],
        mesh.Point(expansionRadius*np.sin(opening), vesselHalfLength, -coreBaseHeight),
        expansionBaseBackRight1.points[6],
        mesh.Point(expansionRadius*np.sin(opening), vesselHalfLength, 0)
    ], ny=expansionBaseBackRight1.nx)

    expansionBaseBackRight4 = crocus.add_left(expansionBaseBackRight3, "expansionBaseBackRight", [
        mesh.Point(expansionRadius*np.sin(opening/2), expansionRadius*np.cos(opening/2), -coreBaseHeight),
        mesh.Point(expansionRadius*np.sin(opening/2), vesselHalfLength, -coreBaseHeight),
        mesh.Point(expansionRadius*np.sin(opening/2), expansionRadius*np.cos(opening/2), 0),
        mesh.Point(expansionRadius*np.sin(opening/2), vesselHalfLength, 0),
    ], nx=expansionBaseBackRight2.ny)


    (
        expansionBackRight1, expansionBackRight2,
        expansionBackRight3, expansionBackRight4
    ) = crocus.extrude_top([
        expansionBaseBackRight1, expansionBaseBackRight2,
        expansionBaseBackRight3, expansionBaseBackRight4
    ], "expansionBackRight", coreHeight, nz=nzExpansionTanks)

    expansionBackRight1.addEdge("arc", 3, 0, x=expansionRadius*np.cos(opening/2), y=expansionRadius*np.sin(opening/2))
    expansionBackRight1.addEdge("arc", 7, 4, x=expansionRadius*np.cos(opening/2), y=expansionRadius*np.sin(opening/2))
    expansionBackRight3.addEdge("arc", 0, 1, x=expansionRadius*np.cos(opening/2), y=expansionRadius*np.sin(opening/2))
    expansionBackRight3.addEdge("arc", 4, 5, x=expansionRadius*np.cos(opening/2), y=expansionRadius*np.sin(opening/2))


    #---------------------------------------------------------------------------
    # Expansion back left


    channelBackLeft1 = crocus.add_left(baseLeft2, "channelBackLeft", [
        mesh.Point(-expansionRadius*np.cos(opening), expansionRadius*np.sin(opening), -coreBaseHeight),
        mesh.Point(-expansionRadius/sqrt2, expansionRadius/sqrt2, -coreBaseHeight),
        mesh.Point(-expansionRadius*np.cos(opening), expansionRadius*np.sin(opening), 0),
        mesh.Point(-expansionRadius/sqrt2, expansionRadius/sqrt2, 0),
    ], nx=channelBackRight1.nx)

    expansionBaseBackLeft1 = crocus.add_left(channelBackLeft1, "expansionBaseBackLeft", [
        mesh.Point(-vesselHalfLength, expansionRadius*np.sin(opening), -coreBaseHeight),
        mesh.Point(-vesselHalfLength, vesselHalfLength, -coreBaseHeight),
        mesh.Point(-vesselHalfLength, expansionRadius*np.sin(opening), 0),
        mesh.Point(-vesselHalfLength, vesselHalfLength, 0),
    ], nx=expansionBaseBackRight1.nx)

    expansionBaseBackLeft2 = crocus.add_front(expansionBaseBackLeft1, "expansionBaseBackLeft", [
        mesh.Point(-vesselHalfLength, expansionRadius*np.sin(opening/2), -coreBaseHeight),
        mesh.Point(-expansionRadius*np.cos(opening/2), expansionRadius*np.sin(opening/2), -coreBaseHeight),
        mesh.Point(-vesselHalfLength, expansionRadius*np.sin(opening/2), 0),
        mesh.Point(-expansionRadius*np.cos(opening/2), expansionRadius*np.sin(opening/2), 0),
    ], ny=expansionBaseBackRight2.ny)


    channelBackLeft2 = crocus.add_back(baseBack1, "channelBackLeft", [
        mesh.Point(-expansionRadius*np.sin(opening), expansionRadius*np.cos(opening), -coreBaseHeight),
        channelBackLeft1.points[3],
        mesh.Point(-expansionRadius*np.sin(opening), expansionRadius*np.cos(opening), 0),
        channelBackLeft1.points[7],
    ], ny=channelBackLeft1.nx)

    expansionBaseBackLeft3 = crocus.add_back(channelBackLeft2, "expansionBaseBackLeft", [
        mesh.Point(-expansionRadius*np.sin(opening), vesselHalfLength, -coreBaseHeight),
        expansionBaseBackLeft1.points[3],
        mesh.Point(-expansionRadius*np.sin(opening), vesselHalfLength, 0),
        expansionBaseBackLeft1.points[7],
    ], ny=expansionBaseBackLeft1.nx)

    expansionBaseBackLeft4 = crocus.add_right(expansionBaseBackLeft3, "expansionBaseBackLeft", [
        mesh.Point(-expansionRadius*np.sin(opening/2), expansionRadius*np.cos(opening/2), -coreBaseHeight),
        mesh.Point(-expansionRadius*np.sin(opening/2), vesselHalfLength, -coreBaseHeight),
        mesh.Point(-expansionRadius*np.sin(opening/2), expansionRadius*np.cos(opening/2), 0),
        mesh.Point(-expansionRadius*np.sin(opening/2), vesselHalfLength, 0),
    ], nx=expansionBaseBackLeft2.ny)



    (
        expansionBackLeft1, expansionBackLeft2,
        expansionBackLeft3, expansionBackLeft4
    ) = crocus.extrude_top([
        expansionBaseBackLeft1, expansionBaseBackLeft2,
        expansionBaseBackLeft3, expansionBaseBackLeft4
    ], "expansionBackLeft", coreHeight, nz=expansionBackRight1.nz)

    expansionBackLeft1.addEdge("arc", 1, 2, x=-expansionRadius*np.sin(opening/2), y=expansionRadius*np.cos(opening/2))
    expansionBackLeft1.addEdge("arc", 5, 6, x=-expansionRadius*np.sin(opening/2), y=expansionRadius*np.cos(opening/2))
    expansionBackLeft3.addEdge("arc", 0, 1, x=-expansionRadius*np.sin(opening/2), y=expansionRadius*np.cos(opening/2))
    expansionBackLeft3.addEdge("arc", 4, 5, x=-expansionRadius*np.sin(opening/2), y=expansionRadius*np.cos(opening/2))


    #---------------------------------------------------------------------------
    # Expansion front left

    channelFrontLeft1 = crocus.add_left(baseLeft1, "channelFrontLeft", [
        mesh.Point(-expansionRadius/sqrt2, -expansionRadius/sqrt2, -coreBaseHeight),
        mesh.Point(-expansionRadius*np.cos(opening), -expansionRadius*np.sin(opening), -coreBaseHeight),
        mesh.Point(-expansionRadius/sqrt2, -expansionRadius/sqrt2, 0),
        mesh.Point(-expansionRadius*np.cos(opening), -expansionRadius*np.sin(opening), 0),
    ], nx=channelBackRight1.nx)

    expansionBaseFrontLeft1 = crocus.add_left(channelFrontLeft1, "expansionBaseFrontLeft", [
        mesh.Point(-vesselHalfLength, -vesselHalfLength, -coreBaseHeight),
        mesh.Point(-vesselHalfLength, -expansionRadius*np.sin(opening), -coreBaseHeight),
        mesh.Point(-vesselHalfLength, -vesselHalfLength, 0),
        mesh.Point(-vesselHalfLength, -expansionRadius*np.sin(opening), 0),
    ], nx=expansionBaseBackRight1.nx)

    expansionBaseFrontLeft2 = crocus.add_back(expansionBaseFrontLeft1, "expansionBaseFrontLeft", [
        mesh.Point(-expansionRadius*np.cos(opening/2), -expansionRadius*np.sin(opening/2), -coreBaseHeight),
        mesh.Point(-vesselHalfLength, -expansionRadius*np.sin(opening/2), -coreBaseHeight),
        mesh.Point(-expansionRadius*np.cos(opening/2), -expansionRadius*np.sin(opening/2), 0),
        mesh.Point(-vesselHalfLength, -expansionRadius*np.sin(opening/2), 0),
    ], ny=expansionBaseBackRight2.ny)


    channelFrontLeft2 = crocus.add_front(baseFront2, "channelFrontLeft", [
        channelFrontLeft1.points[0],
        mesh.Point(-expansionRadius*np.sin(opening), -expansionRadius*np.cos(opening), -coreBaseHeight),
        channelFrontLeft1.points[4],
        mesh.Point(-expansionRadius*np.sin(opening), -expansionRadius*np.cos(opening), 0),
    ], ny=channelFrontLeft1.nx)

    expansionBaseFrontLeft3 = crocus.add_front(channelFrontLeft2, "expansionBaseFrontLeft", [
        expansionBaseFrontLeft1.points[0],
        mesh.Point(-expansionRadius*np.sin(opening), -vesselHalfLength, -coreBaseHeight),
        expansionBaseFrontLeft1.points[4],
        mesh.Point(-expansionRadius*np.sin(opening), -vesselHalfLength, 0)
    ], ny=expansionBaseFrontLeft1.nx)

    expansionBaseFrontLeft4 = crocus.add_right(expansionBaseFrontLeft3, "expansionBaseFrontLeft", [
        mesh.Point(-expansionRadius*np.sin(opening/2), -vesselHalfLength, -coreBaseHeight),
        mesh.Point(-expansionRadius*np.sin(opening/2), -expansionRadius*np.cos(opening/2), -coreBaseHeight),
        mesh.Point(-expansionRadius*np.sin(opening/2), -vesselHalfLength, 0),
        mesh.Point(-expansionRadius*np.sin(opening/2), -expansionRadius*np.cos(opening/2), 0),
    ], nx=expansionBaseFrontLeft2.ny)


    (
        expansionFrontLeft1, expansionFrontLeft2,
        expansionFrontLeft3, expansionFrontLeft4
    ) = crocus.extrude_top([
        expansionBaseFrontLeft1, expansionBaseFrontLeft2,
        expansionBaseFrontLeft3, expansionBaseFrontLeft4
    ], "expansionFrontLeft", coreHeight, nz=expansionBackRight1.nz)

    expansionFrontLeft1.addEdge("arc", 1, 2, x=-expansionRadius*np.cos(opening/2), y=-expansionRadius*np.sin(opening/2))
    expansionFrontLeft1.addEdge("arc", 5, 6, x=-expansionRadius*np.cos(opening/2), y=-expansionRadius*np.sin(opening/2))
    expansionFrontLeft3.addEdge("arc", 2, 3, x=-expansionRadius*np.cos(opening/2), y=-expansionRadius*np.sin(opening/2))
    expansionFrontLeft3.addEdge("arc", 6, 7, x=-expansionRadius*np.cos(opening/2), y=-expansionRadius*np.sin(opening/2))


    #---------------------------------------------------------------------------
    # Expansion front right

    channelFrontRight1 = crocus.add_right(baseRight2, "channelFrontRight", [
        mesh.Point(expansionRadius/sqrt2, -expansionRadius/sqrt2, -coreBaseHeight),
        mesh.Point(expansionRadius*np.cos(opening), -expansionRadius*np.sin(opening), -coreBaseHeight),
        mesh.Point(expansionRadius/sqrt2, -expansionRadius/sqrt2, 0),
        mesh.Point(expansionRadius*np.cos(opening), -expansionRadius*np.sin(opening), 0),
    ], nx=channelBackRight1.nx)

    expansionBaseFrontRight1 = crocus.add_right(channelFrontRight1, "expansionBaseFrontRight", [
        mesh.Point(vesselHalfLength, -vesselHalfLength, -coreBaseHeight),
        mesh.Point(vesselHalfLength, -expansionRadius*np.sin(opening), -coreBaseHeight),
        mesh.Point(vesselHalfLength, -vesselHalfLength, 0),
        mesh.Point(vesselHalfLength, -expansionRadius*np.sin(opening), 0),
    ], nx=expansionBaseBackRight1.nx)

    expansionBaseFrontRight2 = crocus.add_back(expansionBaseFrontRight1, "expansionBaseFrontRight", [
        mesh.Point(vesselHalfLength, -expansionRadius*np.sin(opening/2), -coreBaseHeight),
        mesh.Point(expansionRadius*np.cos(opening/2), -expansionRadius*np.sin(opening/2), -coreBaseHeight),
        mesh.Point(vesselHalfLength, -expansionRadius*np.sin(opening/2), 0),
        mesh.Point(expansionRadius*np.cos(opening/2), -expansionRadius*np.sin(opening/2), 0),
    ], ny=expansionBaseBackRight2.ny)


    channelFrontRight2 = crocus.add_front(baseFront1, "channelFrontRight", [
        mesh.Point(expansionRadius*np.sin(opening), -expansionRadius*np.cos(opening), -coreBaseHeight),
        channelFrontRight1.points[1],
        mesh.Point(expansionRadius*np.sin(opening), -expansionRadius*np.cos(opening), 0),
        channelFrontRight1.points[5],
    ], ny=channelFrontRight1.nx)

    expansionBaseFrontRight3 = crocus.add_front(channelFrontRight2, "expansionBaseFrontRight", [
        mesh.Point(expansionRadius*np.sin(opening), -vesselHalfLength, -coreBaseHeight),
        expansionBaseFrontRight1.points[1],
        mesh.Point(expansionRadius*np.sin(opening), -vesselHalfLength, 0),
        expansionBaseFrontRight1.points[5],
    ], ny=expansionBaseFrontRight1.nx)

    expansionBaseFrontRight4 = crocus.add_left(expansionBaseFrontRight3, "expansionBaseFrontRight", [
        mesh.Point(expansionRadius*np.sin(opening/2), -vesselHalfLength, -coreBaseHeight),
        mesh.Point(expansionRadius*np.sin(opening/2), -expansionRadius*np.cos(opening/2), -coreBaseHeight),
        mesh.Point(expansionRadius*np.sin(opening/2), -vesselHalfLength, 0),
        mesh.Point(expansionRadius*np.sin(opening/2), -expansionRadius*np.cos(opening/2), 0),
    ], nx=expansionBaseFrontRight2.ny)


    (
        expansionFrontRight1, expansionFrontRight2,
        expansionFrontRight3, expansionFrontRight4
    ) = crocus.extrude_top([
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
    wall.add_sub_face(coreFront.frontFace())
    wall.add_sub_face(coreFront1.frontFace())
    wall.add_sub_face(coreFront2.frontFace())
    wall.add_sub_face(coreRight.rightFace())
    wall.add_sub_face(coreRight1.rightFace())
    wall.add_sub_face(coreRight2.rightFace())
    wall.add_sub_face(coreBack.backFace())
    wall.add_sub_face(coreBack1.backFace())
    wall.add_sub_face(coreBack2.backFace())
    wall.add_sub_face(coreLeft.leftFace())
    wall.add_sub_face(coreLeft1.leftFace())
    wall.add_sub_face(coreLeft2.leftFace())
    wall.add_sub_face(channelBackLeft1.topFace())
    wall.add_sub_face(channelBackLeft2.topFace())
    wall.add_sub_face(channelBackRight1.topFace())
    wall.add_sub_face(channelBackRight2.topFace())
    wall.add_sub_face(channelFrontLeft1.topFace())
    wall.add_sub_face(channelFrontLeft2.topFace())
    wall.add_sub_face(channelFrontRight1.topFace())
    wall.add_sub_face(channelFrontRight2.topFace())
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
        wall.add_sub_face(block.bottomFace())
    wall.add_sub_face(baseCenterBack.frontFace())
    wall.add_sub_face(baseCenterLeft.rightFace())
    wall.add_sub_face(baseCenterRight.leftFace())
    wall.add_sub_face(baseFrontLeft.rightFace())
    wall.add_sub_face(baseFrontRight.leftFace())
    wall.add_sub_face(baseBack.backFace())
    wall.add_sub_face(baseLeft.leftFace())
    wall.add_sub_face(baseFront1.leftFace())
    wall.add_sub_face(baseFront2.rightFace())
    wall.add_sub_face(baseRight.rightFace())
    wall.add_sub_face(channelFrontLeft1.backFace())
    wall.add_sub_face(channelFrontLeft2.rightFace())
    wall.add_sub_face(channelBackLeft1.frontFace())
    wall.add_sub_face(channelBackLeft2.rightFace())
    wall.add_sub_face(channelBackRight1.frontFace())
    wall.add_sub_face(channelBackRight2.leftFace())
    wall.add_sub_face(channelFrontRight1.backFace())
    wall.add_sub_face(channelFrontRight2.leftFace())
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
        wall.add_sub_face(block.leftFace())
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
        wall.add_sub_face(block.rightFace())
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
        wall.add_sub_face(block.backFace())
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
        wall.add_sub_face(block.frontFace())

    atmosphereCore = ffn.Face("atmosphereCore", boundaryType="patch")
    atmosphereCore.add_sub_face(coreCenter.topFace())
    atmosphereCore.add_sub_face(coreCenterLeft.topFace())
    atmosphereCore.add_sub_face(coreCenterFront.topFace())
    atmosphereCore.add_sub_face(coreCenterRight.topFace())
    atmosphereCore.add_sub_face(coreCenterBack.topFace())
    atmosphereCore.add_sub_face(coreLeft.topFace())
    atmosphereCore.add_sub_face(coreLeft1.topFace())
    atmosphereCore.add_sub_face(coreLeft2.topFace())
    atmosphereCore.add_sub_face(coreRight.topFace())
    atmosphereCore.add_sub_face(coreRight1.topFace())
    atmosphereCore.add_sub_face(coreRight2.topFace())
    atmosphereCore.add_sub_face(coreFront.topFace())
    atmosphereCore.add_sub_face(coreFront1.topFace())
    atmosphereCore.add_sub_face(coreFront2.topFace())
    atmosphereCore.add_sub_face(coreFrontLeft.topFace())
    atmosphereCore.add_sub_face(coreFrontRight.topFace())
    atmosphereCore.add_sub_face(coreBack.topFace())
    atmosphereCore.add_sub_face(coreBackLeft.topFace())
    atmosphereCore.add_sub_face(coreBackRight.topFace())
    atmosphereCore.add_sub_face(coreBack1.topFace())
    atmosphereCore.add_sub_face(coreBack2.topFace())


    atmosphereExpansionFrontRight = ffn.Face("atmosphereExpansionFrontRight", boundaryType="patch")
    atmosphereExpansionFrontRight.add_sub_face(expansionFrontRight1.topFace())
    atmosphereExpansionFrontRight.add_sub_face(expansionFrontRight2.topFace())
    atmosphereExpansionFrontRight.add_sub_face(expansionFrontRight3.topFace())
    atmosphereExpansionFrontRight.add_sub_face(expansionFrontRight4.topFace())

    atmosphereExpansionBackRight = ffn.Face("atmosphereExpansionBackRight", boundaryType="patch")
    atmosphereExpansionBackRight.add_sub_face(expansionBackRight1.topFace())
    atmosphereExpansionBackRight.add_sub_face(expansionBackRight2.topFace())
    atmosphereExpansionBackRight.add_sub_face(expansionBackRight3.topFace())
    atmosphereExpansionBackRight.add_sub_face(expansionBackRight4.topFace())

    atmosphereExpansionBackLeft = ffn.Face("atmosphereExpansionBackLeft", boundaryType="patch")
    atmosphereExpansionBackLeft.add_sub_face(expansionBackLeft1.topFace())
    atmosphereExpansionBackLeft.add_sub_face(expansionBackLeft2.topFace())
    atmosphereExpansionBackLeft.add_sub_face(expansionBackLeft3.topFace())
    atmosphereExpansionBackLeft.add_sub_face(expansionBackLeft4.topFace())

    atmosphereExpansionFrontLeft = ffn.Face("atmosphereExpansionFrontLeft", boundaryType="patch")
    atmosphereExpansionFrontLeft.add_sub_face(expansionFrontLeft1.topFace())
    atmosphereExpansionFrontLeft.add_sub_face(expansionFrontLeft2.topFace())
    atmosphereExpansionFrontLeft.add_sub_face(expansionFrontLeft3.topFace())
    atmosphereExpansionFrontLeft.add_sub_face(expansionFrontLeft4.topFace())

    crocus.add_boundary(wall)
    crocus.add_boundary(atmosphereCore)
    crocus.add_boundary(atmosphereExpansionFrontRight)
    crocus.add_boundary(atmosphereExpansionBackRight)
    crocus.add_boundary(atmosphereExpansionBackLeft)
    crocus.add_boundary(atmosphereExpansionFrontLeft)

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
