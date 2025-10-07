#!/usr/bin/env python
"""
Script to generate an hexagonal mesh for reactor core. The script uses the
Salome API to generate UNV file and relies on a lattice map file in the standard
MCNP format.

-----------
HOW TO USE:
-----------
Go to the section "Beginning of MAIN for Users" and uses all the functions and
objects developed to ease the generation of meshes. The common approach is to
use the "LargeMesh" object. Then create a first sub-mesh using the
"createSubMesh()" method.


Author: Eymeric Simonnot, EPFL,   04/2023, First version
        Thomas Guilbaud, EPFL, 29/06/2023, Refactoring
                               16/11/2023, Refactoring, update to MCNP-like
                                 lattice format, improve usability
"""

#=============================================================================*
# Imports
#=============================================================================*

import sys
import salome
import numpy as np

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()

from SketchAPI import *
from salome.shaper import model

import SMESH, SALOMEDS
from salome.smesh import smeshBuilder
import SHAPERSTUDY


#=============================================================================*
# Useful Objects
#=============================================================================*

class SubMesh:
    """ Sub mesh class """

    def __init__(
        self,
        name: str,
        faceShape,
        latticeKey: str,
        latticeMapFilename: str,
        flatToFlat: float,
        zCenterPosition: float,
        height: float,
        nz: int=1,
        nSegments: int=1,
        isMergeLatticeElementsInOneMesh: bool=False,
        keyToNotInclude: list=[],
        isRemoveAllBaffles: bool=False
    ) -> None:
        self.name = name
        self.latticeKey = latticeKey
        self.faceShape = faceShape

        self.zBottomPosition = zCenterPosition - height/2
        self.zCenterPosition = zCenterPosition
        self.height = height
        self.nz = nz

        self.elementMesh = smesh.Mesh(faceShape, name)
        self.elementMesh = meshElement(self.elementMesh, height, nz, nSegments=nSegments)

        self.mesh = createMeshWithLattice(
            mesh=self.elementMesh,
            latticeMapFilename=latticeMapFilename,
            flatToFlat=flatToFlat,
            z=self.zBottomPosition,
            key=latticeKey,
            name=name,
            isMergeLatticeElementsInOneMesh=isMergeLatticeElementsInOneMesh,
            keyToNotInclude=keyToNotInclude,
            isRemoveAllBaffles=isRemoveAllBaffles
        )


class LargeMesh:
    """
    Main object that store and manipulate the final mesh, ready to be exported.
    Works only for hexagonal meshes.
    """
    def __init__(self, name: str, flatToFlat: float, faceShape) -> None:
        """
        Parameters
        ----------
        name: Name of the final mesh.
        flatToFlat: Flat-to-flat distance between the elements (also called picth).
        faceShape: Face object extracted from the Salome Shaper study.
        """
        self.name = name

        self.flatToFlat = flatToFlat
        self.summitToSummit = self.flatToFlat * 2/np.sqrt(3)
        self.side = self.summitToSummit/2

        self.faceShape = faceShape

        self.listSubMesh = []


    def createSubMesh(
        self,
        name: str,
        latticeMapFilename: str,
        height: float,
        nz: int=1,
        nSegments: int=1,
        zCenterPosition: float=0,
        latticeKey: str="F",
        isMergeLatticeElementsInOneMesh: bool=False,
        keyToNotInclude: list=[],
        topPatchName: str="",
        bottomPatchName: str="",
        isRemoveAllBaffles: bool=False
    ) -> None:
        """
        Create a sub-mesh.

        name: Name of the sub-mesh.
        latticeMapFilename: Lattice map filename.
        height: Height of the sub-mesh.
        nz: Number of axial cells.
        zCenterPosition:
        latticeKey: Filter by lattice keys (string or list of strings).
        isMergeLatticeElementsInOneMesh: If true, places a mesh on each non-0
            element in latticeMapFilename. Else, places a mesh on the specified key.
        keyToNotInclude: Lattice keys to not include. Can be used to reduce the number of baffles.
        topPatchName: Name top patch of the sub-mesh.
        bottomPatchName: Name bottom patch of the sub-mesh.
        isRemoveAllBaffles: If true, remove all the internal lateral faces.
        """
        # Check lattice key is present in lattice map
        with open(latticeMapFilename, 'r') as file:
            if (type(latticeKey) == str):
                if (latticeKey not in file.read()):
                    print(f"Lattice key '{latticeKey}' is not present in {latticeMapFilename} file. Mesh not generated")
                    return

        # Create sub mesh
        subMesh = SubMesh(
            name=name,
            faceShape=self.faceShape,
            latticeKey=latticeKey,
            latticeMapFilename=latticeMapFilename,
            flatToFlat=self.flatToFlat,
            zCenterPosition=zCenterPosition,
            height=height,
            nz=nz,
            nSegments=nSegments,
            isMergeLatticeElementsInOneMesh=isMergeLatticeElementsInOneMesh,
            keyToNotInclude=keyToNotInclude,
            isRemoveAllBaffles=isRemoveAllBaffles
        )

        # Replace or remove top patch
        if (topPatchName != ""):
            renameGroupByName(subMesh.mesh, topPatchName, "inlet")
        else:
            removeGroupByName(subMesh.mesh, "inlet")

        # Replace or remove bottom patch
        if (bottomPatchName != ""):
            renameGroupByName(subMesh.mesh, bottomPatchName, "outlet")
        else:
            removeGroupByName(subMesh.mesh, "outlet")

        self.listSubMesh.append(subMesh)


    def addBottom(
        self,
        name: str,
        targetBottomMeshName: str,
        latticeMapFilename: str,
        height: float,
        nz: int=1,
        nSegments: int=1,
        latticeKey: str="F",
        isMergeLatticeElementsInOneMesh: bool=False,
        keyToNotInclude: list=[],
        bottomPatchName: str="",
        isRemoveAllBaffles: bool=False
    ) -> None:
        """
        Add a sub mesh on bottom of the target mesh.

        name: Name of the sub-mesh.
        targetTopMeshName: Name of the target mesh to place the new sub-mesh on bottom.
        latticeMapFilename: Lattice map filename.
        height: Height of the sub-mesh.
        nz: Number of axial cells.
        latticeKey: Filter by lattice keys (string or list of strings).
        isMergeLatticeElementsInOneMesh: If true, places a mesh on each non-0
            element in latticeMapFilename. Else, places a mesh on the specified key.
        keyToNotInclude: Lattice keys to not include. Can be used to reduce the number of baffles.
        bottomPatchName: Name bottom patch of the sub-mesh.
        isRemoveAllBaffles: If true, remove all the internal lateral faces.
        """
        # Found sub mesh with target name
        for subMesh in self.listSubMesh:
            if (subMesh.name == targetBottomMeshName):
                self.listSubMesh.append(
                    SubMesh(
                        name=name,
                        faceShape=self.faceShape,
                        latticeKey=latticeKey,
                        latticeMapFilename=latticeMapFilename,
                        flatToFlat=self.flatToFlat,
                        zCenterPosition=subMesh.zBottomPosition-height/2,
                        height=height,
                        nz=nz,
                        nSegments=nSegments,
                        isMergeLatticeElementsInOneMesh=isMergeLatticeElementsInOneMesh,
                        keyToNotInclude=keyToNotInclude,
                        isRemoveAllBaffles=isRemoveAllBaffles
                    )
                )

                # Remove double faces to avoid overlaping in UNV
                removeGroupByName(self.listSubMesh[-1].mesh, "inlet")

                # Replace or remove bottom patch
                if (bottomPatchName != ""):
                    renameGroupByName(self.listSubMesh[-1].mesh, bottomPatchName, "outlet")

                break


    def addTop(
        self,
        name: str,
        targetTopMeshName: str,
        latticeMapFilename: str,
        height: float,
        nz: int=1,
        nSegments: int=1,
        latticeKey: str="F",
        isMergeLatticeElementsInOneMesh: bool=False,
        keyToNotInclude: list=[],
        topPatchName: str="",
        isRemoveAllBaffles: bool=False
    ) -> None:
        """
        Add a sub mesh on top of the target mesh.

        name: Name of the sub-mesh.
        targetTopMeshName: Name of the target mesh to place the new sub-mesh on top.
        latticeMapFilename: Lattice map filename.
        height: Height of the sub-mesh.
        nz: Number of axial cells.
        latticeKey: Filter by lattice keys (string or list of strings).
        isMergeLatticeElementsInOneMesh: If true, places a mesh on each non-0
            element in latticeMapFilename. Else, places a mesh on the specified key.
        keyToNotInclude: Lattice keys to not include. Can be used to reduce the number of baffles.
        topPatchName: Name top patch of the sub-mesh.
        isRemoveAllBaffles: If true, remove all the internal lateral faces.
        """
        # Found sub mesh with target name
        for subMesh in self.listSubMesh:
            if (subMesh.name == targetTopMeshName):
                self.listSubMesh.append(
                    SubMesh(
                        name=name,
                        faceShape=self.faceShape,
                        latticeKey=latticeKey,
                        latticeMapFilename=latticeMapFilename,
                        flatToFlat=self.flatToFlat,
                        zCenterPosition=subMesh.zCenterPosition+subMesh.height/2+height/2,
                        height=height,
                        nz=nz,
                        nSegments=nSegments,
                        isMergeLatticeElementsInOneMesh=isMergeLatticeElementsInOneMesh,
                        keyToNotInclude=keyToNotInclude,
                        isRemoveAllBaffles=isRemoveAllBaffles
                    )
                )

                # Remove double faces to avoid overlaping in UNV
                removeGroupByName(self.listSubMesh[-1].mesh, "outlet")

                # Replace or remove top patch
                if (topPatchName != ""):
                    renameGroupByName(self.listSubMesh[-1].mesh, topPatchName, "inlet")

                break

    def generateFinalMesh(self) -> None:
        """
        Generate the final mesh for export. Clean the mesh tree in Salome GUI.
        """
        # Create the final mesh
        coreMesh = smesh.Concatenate(
            [subMesh.mesh for subMesh in self.listSubMesh], True, 1, 1e-05, False, name=self.name
        )

        isDone = coreMesh.Compute()

        removeGroupByName(coreMesh, "edges_outlet")
        removeGroupByName(coreMesh, "edges_outlet_top")
        removeGroupByName(coreMesh, "toBeRemoved")

        # Remove temporary meshes
        for subMesh in self.listSubMesh:
            try:
                smesh.RemoveMesh(subMesh.mesh)
                smesh.RemoveMesh(subMesh.elementMesh)
            except:
                None


#=============================================================================*
# Useful functions
#=============================================================================*

def renameGroupByName(mesh, newName: str, oldName: str) -> None:
    """
    Rename all groups corresponding to a name

    Parameters
    ----------
        mesh: Target mesh
        newName: New name
        oldName: Name of the group to change
    """
    for group in mesh.GetGroups():
        if (group.GetName() == oldName):
            group.SetName(newName)


def removeGroupByName(mesh, name: str) -> None:
    """
    Remove all groups corresponding to a name

    mesh: Target mesh
    name: Name of the group to remove
    """
    for group in mesh.GetGroups():
        if (group.GetName() == name):
            mesh.RemoveGroup(group)


def meshElement(mesh, lengthElement: float, nzElement: int, nSegments: int=1) -> None:
    """
    mesh: Face of the mesh to be extruded
    lengthElement: Length of the part of the model in the z-direction
    nzElement: Number of meshes in the z-direction
    nSegments: Number of meshes on a segment of the hexagon

    return an extruded mesh along the Z-axis
    """
    NETGEN_1D_2D = mesh.Triangle(algo=smeshBuilder.NETGEN_1D2D)
    NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
    NETGEN_2D_Parameters_1.SetMaxSize( 0.6 )
    NETGEN_2D_Parameters_1.SetMinSize( 0.3 )
    NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
    NETGEN_2D_Parameters_1.SetOptimize( 1 )
    NETGEN_2D_Parameters_1.SetFineness( 1 )
    NETGEN_2D_Parameters_1.SetChordalError( -1 )
    NETGEN_2D_Parameters_1.SetChordalErrorEnabled( 0 )
    NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
    NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
    NETGEN_2D_Parameters_1.SetWorstElemMeasure( 0 )

    outlet = mesh.GroupOnGeom(Group_1_1, 'outlet', SMESH.FACE)
    edges_outlet = mesh.GroupOnGeom(Group_2_1, 'edges_outlet', SMESH.EDGE)
    Regular_1D = mesh.Segment(geom=Group_2_1)
    Number_of_Segments_1 = Regular_1D.NumberOfSegments(nSegments)
    isDone = mesh.Compute()

    [ outlet, edges_outlet ] = mesh.GetGroups()

    [ extrusion, extruded_Inlet, inlet, edges_inlet ] = mesh.ExtrusionSweepObjects(
        [ mesh ],
        [ mesh ],
        [ mesh ],
        [ 0, 0, lengthElement/nzElement ],
        nzElement, True, [  ], False, [  ], [  ], False
    )

    [ outlet, edges_outlet, extrusion, inlet, extruded_Inlet, edges_inlet ] = mesh.GetGroups()
    Sub_mesh_1 = Regular_1D.GetSubMesh()

    ## Set names of Mesh objects
    extruded_Inlet.SetName('baffle')
    inlet.SetName('inlet')

    return(mesh)


def isElementAtBoundary(latticeMapFilename: str, lineIdx: int, colIdx: int, keyToNotInclude: list=[]) -> bool:
    """
    Test if an element of the lattice is at the lateral boundary, touch a "0"
    element or at the border of the lattice map.

    latticeMapFilename: Lattice map file name
    lineIdx: Line index of the element to test
    colIdx: Column index of the element to test

    return True if the element at indices (lineIdx, colIdx) is at the boundary
    """
    # Hexagonal positions of the neighbors
    lineIdxs = [lineIdx-1, lineIdx-1, lineIdx, lineIdx, lineIdx+1, lineIdx+1]
    colIdxs = [colIdx, colIdx+1, colIdx-1, colIdx+1, colIdx-1, colIdx]

    # Open the lattice
    with open(latticeMapFilename, 'r') as file:
        lattice = [line.split() for line in file.readlines()]

        # Edge cases
        if (lineIdx == 0 or len(lattice) == lineIdx+1):
            return(True)
        if (colIdx == 0 or len(lattice[0]) == colIdx+1):
            return(True)

        elements = [lattice[row][col] for row, col in zip(lineIdxs, colIdxs)]

        for boundaryKey in keyToNotInclude+["0"]:
            if (boundaryKey in elements):
                return(True)

        return(False)


def createMeshWithLattice(
    mesh,
    latticeMapFilename: str,
    flatToFlat: float,
    z: float=0,
    key: str='F',
    isMergeLatticeElementsInOneMesh: bool=False,
    keyToNotInclude: list=[],
    name: str="mesh",
    isRemoveAllBaffles: bool=False
):
    """
    mesh: part of the model one want to mesh
    latticeMapFilename: Lattice map file name
    flatToFlat: Flat-to-flat distance between the elements (also called picth).
    z: z-coordinate of the bottom of the mesh
    key: Filter by lattice keys (string or list of strings).
    isMergeLatticeElementsInOneMesh: If true, places a mesh on each non-0
        element in latticeMapFilename. Else, places a mesh on the specified key.
    name: Name of the mesh
    isRemoveAllBaffles: If true, remove all the lateral patches
    return concatenated mesh
    """
    summitToSummit = flatToFlat * 2/np.sqrt(3)
    side = summitToSummit/2

    countElements = 0
    meshElements = []  # Creating list to store each type of element

    # Opening of the lattice map
    with open(latticeMapFilename, 'r') as file:
        lines = file.readlines()
        # First position in latticeMapFilename corresponds to (x,y,z) in Salome
        nx, ny = len(lines), len(lines[0].split())
        x = -1.5*side * (nx-1)/2

        # Loop through the lines
        for lineIdx, line in enumerate(lines):

            # To place the first element of a line according to its location in latticeMapFilename
            y = flatToFlat/2 * (lineIdx - 1.5*(ny-1))

            # Loop over the columns
            for colIdx, element in enumerate(line.split()):
                # Never place 0 elements
                if (
                    element != "0"
                    and (
                        (not isMergeLatticeElementsInOneMesh and element in key)
                        or
                        (isMergeLatticeElementsInOneMesh and element not in keyToNotInclude)
                    )
                ):
                    meshTranslated = mesh.TranslateObjectMakeMesh(
                        mesh, [x, y, z], 1, name
                    )
                    countElements += 1

                    # Special name to remove the baffles
                    if (isElementAtBoundary(latticeMapFilename, lineIdx, colIdx, keyToNotInclude) or isRemoveAllBaffles):
                        # Removing lateral face at the boundary of the mesh to
                        # make it work with ideasUNVToFoam
                        removeGroupByName(meshTranslated, "baffle")

                    # Append in the list for concatenation
                    meshElements.append(meshTranslated.GetMesh())

                # Updating y for next element
                y += flatToFlat

            # End of line, updating x
            x += summitToSummit/2 + side/2

    # Concatenate only non-empty list
    if (len(meshElements) > 0):
        # Merging of all the groups stored in the meshElements list
        meshConcat = smesh.Concatenate(
            meshElements, True, True, 1e-05, False, name=name
        )

        # Rename cellZone
        renameGroupByName(meshConcat, name, "outlet_extruded")

        # Remove translated meshes to clean Salome mesh tree in GUI
        for part in meshElements:
            smesh.RemoveMesh(part)

        print(f"Number of elements placed for key {key}: {countElements}")

        return(meshConcat.GetMesh())


def createHexagon(flatToFlat: float):
    """
    Create hexagonal shape in Shaper workbench of Salome

    flatToFlat: Flat-to-flat distance between the elements (also called picth).
    """
    model.begin()
    partSet = model.moduleDocument()

    ### Create Part
    Part_1 = model.addPart(partSet)
    Part_1_doc = Part_1.document()

    ### Create Sketch
    Sketch_1 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

    SketchProjection_1 = Sketch_1.addProjection(model.selection("VERTEX", "PartSet/Origin"), False)
    SketchPoint_1 = SketchProjection_1.createdFeature()

    ### Create SketchCircle
    SketchCircle_1 = Sketch_1.addCircle(0, 0, 71e-2) # X, Y, Radius
    SketchCircle_1.setAuxiliary(True)
    Sketch_1.setCoincident(SketchCircle_1.center(), SketchAPI_Point(SketchPoint_1).coordinates())

    ### Create polyline representing the hexagon
    SketchLine_1 = Sketch_1.addLine(35e-2, 61e-2, -35e-2, 61e-2)
    Sketch_1.setCoincident(SketchLine_1.startPoint(), SketchCircle_1.results()[1])
    Sketch_1.setCoincident(SketchLine_1.endPoint(), SketchCircle_1.results()[1])
    SketchLine_2 = Sketch_1.addLine(-35e-2, 61e-2, -71e-2, 0)
    Sketch_1.setCoincident(SketchLine_1.endPoint(), SketchLine_2.startPoint())
    Sketch_1.setCoincident(SketchLine_2.endPoint(), SketchCircle_1.results()[1])
    SketchLine_3 = Sketch_1.addLine(-71e-2, 0, -35e-2, -61e-2)
    Sketch_1.setCoincident(SketchLine_2.endPoint(), SketchLine_3.startPoint())
    Sketch_1.setCoincident(SketchLine_3.endPoint(), SketchCircle_1.results()[1])
    SketchLine_4 = Sketch_1.addLine(-35e-2, -61e-2, 35e-2, -61e-2)
    Sketch_1.setCoincident(SketchLine_3.endPoint(), SketchLine_4.startPoint())
    Sketch_1.setCoincident(SketchLine_4.endPoint(), SketchCircle_1.results()[1])
    SketchLine_5 = Sketch_1.addLine(35e-2, -61e-2, 71e-2, -0)
    Sketch_1.setCoincident(SketchLine_4.endPoint(), SketchLine_5.startPoint())
    Sketch_1.setCoincident(SketchLine_5.endPoint(), SketchCircle_1.results()[1])
    SketchLine_6 = Sketch_1.addLine(71e-2, 0, 35e-2, 61e-2)
    Sketch_1.setCoincident(SketchLine_5.endPoint(), SketchLine_6.startPoint())
    Sketch_1.setCoincident(SketchLine_1.startPoint(), SketchLine_6.endPoint())

    ### Set constraints

    # Horizontal
    Sketch_1.setHorizontal(SketchLine_1.result())

    # Set equal distance of the edges
    Sketch_1.setEqual(SketchLine_1.result(), SketchLine_6.result())
    Sketch_1.setEqual(SketchLine_1.result(), SketchLine_5.result())
    Sketch_1.setEqual(SketchLine_1.result(), SketchLine_4.result())
    Sketch_1.setEqual(SketchLine_1.result(), SketchLine_3.result())
    Sketch_1.setEqual(SketchLine_1.result(), SketchLine_2.result())

    # Set flat to flat distance
    Sketch_1.setDistance(SketchLine_1.result(), SketchLine_3.endPoint(), flatToFlat, True)

    model.do()


    #=============================================================================*

    ### Create Face
    Face_1 = model.addFace(Part_1_doc, [model.selection("FACE", "Sketch_1/Face-SketchLine_1f-SketchLine_2f-SketchLine_3f-SketchLine_4f-SketchLine_5f-SketchLine_6f")])

    ### Create Group
    Group_1 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Face_1_1")])

    ### Create Group
    Group_2_objects = [
        model.selection("EDGE", "Face_1_1/Modified_Edge&Sketch_1/SketchLine_6"),
        model.selection("EDGE", "Face_1_1/Modified_Edge&Sketch_1/SketchLine_5"),
        model.selection("EDGE", "Face_1_1/Modified_Edge&Sketch_1/SketchLine_1"),
        model.selection("EDGE", "Face_1_1/Modified_Edge&Sketch_1/SketchLine_4"),
        model.selection("EDGE", "Face_1_1/Modified_Edge&Sketch_1/SketchLine_2"),
        model.selection("EDGE", "Face_1_1/Modified_Edge&Sketch_1/SketchLine_3")
    ]
    Group_2 = model.addGroup(Part_1_doc, "Edges", Group_2_objects)

    model.end()

    model.publishToShaperStudy()
    Face_1_1, Group_1_1, Group_2_1, = SHAPERSTUDY.shape(model.featureStringId(Face_1))

    return(Face_1_1, Group_1_1, Group_2_1)


smesh = smeshBuilder.New()


#=============================================================================*
#=============================================================================*
# Beginning of MAIN for Users
#=============================================================================*

# Units
inch = 2.54e-2 # m

coreDiameter = 35 * inch

# Dimensions of the hexagonal face of the element
flatToFlat = 0.75 * inch # 1.905 cm
summitToSummit = flatToFlat * 2/np.sqrt(3)
side = summitToSummit/2

### Number of meshes in z-direction for each par of the model
nzFuelElement   = 10 * 2 # Need to increase to 20 axial cells for journal paper
nzNozzleChamber =  2 * 4
nzCoreSupport   =  1 * 2
nzSupportPlate  =  1 * 2
nzUpperPlenum   =  2 * 4
nzSpacer        =  1
# nzGap           =  1 * 4

nSegmentsNeutro = 1
nSegmentsFluid  = 1

### Dimensions found in : 1965 - Survey description of the design and testing of Kiwi-B-4E-301
lengthFuelElement   = 52.0 * inch
lengthNozzleChamber = 10.0 * inch  ## arbitrary length to study what happens in this zone
lengthCoreSupportGr =  1.5 * inch  ## Graphite
lengthSupportPlateAl=  4.5 * inch  ## Aluminum
lengthSpacer        =  5.812 * inch - lengthSupportPlateAl
lengthUpperPlenum   =  5.0 * inch  ## arbitrary length to study what happens in this zone
# lengthGap           =  1.2 * inch


# Absolute path of the lattice map file. This file contains the pattern of the
# elements
latticeMapFilename = "/home/thomas-guilbaud/Simulation/GeN-Foam/gitlab-nerva/nerva/KiwiB4E/core/GeN-Foam/mesh/coreMap.txt"


hexagonFace, Group_1_1, Group_2_1 = createHexagon(flatToFlat)

# Create mesh final object
neutronicsMesh = LargeMesh(name="neutroMesh", flatToFlat=flatToFlat, faceShape=hexagonFace)
thermalHydraulicsMesh = LargeMesh(name="fluidMesh", flatToFlat=flatToFlat, faceShape=hexagonFace)

# Add sub-meshes for neutronics region
for name, latticeKey in zip(
    ["45", "39", "34", "30", "26", "22", "19", "17", "15", "13"],
    ["F", "9", "8", "7", "6", "5", "4", "3", "2", "1"]
):
    neutronicsMesh.createSubMesh(
        name="fuelElement_"+name,
        latticeMapFilename=latticeMapFilename,
        height=lengthFuelElement,
        nz=nzFuelElement,
        nSegments=nSegmentsNeutro,
        latticeKey=latticeKey,
        keyToNotInclude=[],
        isRemoveAllBaffles=True
    )

# Unloaded elements
neutronicsMesh.createSubMesh(
    name="unloadedElementsTypeC",
    latticeMapFilename=latticeMapFilename,
    height=lengthFuelElement,
    nz=nzFuelElement,
    nSegments=nSegmentsNeutro,
    latticeKey="C",
    keyToNotInclude=[],
    isRemoveAllBaffles=True
)
neutronicsMesh.createSubMesh(
    name="unloadedElementsTypeT",
    latticeMapFilename=latticeMapFilename,
    height=lengthFuelElement,
    nz=nzFuelElement,
    nSegments=nSegmentsNeutro,
    latticeKey="T",
    keyToNotInclude=[],
    isRemoveAllBaffles=True
)

# Reflectors
neutronicsMesh.createSubMesh(
    name="reflector",
    latticeMapFilename=latticeMapFilename,
    height=lengthFuelElement,
    nz=nzFuelElement,
    nSegments=nSegmentsNeutro,
    latticeKey="g",
    keyToNotInclude=[],
    isRemoveAllBaffles=True
)
neutronicsMesh.createSubMesh(
    name="reflectorSystemBe",
    latticeMapFilename=latticeMapFilename,
    height=lengthFuelElement,
    nz=nzFuelElement,
    nSegments=nSegmentsNeutro,
    latticeKey="R",
    keyToNotInclude=[],
    isRemoveAllBaffles=True,
    # topPatchName="topReflectorSystem",
    # bottomPatchName="bottomReflectorSystem"
)
neutronicsMesh.createSubMesh(
    name="controlDrum",
    latticeMapFilename=latticeMapFilename,
    height=lengthFuelElement,
    nz=nzFuelElement,
    nSegments=nSegmentsNeutro,
    latticeKey="D",
    keyToNotInclude=[],
    isRemoveAllBaffles=True,
    # topPatchName="topReflectorSystem",
    # bottomPatchName="bottomReflectorSystem"
)

# Bottom core
neutronicsMesh.addBottom(
    name='coreSupport',
    targetBottomMeshName="fuelElement_45",
    latticeMapFilename=latticeMapFilename,
    height=lengthCoreSupportGr,
    nz=nzCoreSupport,
    nSegments=nSegmentsNeutro,
    isMergeLatticeElementsInOneMesh=True,
    keyToNotInclude=["R", "D"],
    bottomPatchName="toBeRemoved",
    isRemoveAllBaffles=True
)
neutronicsMesh.addBottom(
    name='nozzleChamber',
    targetBottomMeshName="coreSupport",
    latticeMapFilename=latticeMapFilename,
    height=lengthNozzleChamber,
    nz=nzNozzleChamber,
    nSegments=nSegmentsNeutro,
    isMergeLatticeElementsInOneMesh=True,
    keyToNotInclude=["R", "D"],
    isRemoveAllBaffles=True
)

# Inlet hydrogen
neutronicsMesh.addBottom(
    name='reflectorInletPropellant',
    targetBottomMeshName="reflectorSystemBe",
    latticeMapFilename=latticeMapFilename,
    height=lengthCoreSupportGr,
    nz=nzCoreSupport,
    nSegments=nSegmentsNeutro,
    latticeKey=["R", "D"],
    bottomPatchName="toBeRemoved",
    isRemoveAllBaffles=True
)
neutronicsMesh.addBottom(
    name='reflectorInletPropellant',
    targetBottomMeshName="reflectorInletPropellant",
    latticeMapFilename=latticeMapFilename,
    height=lengthNozzleChamber,
    nz=nzNozzleChamber,
    nSegments=nSegmentsNeutro,
    latticeKey=["R", "D"],
    isRemoveAllBaffles=True
)

# Hydrogen U-flow
neutronicsMesh.addTop(
    name='inletPropellant',
    targetTopMeshName="fuelElement_45",
    latticeMapFilename=latticeMapFilename,
    height=lengthSpacer,
    nz=nzSpacer,
    nSegments=nSegmentsNeutro,
    isMergeLatticeElementsInOneMesh=True,
    keyToNotInclude=["g", "R", "D"],
    topPatchName="toBeRemoved",
    isRemoveAllBaffles=True
)
neutronicsMesh.addTop(
    name='supportPlate',
    targetTopMeshName="inletPropellant",
    latticeMapFilename=latticeMapFilename,
    height=lengthSupportPlateAl,
    nz=nzSupportPlate,
    nSegments=nSegmentsNeutro,
    isMergeLatticeElementsInOneMesh=True,
    keyToNotInclude=["g", "R", "D"],
    topPatchName="toBeRemoved",
    isRemoveAllBaffles=True
)
neutronicsMesh.addTop(
    name='upperPropellantPlenum',
    targetTopMeshName="supportPlate",
    latticeMapFilename=latticeMapFilename,
    height=lengthSupportPlateAl,
    nz=nzSupportPlate,
    nSegments=nSegmentsNeutro,
    isMergeLatticeElementsInOneMesh=True,
    isRemoveAllBaffles=True
)
neutronicsMesh.addTop(
    name='reflectorCornerPropellant',
    targetTopMeshName="reflectorSystemBe",
    latticeMapFilename=latticeMapFilename,
    height=lengthSpacer,
    nz=nzSpacer,
    nSegments=nSegmentsNeutro,
    latticeKey=["R", "D"],
    topPatchName="toBeRemoved",
    isRemoveAllBaffles=True
)
neutronicsMesh.addTop(
    name='inletSpacer',
    targetTopMeshName="reflector",
    latticeMapFilename=latticeMapFilename,
    height=lengthSpacer,
    nz=nzSpacer,
    nSegments=nSegmentsNeutro,
    latticeKey="g",
    topPatchName="toBeRemoved",
    isRemoveAllBaffles=True
)
neutronicsMesh.addTop(
    name='reflectorCornerPropellant',
    targetTopMeshName="inletSpacer",
    latticeMapFilename=latticeMapFilename,
    height=lengthSupportPlateAl,
    nz=nzSupportPlate,
    nSegments=nSegmentsNeutro,
    latticeKey=["g", "R", "D"],
    topPatchName="toBeRemoved",
    isRemoveAllBaffles=True
)



# Add sub-meshes for thermal-hydraulics region
thermalHydraulicsMesh.createSubMesh(
    name='fuelElements',
    latticeMapFilename=latticeMapFilename,
    height=lengthFuelElement,
    nz=nzFuelElement,
    nSegments=nSegmentsFluid,
    latticeKey=["F", "1", "2", "3", "4", "5", "6", "7", "8", "9"],
    keyToNotInclude=["g", "R", "D"],
    topPatchName="inletFuelElement",
    bottomPatchName="outletFuelElement"
)
thermalHydraulicsMesh.createSubMesh(
    name="unloadedElements",
    latticeMapFilename=latticeMapFilename,
    height=lengthFuelElement,
    nz=nzFuelElement,
    nSegments=nSegmentsFluid,
    latticeKey=["C", "T"],
    keyToNotInclude=["g", "R", "D"],
    topPatchName="inletCentralUnloaded",
    bottomPatchName="outletCentralUnloaded"
)
thermalHydraulicsMesh.addTop(
    name='inletSpacer',
    targetTopMeshName="fuelElements",
    latticeMapFilename=latticeMapFilename,
    height=lengthSpacer,
    nz=nzSpacer,
    nSegments=nSegmentsFluid,
    isMergeLatticeElementsInOneMesh=True,
    keyToNotInclude=["g", "R", "D"],
    isRemoveAllBaffles=True
)
thermalHydraulicsMesh.addBottom(
    name='coreSupport',
    targetBottomMeshName="fuelElements",
    latticeMapFilename=latticeMapFilename,
    height=lengthCoreSupportGr,
    nz=nzCoreSupport,
    nSegments=nSegmentsFluid,
    isMergeLatticeElementsInOneMesh=True,
    keyToNotInclude=["g", "R", "D"],
    bottomPatchName="toBeRemoved",
    isRemoveAllBaffles=False
)
thermalHydraulicsMesh.addBottom(
    name='nozzleChamber',
    targetBottomMeshName="coreSupport",
    latticeMapFilename=latticeMapFilename,
    height=lengthNozzleChamber,
    nz=nzNozzleChamber,
    nSegments=nSegmentsFluid,
    isMergeLatticeElementsInOneMesh=True,
    keyToNotInclude=["g", "R", "D"],
    isRemoveAllBaffles=True
)

# Generate the final mesh
# neutronicsMesh.generateFinalMesh()
thermalHydraulicsMesh.generateFinalMesh()


#=============================================================================*
# End of MAIN
#=============================================================================*
#=============================================================================*

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()

#=============================================================================*
