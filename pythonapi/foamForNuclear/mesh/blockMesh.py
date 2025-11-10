import json
import numpy as np
import tqdm

from foamForNuclear.checkvalue import check_positive, check_type, check_value

from .mesh import _LATTICE_TYPES, Mesh
from foamForNuclear.openfoamFile import OpenFOAMFile
from foamForNuclear.common import *
import math

_CYCLIC_PATCH_TRANSFORM_TYPES = {
    "unknown", "rotational", "translational", "coincidentFullMatch", "noOrdering"
}
_EDGE_TYPES = {"arc", "polyLine"}
_FACE_NAME_TYPES = {'top', 'bottom', 'left', 'right', 'front', 'back'}

_SAMPLE_MODE_TYPES = {
    "nearestCell", "nearestOnlyCell", "nearestPatchFace", "nearestPatchFaceAMI",
    "nearestFace", "nearestPatchPoint"
}
_OFFSET_MODE_TYPES = {"uniform", "nonuniform", "normal"}

idPoint = 0


class Point:
    def __init__(self, x: float, y: float, z: float, isIndexed: bool=True):
        self.x: float = x
        self.y: float = y
        self.z: float = z
        self.id: int = -1
        self.isIndexed = isIndexed
        if (isIndexed):
            self.setIndex()

    def __repr__(self):
        return((f"/* {self.id:3} */ " if self.id >= 0 else "") + f"( {self.x} {self.y} {self.z} )")

    def __eq__(self, value):
        eps = 1e-6
        return(
            abs(self.x - value.x) <= eps and
            abs(self.y - value.y) <= eps and
            abs(self.z - value.z) <= eps
        )

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, x) -> None:
        check_type("x", x, (float, int))
        self._x = x

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, y) -> None:
        check_type("y", y, (float, int))
        self._y = y

    @property
    def z(self):
        return self._z

    @z.setter
    def z(self, z) -> None:
        check_type("z", z, (float, int))
        self._z = z

    @property
    def isIndexed(self):
        return self._isIndexed

    @isIndexed.setter
    def isIndexed(self, isIndexed) -> None:
        check_type("isIndexed", isIndexed, bool)
        self._isIndexed = isIndexed

    def setIndex(self) -> None:
        global idPoint
        self.isIndexed = True
        self.id = idPoint if self.isIndexed else -1
        if (self.isIndexed):
            idPoint += 1

    def translate(self, dx: float=0, dy: float=0, dz: float=0) -> None:
        self.x += dx
        self.y += dy
        self.z += dz

    def rotateX(self, theta: float=0) -> None:
        py = self.y
        self.y =  py * np.cos(theta) - self.z * np.sin(theta)
        self.z =  py * np.sin(theta) + self.z * np.cos(theta)

    def rotateY(self, theta: float=0) -> None:
        px = self.x
        self.x =  px * np.cos(theta) - self.z * np.sin(theta)
        self.z =  px * np.sin(theta) + self.z * np.cos(theta)

    def rotateZ(self, theta: float=0) -> None:
        px = self.x
        self.x =  px * np.cos(theta) - self.y * np.sin(theta)
        self.y =  px * np.sin(theta) + self.y * np.cos(theta)


def isEqualSubFace(face1: list[Point], face2: list[Point]) -> bool:
    n = len(face1)
    count = 0
    for p1 in face1:
        for p2 in face2:
            if (p1 == p2):
                count += 1
                break
    return(n == count)


class Edge:
    def __init__(
            self,
            edgeType: str,
            pointId1: Point,
            pointId2: Point,
        ):
        self.edgeType: str = edgeType
        self.pointId1: Point = pointId1
        self.pointId2: Point = pointId2


    @property
    def edgeType(self):
        return self._edgeType

    @edgeType.setter
    def edgeType(self, edgeType) -> None:
        check_type("edgeType", edgeType, str)
        check_value("edgeType", edgeType, _EDGE_TYPES)
        self._edgeType = edgeType

    @property
    def pointId1(self):
        return self._pointId1

    @pointId1.setter
    def pointId1(self, pointId1: Point) -> None:
        check_type("pointId1", pointId1, Point)
        self._pointId1 = pointId1

    @property
    def pointId2(self):
        return self._pointId2

    @pointId2.setter
    def pointId2(self, pointId2: Point) -> None:
        check_type("pointId2", pointId2, Point)
        self._pointId2 = pointId2



class EdgeArc(Edge):
    def __init__(self, pointId1, pointId2, midPoint, isOrigin = False):
        super().__init__("arc", pointId1, pointId2)

        self.midPoint: Point = midPoint
        self.isOrigin: bool = isOrigin

    def __repr__(self):
        return(f"{self.edgeType} {self.pointId1.id} {self.pointId2.id} {'origin' if self.isOrigin else ''} {self.midPoint}")

    @property
    def midPoint(self):
        return self._midPoint

    @midPoint.setter
    def midPoint(self, midPoint) -> None:
        check_type("midPoint", midPoint, Point)
        self._midPoint = midPoint

    @property
    def isOrigin(self):
        return self._isOrigin

    @isOrigin.setter
    def isOrigin(self, isOrigin) -> None:
        check_type("isOrigin", isOrigin, bool)
        self._isOrigin = isOrigin


class EdgePolyLine(Edge):
    def __init__(self, pointId1, pointId2, points: list[Vector]):
        super().__init__("polyLine", pointId1, pointId2)

        self.points = points

    def __repr__(self):
        text = f"{self.edgeType} {self.pointId1.id} {self.pointId2.id}\n"
        text += tab+"(\n"
        for point in self.points:
            text += f"{2*tab}{point}\n"
        text += tab+")"
        return(text)


class Face:
    def __init__(
            self,
            name: str,
            boundaryType: str="patch",
            isPrint: bool=True,
            extraParameters: dict={},
            inGroups: list=[]
        ):
        self.name: str = name
        self.inGroups = inGroups
        self.boundaryType: str = boundaryType
        self.faces: list[list[Point]] = []
        self.isPrint = isPrint
        self.extraParameters = extraParameters

    def __repr__(self):
        txt = f"{tab}{self.name}\n"
        txt += tab+"{\n"
        txt +=f"{2*tab}type {self.boundaryType};\n"

        if (len(self.inGroups) > 0):
            txt +=f"{2*tab}inGroups ({' '.join(self.inGroups)});\n"

        for key, item in self.extraParameters.items():
            txt +=f"{2*tab}{key:15} {item};\n"

        txt += f"{2*tab}faces\n"
        txt += f"{2*tab}(\n"
        for face in self.faces:
            txt += f"{3*tab}( {' '.join([str(e.id) for e in face])} )\n"
        txt += f"{2*tab});\n"
        txt += tab + "}"
        return(txt)

    @property
    def is_empty(self):
        return(len(self.faces) == 0)

    def addSubFace(self, face: list[Point]) -> None:
        self.faces.append(face)

    def getCommonSubFaces(self, targetFace):
        common = []
        for facei in self.faces:
            for facej in targetFace.faces:
                if (self.isEqualSubFace(facei, facej)):
                    common.append((facei, facej))
        return(common)


class FaceCyclic(Face):
    """
    Parameters
    ----------
    transform : str
        Options: "unknown", "rotational", "translational", "coincidentFullMatch",
        "noOrdering"
    separationVector : Vector
        Only used when transform == "translational"
    """

    def __init__(
            self,
            name: str,
            neighbourPatch: str,
            transform: str,
            rotationAxis: Vector=None,
            rotationCentre: Vector=None,
            separationVector: Vector=None,
            isPrint: bool=True,
            inGroups: list=[]
        ):
        super().__init__(name, "discontinuousCyclicAMI", isPrint, {}, inGroups)

        self.neighbourPatch = neighbourPatch
        self.transform = transform
        self.rotationAxis = rotationAxis
        self.rotationCentre = rotationCentre
        self.separationVector = separationVector
        self.extraParameters['matchTolerance'] = 0.1

    @property
    def neighbourPatch(self):
        return self._neighbourPatch

    @neighbourPatch.setter
    def neighbourPatch(self, neighbourPatch) -> None:
        check_type("neighbourPatch", neighbourPatch, str)
        self._neighbourPatch = neighbourPatch
        self.extraParameters['neighbourPatch'] = neighbourPatch

    @property
    def transform(self):
        return self._transform

    @transform.setter
    def transform(self, transform) -> None:
        check_type("transform", transform, str)
        check_value("transform", transform, _CYCLIC_PATCH_TRANSFORM_TYPES)
        self._transform = transform
        self.extraParameters['transform'] = transform

    @property
    def rotationAxis(self):
        return self._rotationAxis

    @rotationAxis.setter
    def rotationAxis(self, rotationAxis) -> None:
        check_type("rotationAxis", rotationAxis, Vector, none_ok=True)
        self._rotationAxis = rotationAxis
        # Only add if rotational
        if self.transform == "rotational" and rotationAxis is not None:
            self.extraParameters['rotationAxis'] = f"{rotationAxis}"
        elif 'rotationAxis' in self.extraParameters:
            del self.extraParameters['rotationAxis']

    @property
    def rotationCentre(self):
        return self._rotationCentre

    @rotationCentre.setter
    def rotationCentre(self, rotationCentre) -> None:
        check_type("rotationCentre", rotationCentre, Vector, none_ok=True)
        self._rotationCentre = rotationCentre
        # Only add if rotational
        if self.transform == "rotational" and rotationCentre is not None:
            self.extraParameters['rotationCentre'] = f"{rotationCentre}"
        elif 'rotationCentre' in self.extraParameters:
            del self.extraParameters['rotationCentre']

    @property
    def separationVector(self):
        return self._separationVector

    @separationVector.setter
    def separationVector(self, separationVector) -> None:
        # ✅ Only meaningful for translational
        check_type("separationVector", separationVector, Vector, none_ok=True)
        self._separationVector = separationVector
        if self.transform == "translational" and separationVector is not None:
            self.extraParameters['separationVector'] = f"{separationVector}"
        elif 'separationVector' in self.extraParameters:
            del self.extraParameters['separationVector']

class FaceMappedPatch(Face):
    """
    Determines a mapping between patch face centres and mesh cell or face
    centres and processors they're on.

    Parameters
    ----------
    name : str
        Name of the patch.
    sampleRegion : str
        Name of the region to sample
    sampleMode : str
        - `nearestCell`: sample cell containing point
        - `nearestOnlyCell`: nearest sample cell (even if not containing point)
        - `nearestPatchFace`: nearest face on selected patch
        - `nearestPatchFaceAMI`: nearest face on selected patch (patches need
            not conform; uses AMI interpolation)
        - `nearestFace`: nearest boundary face on any patch
        - `nearestPatchPoint`: nearest patch point (for coupled points this
            might be any of the points so you have to guarantee the point data
            is synchronised beforehand)
    samplePatch  : str
        If sampleMode is `nearestPatchFace`: patch to find faces of
    offsetMode : str
        Default `uniform`.
        How to supply offset (w.r.t. my patch face centres):
        - `uniform`: single offset vector
        - `nonuniform`: per-face offset vector
        - `normal`: using supplied distance and face normal
    offset : Vector
        According to `offsetMode` (see above) supply one of offset, offsets or
        distance (default `Vector(0, 0, 0)`)
    coupleGroup : str
        If sampleMode is `nearestPatchFace`: specify patchgroup to find
        `samplePatch` and `sampleRegion` (if not provided)
    """
    def __init__(
            self,
            name: str,
            sampleRegion: str,
            sampleMode: str,
            samplePatch: str,
            offsetMode: str='uniform',
            offset: Vector=Vector(0, 0, 0),
            coupleGroup: str=None,
            isPrint: bool=True,
            inGroups: list=[]
        ):
        super().__init__(name, "mappedPatch", isPrint, {}, inGroups)

        self.sampleRegion = sampleRegion
        self.sampleMode = sampleMode
        self.samplePatch = samplePatch
        self.coupleGroup = coupleGroup
        self.offsetMode = offsetMode
        self.offset = offset

    @property
    def sampleRegion(self):
        return self._sampleRegion

    @sampleRegion.setter
    def sampleRegion(self, sampleRegion) -> None:
        check_type("sampleRegion", sampleRegion, str)
        self._sampleRegion = sampleRegion
        if (sampleRegion is not None):
            self.extraParameters['sampleRegion'] = sampleRegion

    @property
    def sampleMode(self):
        return self._sampleMode

    @sampleMode.setter
    def sampleMode(self, sampleMode) -> None:
        check_type("sampleMode", sampleMode, str)
        self._sampleMode = sampleMode
        if (sampleMode is not None):
            check_value("sampleMode", sampleMode, _SAMPLE_MODE_TYPES)
            self.extraParameters['sampleMode'] = sampleMode

    @property
    def samplePatch(self):
        return self._samplePatch

    @samplePatch.setter
    def samplePatch(self, samplePatch) -> None:
        check_type("samplePatch", samplePatch, str)
        self._samplePatch = samplePatch
        if (samplePatch is not None):
            self.extraParameters['samplePatch'] = samplePatch

    @property
    def coupleGroup(self):
        return self._coupleGroup

    @coupleGroup.setter
    def coupleGroup(self, coupleGroup) -> None:
        check_type("coupleGroup", coupleGroup, str, none_ok=True)
        self._coupleGroup = coupleGroup
        if (coupleGroup is not None):
            self.extraParameters['coupleGroup'] = coupleGroup

    @property
    def offsetMode(self):
        return self._offsetMode

    @offsetMode.setter
    def offsetMode(self, offsetMode) -> None:
        check_type("offsetMode", offsetMode, str, none_ok=True)
        self._offsetMode = offsetMode
        if (offsetMode is not None):
            check_value("offsetMode", offsetMode, _OFFSET_MODE_TYPES)
            self.extraParameters['offsetMode'] = offsetMode

    @property
    def offset(self):
        return self._offset

    @offset.setter
    def offset(self, offset) -> None:
        check_type("offset", offset, Vector, none_ok=True)
        self._offset = offset
        if (offset is not None):
            self.extraParameters['offset'] = f"{offset}"


class Block:
    """
    A Block is defined by 8 points forming a hexahedron.

    Parameters
    ----------
    name : str
        Name of the Block
    points : list[Point]
        List of Point objects
    nx : int
        Number of discretization along the X-direction (default `1`).
    ny : int=1
        Number of discretization along the Y-direction (default `1`).
    nz : int=1
        Number of discretization along the Z-direction (default `1`).
    isPrint : bool
        Flag to allow printing the object in the final blockMeshDict
        (default `True`).

    Attributes
    ----------
    edges : list[Edge]
        List of Edge objects.
    faceProjection : list
        List of face projection. Updated using `addFaceProjection`
    edgeProjection : list
        List of edge projection. Updated using `addEdgeProjection`
    """
    def __init__(
            self,
            name: str,
            points: list[Point],
            nx: int=1,
            ny: int=1,
            nz: int=1,
            gradx: float=1,
            grady: float=1,
            gradz: float=1,
            isPrint: bool=True
        ):
        self.name: str = name
        self.points: list[Point] = points
        self.edges: list[Edge] = []
        self.faceProjection: list = []
        self.edgeProjection: list = []
        self.nx: int = nx
        self.ny: int = ny
        self.nz: int = nz
        self.gradx: float = gradx
        self.grady: float = grady
        self.gradz: float = gradz
        self.isPrint = isPrint

    def __repr__(self):
        return(f"hex ({' '.join([str(p.id) for p in self.points])}) {self.name} ({self.nx} {self.ny} {self.nz}) simpleGrading ({self.gradx} {self.grady} {self.gradz})")

    def addEdgeArc(
            self,
            pointIdx1: int, pointIdx2: int,
            x: float=None, y: float=None, z: float=None,
            isMidPointIndexed: bool=False,
            isOrigin: bool=False
        ):
        """
        Deform an edge to be an arc.

        Parameters
        ----------
        x : float
            X-position, if `x` is `None`, use X-coordinates of `pointIdx1`.
        y : float
            Y-position, if `y` is `None`, use Y-coordinates of `pointIdx1`.
        z : float
            Z-position, if `z` is `None`, use Z-coordinates of `pointIdx1`.
        isOrigin : bool
            If true, use `x`, `y`, and `z` as the coordinates of the arc center
            of rotation. If false, use `x`, `y`, and `z` as the coordinates of
            the third point where the arc needs to pass (default `False`).
        """
        check_type("pointIdx1", pointIdx1, int)
        check_type("pointIdx2", pointIdx2, int)

        if (x == None):
            x = self.points[pointIdx1].x
        if (y == None):
            y = self.points[pointIdx1].y
        if (z == None):
            z = self.points[pointIdx1].z

        self.edges.append(EdgeArc(
            pointId1=self.points[pointIdx1],
            pointId2=self.points[pointIdx2],
            midPoint=Point(x, y, z, isIndexed=isMidPointIndexed),
            isOrigin=isOrigin
        ))

    def addEdgePolyLine(
            self,
            pointIdx1: int,
            pointIdx2: int,
            points: list[Vector]
        ):
        """
        Deform an edge to follow a polyline
        """
        check_type("pointIdx1", pointIdx1, int)
        check_type("pointIdx2", pointIdx2, int)

        self.edges.append(EdgePolyLine(
            pointId1=self.points[pointIdx1],
            pointId2=self.points[pointIdx2],
            points=points
        ))

    def addFaceProjection(self, faceName: int, geometryName: str) -> None:
        """
        Deform a face using a projection geometry.
        """
        faceNameToId = {
            "left": 0, "right": 1, "front": 2, "back": 3, "bottom": 4, "top": 5
        }
        self.faceProjection.append({
            "faceId": faceNameToId[faceName], "geometryName": geometryName
        })

    def addEdgeProjection(
            self,
            verticeIdx1: int,
            verticeIdx2: int,
            geometryNames: list[str]
        ) -> None:
        """
        Deform an edge using a projection geometry.
        """
        check_type("geometryNames", geometryNames, list)

        self.edgeProjection.append({
            'verticeIdx1': self.points[verticeIdx1],
            'verticeIdx2': self.points[verticeIdx2],
            'geometryNames': ' '.join(geometryNames)
        })

    def translate(self, dx: float=0, dy: float=0, dz: float=0) -> None:
        for point in self.points:
            point.translate(dx, dy, dz)

    def rotateZ(self, theta: float=0) -> None:
        """
        Rotation of all the points coordinates along the Z-axis.

        Parameters
        ----------
        theta : float
            Rotation angle in rad (default `0`).
        """
        for point in self.points:
            point.rotateZ(theta)

    def frontFace(self) -> list[Point]:
        """
        Points [0, 1, 5, 4]
        """
        return([self.points[0], self.points[1], self.points[5], self.points[4]])

    def backFace(self) -> list[Point]:
        """
        Points [3, 7, 6, 2]
        """
        # return([self.points[3], self.points[2], self.points[6], self.points[7]])
        return([self.points[3], self.points[7], self.points[6], self.points[2]])

    def topFace(self) -> list[Point]:
        """
        Points [4, 5, 6, 7]
        """
        return([self.points[4], self.points[5], self.points[6], self.points[7]])

    def bottomFace(self) -> list[Point]:
        """
        Points [0, 3, 2, 1]
        """
        return([self.points[0], self.points[3], self.points[2], self.points[1]])

    def leftFace(self) -> list[Point]:
        """
        Points [0, 4, 7, 3]
        """
        return([self.points[0], self.points[4], self.points[7], self.points[3]])

    def rightFace(self) -> list[Point]:
        """
        Points [1, 2, 6, 5]
        """
        return([self.points[1], self.points[2], self.points[6], self.points[5]])

    def getFace(self, faceName: str) -> list[Point]:
        """
        Return the face as a list of Point. Equivalent to :

            block.getFace('top') == block.topFace()
        """
        check_type("faceName", faceName, str)
        check_value("faceName", faceName, _FACE_NAME_TYPES)
        if (faceName == "top"):
            return(self.topFace())
        elif (faceName == "bottom"):
            return(self.bottomFace())
        elif (faceName == "left"):
            return(self.leftFace())
        elif (faceName == "right"):
            return(self.rightFace())
        elif (faceName == "front"):
            return(self.frontFace())
        elif (faceName == "back"):
            return(self.backFace())
        raise ValueError(f"{faceName} does not exists")

    def getOppositeFacename(self, faceName: str) -> str:
        """
        Return the face name at the opposite of the block.
        - top -> bottom
        - bottom -> top
        - left -> right
        - right -> left
        - front -> back
        - back -> front

        Equivalent to :

            block.getOppositeFace('top') == block.bottomFace()
        """
        check_type("faceName", faceName, str)
        check_value("faceName", faceName, _FACE_NAME_TYPES)
        if (faceName == "bottom"):
            return("top")
        elif (faceName == "top"):
            return("bottom")
        elif (faceName == "right"):
            return("left")
        elif (faceName == "left"):
            return("right")
        elif (faceName == "back"):
            return("front")
        elif (faceName == "front"):
            return("back")
        raise ValueError(f"{faceName} does not exists")

    def getOppositeFace(self, faceName: str) -> list[Point]:
        """
        Return the face at the opposite of the block as a list of Point.
        - top -> bottom
        - bottom -> top
        - left -> right
        - right -> left
        - front -> back
        - back -> front

        Equivalent to :

            block.getOppositeFace('top') == block.bottomFace()
        """
        check_type("faceName", faceName, str)
        check_value("faceName", faceName, _FACE_NAME_TYPES)
        if (faceName == "bottom"):
            return(self.topFace())
        elif (faceName == "top"):
            return(self.bottomFace())
        elif (faceName == "right"):
            return(self.leftFace())
        elif (faceName == "left"):
            return(self.rightFace())
        elif (faceName == "back"):
            return(self.frontFace())
        elif (faceName == "front"):
            return(self.backFace())
        raise ValueError(f"{faceName} does not exists")

    def getFaceBarycenter(self, face: list[Point]) -> Vector:
        """
        Compute the barycenter position of the face
        """
        nPoints = len(face)
        return(Vector(
            sum([p.x for p in face]) / nPoints,
            sum([p.y for p in face]) / nPoints,
            sum([p.z for p in face]) / nPoints
        ))

    def getFaceNormal(self, faceName: str) -> Vector:
        check_type("faceName", faceName, str)
        check_value("faceName", faceName, _FACE_NAME_TYPES)

        face = self.getFace(faceName)
        p1, p2, p3 = face[0], face[1], face[2]
        vec1 = Vector(p1.x-p2.x, p1.y-p2.y, p1.z-p2.z)
        vec2 = Vector(p2.x-p3.x, p2.y-p3.y, p2.z-p3.z)
        # if (faceName in ['front', 'top', 'bottom', 'left', 'right']):
        return(vec1.cross(vec2))
        return(vec2.cross(vec1))



class BlockMesh(OpenFOAMFile, Mesh):
    """
    Meshing object using the `blockMesh <https://www.openfoam.com/documentation/user-guide/4-mesh-generation-and-conversion/4.3-mesh-generation-with-the-blockmesh-utility>`_
    utility of OpenFOAM.

    Parameters
    ----------
    region : str
        Name of the region (default `""`).
    scale : float
        Scaling factor of the mesh (default 1).

    Attributes
    ----------
    scale : float
        Scaling factor of the mesh.
    blocks : list[Block]
        List of blocks
    faces : list[Face]
        List of faces
    mergePatchPairs : list[tuple[Face]]
        List of merged patch pairs faces. Used to stitch `Block`.
    baffleFaces : list[tuple[Face]]
        List of baffle patch pairs faces. Used to create baffle faces
    spheres : list[dict]
        List of spheres geometries usually used for point/edge/face projections
    cylinders : list[dict]
        List of cylinders geometries usually used for point/edge/face
        projections
    cones : list[dict]
        List of cones geometries usually used for point/edge/face projections
    pipeWallBC : Face
        Pipe wall boundary condition face named `pipeWall`.
    """
    def __init__(
            self,
            region: str="",
            scale: float=1
        ):
        OpenFOAMFile.__init__(self, name="blockMeshDict", folder="system", region=region)
        Mesh.__init__(self, region=region)

        self.scale = scale
        self.isReducedCells: bool = False
        self.blocks: list[Block] = []
        self.faces: list[Face] = []
        self.mergePatchPairs: list[tuple[Face]] = []
        self.baffleFaces: list[tuple[Face]] = []
        # List of basic geometry for face projection
        self.spheres: list[dict] = []
        self.cylinders: list[dict] = []
        self.cones: list[dict] = []

        self.pipeWallBC = Face("pipeWall", boundaryType="wall")
        self.addBoundary(self.pipeWallBC)


    @property
    def cellZones(self) -> list[str]:
        """
        Return the list of cellZone names.
        """
        zones = []
        for block in self.blocks:
            if (block.name not in zones):
                zones.append(block.name)
        return(zones)


    def addSphere(
            self,
            name: str,
            radius: float,
            origin: tuple=(0, 0, 0)
        ) -> None:
        """
        Create a sphere geometry shape that can be used for point, edge and
        face projection.
        """
        self.spheres.append({
            'name': name,
            'origin': origin,
            'radius': radius
        })


    def addCylinder(
            self,
            name: str,
            radius: float,
            point1: tuple=(0, 0, 0), point2: tuple=(0, 0, 0)
        ) -> None:
        """
        Create a cylinder geometry shape that can be used for point, edge and
        face projection.
        """
        self.cylinders.append({
            'name': name,
            'point1': point1,
            'point2': point2,
            'radius': radius
        })


    def addCone(
            self,
            name: str,
            radius1: float,
            radius2: float,
            point1: tuple=(0, 0, 0), point2: tuple=(0, 0, 0)
        ) -> None:
        """
        Create a cone geometry shape that can be used for point, edge and
        face projection.
        """
        self.cones.append({
            'name': name,
            'point1': point1,
            'point2': point2,
            'radius1': radius1,
            'radius2': radius2
        })


    def createBlock(
            self,
            name: str,
            points: list[Point],
            nx: int, ny: int, nz: int,
            gradx: float=1, grady: float=1, gradz: float=1
        ) -> Block:
        """
        Create a block using 8 explicit Point.

        Return
        ------
            (newBlock)
        """
        check_type("name", name, str)
        check_type("points", points, list)
        check_type("nx", nx, int)
        check_type("ny", ny, int)
        check_type("nz", nz, int)
        check_type("gradx", gradx, (float, int))
        check_type("grady", grady, (float, int))
        check_type("gradz", gradz, (float, int))
        newBlock = Block(name, points, nx, ny, nz, gradx=gradx, grady=grady, gradz=gradz)
        self.blocks.append(newBlock)
        return(newBlock)


    def createWedge(
            self,
            name: str,
            innerRadius: float, outerRadius: float, lowZ: float, highZ: float,
            wedgeAngle: float,
            nr: int, nz: int,
            gradr: float=1, gradz: float=1
        ) -> Block:
        """
        Parameters
        ----------
        wedgeAngle : float
            Total wedge aperture angle in degree

        Return
        ------
            (newBlock)
        """
        deg = np.pi/180
        irX = innerRadius * np.cos(wedgeAngle/2 * deg)
        irY = innerRadius * np.sin(wedgeAngle/2 * deg)
        orX = outerRadius * np.cos(wedgeAngle/2 * deg)
        orY = outerRadius * np.sin(wedgeAngle/2 * deg)

        newBlock = self.createBlock(name, [
            Point(irX, -irY, lowZ),
            Point(orX, -orY, lowZ),
            Point(orX, orY, lowZ),
            Point(irX, irY, lowZ),
            Point(irX, -irY, highZ),
            Point(orX, -orY, highZ),
            Point(orX, orY, highZ),
            Point(irX, irY, highZ),
        ], nr, 1, nz, gradx=gradr, gradz=gradz)

        return(newBlock)


    def createWedgeConical(
            self,
            name: str,
            innerRadiusBottom: float,
            outerRadiusBottom: float,
            innerRadiusTop: float,
            outerRadiusTop: float,
            lowZ: float,
            highZ: float,
            wedgeAngle: float,
            nr: int, nz: int,
            gradr: float=1, gradz: float=1
        ) -> Block:
        """
        Create a conical wedge along the Z-axis.

        Parameters
        ----------
        wedgeAngle : float
            Total wedge aperture angle in degree

        Return
        ------
            (newBlock)
        """
        deg = np.pi/180
        irBotX = innerRadiusBottom * np.cos(wedgeAngle/2 * deg)
        irBotY = innerRadiusBottom * np.sin(wedgeAngle/2 * deg)
        orBotX = outerRadiusBottom * np.cos(wedgeAngle/2 * deg)
        orBotY = outerRadiusBottom * np.sin(wedgeAngle/2 * deg)
        irTopX = innerRadiusTop * np.cos(wedgeAngle/2 * deg)
        irTopY = innerRadiusTop * np.sin(wedgeAngle/2 * deg)
        orTopX = outerRadiusTop * np.cos(wedgeAngle/2 * deg)
        orTopY = outerRadiusTop * np.sin(wedgeAngle/2 * deg)

        newBlock = self.createBlock(name, [
            Point(irBotX, -irBotY, lowZ),
            Point(orBotX, -orBotY, lowZ),
            Point(orBotX, orBotY, lowZ),
            Point(irBotX, irBotY, lowZ),
            Point(irTopX, -irTopY, highZ),
            Point(orTopX, -orTopY, highZ),
            Point(orTopX, orTopY, highZ),
            Point(irTopX, irTopY, highZ),
        ], nr, 1, nz, gradx=gradr, gradz=gradz)

        return(newBlock)


    def createCube(
            self, name: str,
            lowX: float, lowY: float, lowZ: float,
            highX: float, highY: float, highZ: float,
            nx: int=1, ny: int=1, nz: int=1,
            gradx: float=1, grady: float=1, gradz: float=1,
            isAddBoundaryConditions: bool=False
        ) -> Block:
        """
        Return
        ------
            (newBlock)
        """
        newBlock = self.createBlock(name, [
            Point(lowX, lowY, lowZ),
            Point(highX, lowY, lowZ),
            Point(highX, highY, lowZ),
            Point(lowX, highY, lowZ),
            Point(lowX, lowY, highZ),
            Point(highX, lowY, highZ),
            Point(highX, highY, highZ),
            Point(lowX, highY, highZ),
        ], nx, ny, nz, gradx=gradx, grady=grady, gradz=gradz)

        if (isAddBoundaryConditions):
            idxFace = len(self.faces)

            topFace = Face(f"{name}Top_{idxFace}")
            topFace.addSubFace(newBlock.topFace())
            botFace = Face(f"{name}Bottom_{idxFace}")
            botFace.addSubFace(newBlock.bottomFace())

            wallFront = Face(f"{name}WallFront_{idxFace}", boundaryType='wall')
            wallFront.addSubFace(newBlock.frontFace())
            wallLeft = Face(f"{name}WallLeft_{idxFace}", boundaryType='wall')
            wallLeft.addSubFace(newBlock.leftFace())
            wallRight = Face(f"{name}WallRight_{idxFace}", boundaryType='wall')
            wallRight.addSubFace(newBlock.rightFace())
            wallBack = Face(f"{name}WallBack_{idxFace}", boundaryType='wall')
            wallBack.addSubFace(newBlock.backFace())

            self.addBoundary(topFace)
            self.addBoundary(botFace)
            self.addBoundary(wallFront)
            self.addBoundary(wallLeft)
            self.addBoundary(wallRight)
            self.addBoundary(wallBack)

        return(newBlock)


    def createSphere(
            self,
            name: str,
            radius: float,
            x: float=0,
            y: float=0,
            z: float=0,
            nCenter: int=1,
            nBorder: int=1,
            gradCenter: float=1,
            gradBorder: float=1,
            isAddBoundaryConditions: bool=False
        ):
        """
        Return
        ------
            (ballCenter, ballTop, ballRight, ballBack, ballLeft, ballFront, ballBottom)
        """
        sqrt3 = np.sqrt(3)

        self.addSphere(name, radius, (x, y, z))

        ballCenter = self.createCube(
            name,
            x-radius/3, y-radius/3, z-radius/3,
            x+radius/3, y+radius/3, z+radius/3,
            nCenter, nCenter, nCenter,
            gradx=gradCenter, grady=gradCenter, gradz=gradCenter
        )

        ballRight = self.addRight(ballCenter, name, [
            Point(x+radius/sqrt3, y-radius/sqrt3, z-radius/sqrt3),
            Point(x+radius/sqrt3, y+radius/sqrt3, z-radius/sqrt3),
            Point(x+radius/sqrt3, y-radius/sqrt3, z+radius/sqrt3),
            Point(x+radius/sqrt3, y+radius/sqrt3, z+radius/sqrt3),
        ], nx=nBorder, gradx=gradBorder)

        ballRight.addEdgeArc(1, 5, x=x, y=y, z=z, isOrigin=True)
        ballRight.addEdgeArc(5, 6, x=x, y=y, z=z, isOrigin=True)
        ballRight.addEdgeArc(6, 2, x=x, y=y, z=z, isOrigin=True)
        ballRight.addEdgeArc(1, 2, x=x, y=y, z=z, isOrigin=True)

        ballRight.addFaceProjection("right", name)

        ballLeft = self.addLeft(ballCenter, name, [
            Point(x-radius/sqrt3, y-radius/sqrt3, z-radius/sqrt3),
            Point(x-radius/sqrt3, y+radius/sqrt3, z-radius/sqrt3),
            Point(x-radius/sqrt3, y-radius/sqrt3, z+radius/sqrt3),
            Point(x-radius/sqrt3, y+radius/sqrt3, z+radius/sqrt3),
        ], ballRight.nx, gradx=1/gradBorder)

        ballLeft.addEdgeArc(0, 4, x=x, y=y, z=z, isOrigin=True)
        ballLeft.addEdgeArc(4, 7, x=x, y=y, z=z, isOrigin=True)
        ballLeft.addEdgeArc(7, 3, x=x, y=y, z=z, isOrigin=True)
        ballLeft.addEdgeArc(0, 3, x=x, y=y, z=z, isOrigin=True)

        ballLeft.addFaceProjection("left", name)

        ballFront = self.addFront(ballCenter, name, [
            ballLeft.points[0],
            ballRight.points[1],
            ballLeft.points[4],
            ballRight.points[5],
        ], ballRight.nx, grady=1/gradBorder)

        ballFront.addEdgeArc(4, 5, x=x, y=y, z=z, isOrigin=True)
        ballFront.addEdgeArc(0, 1, x=x, y=y, z=z, isOrigin=True)

        ballFront.addFaceProjection("front", name)

        ballBack = self.addBack(ballCenter, name, [
            ballRight.points[2],
            ballLeft.points[3],
            ballRight.points[6],
            ballLeft.points[7],
        ], ballRight.nx, grady=gradBorder)

        ballBack.addEdgeArc(6, 7, x=x, y=y, z=z, isOrigin=True)
        ballBack.addEdgeArc(3, 2, x=x, y=y, z=z, isOrigin=True)

        ballBack.addFaceProjection("back", name)

        ballTop = self.addTop(ballCenter, name, [
            ballLeft.points[4],
            ballRight.points[5],
            ballRight.points[6],
            ballLeft.points[7],
        ], ballRight.nx, gradz=gradBorder)

        ballTop.addFaceProjection("top", name)

        ballBottom = self.addBottom(ballCenter, name, [
            ballLeft.points[0],
            ballRight.points[1],
            ballRight.points[2],
            ballLeft.points[3],
        ], ballRight.nx, gradz=1/gradBorder)

        ballBottom.addFaceProjection("bottom", name)


        if (isAddBoundaryConditions):
            idxFace = len(self.faces)

            wallFaceFront = Face(f"{name}WallFront_{idxFace}", boundaryType='wall')
            wallFaceFront.addSubFace(ballFront.frontFace())
            wallFaceBack = Face(f"{name}WallBack_{idxFace}", boundaryType='wall')
            wallFaceBack.addSubFace(ballBack.backFace())
            wallFaceRight = Face(f"{name}WallRight_{idxFace}", boundaryType='wall')
            wallFaceRight.addSubFace(ballRight.rightFace())
            wallFaceLeft = Face(f"{name}WallLeft_{idxFace}", boundaryType='wall')
            wallFaceLeft.addSubFace(ballLeft.leftFace())
            wallFaceTop = Face(f"{name}WallTop_{idxFace}", boundaryType='wall')
            wallFaceTop.addSubFace(ballTop.topFace())
            wallFaceBottom = Face(f"{name}WallBottom_{idxFace}", boundaryType='wall')
            wallFaceBottom.addSubFace(ballBottom.bottomFace())

            self.addBoundary(wallFaceFront)
            self.addBoundary(wallFaceBack)
            self.addBoundary(wallFaceRight)
            self.addBoundary(wallFaceLeft)
            self.addBoundary(wallFaceTop)
            self.addBoundary(wallFaceBottom)

        return(ballCenter, ballTop, ballRight, ballBack, ballLeft, ballFront, ballBottom)


    def createHalfSphere(
            self,
            name: str,
            radius: float,
            x: float=0,
            y: float=0,
            z: float=0,
            nCenter: int=1,
            nBorder: int=1,
            gradCenter: float=1,
            gradBorder: float=1,
            isAddBoundaryConditions: bool=False
        ):
        """
        Return
        ------
            (ballCenter, ballRight, ballBack, ballLeft, ballFront, ballBottom)
        """
        sqrt2 = np.sqrt(2)
        sqrt3 = np.sqrt(3)

        self.addSphere(name, radius, (x, y, z))

        ballCenter = self.createCube(
            name,
            x-radius/3, y-radius/3, z-radius/3,
            x+radius/3, y+radius/3, z,
            nCenter, nCenter, nCenter,
            gradx=gradCenter, grady=gradCenter, gradz=gradCenter
        )

        ballRight = self.addRight(ballCenter, name, [
            Point(x+radius/sqrt3, y-radius/sqrt3, z-radius/sqrt3),
            Point(x+radius/sqrt3, y+radius/sqrt3, z-radius/sqrt3),
            Point(x+radius/sqrt2, y-radius/sqrt2, z),
            Point(x+radius/sqrt2, y+radius/sqrt2, z),
        ], nx=nBorder, gradx=gradBorder)

        ballRight.addEdgeArc(1, 5, x=x, y=y, z=z, isOrigin=True)
        ballRight.addEdgeArc(5, 6, x=x, y=y, z=z, isOrigin=True)
        ballRight.addEdgeArc(6, 2, x=x, y=y, z=z, isOrigin=True)
        ballRight.addEdgeArc(1, 2, x=x, y=y, z=z, isOrigin=True)

        ballRight.addFaceProjection("right", name)

        ballLeft = self.addLeft(ballCenter, name, [
            Point(x-radius/sqrt3, y-radius/sqrt3, z-radius/sqrt3),
            Point(x-radius/sqrt3, y+radius/sqrt3, z-radius/sqrt3),
            Point(x-radius/sqrt2, y-radius/sqrt2, z),
            Point(x-radius/sqrt2, y+radius/sqrt2, z),
        ], ballRight.nx, gradx=1/gradBorder)

        ballLeft.addEdgeArc(0, 4, x=x, y=y, z=z, isOrigin=True)
        ballLeft.addEdgeArc(4, 7, x=x, y=y, z=z, isOrigin=True)
        ballLeft.addEdgeArc(7, 3, x=x, y=y, z=z, isOrigin=True)
        ballLeft.addEdgeArc(0, 3, x=x, y=y, z=z, isOrigin=True)

        ballLeft.addFaceProjection("left", name)

        ballFront = self.addFront(ballCenter, name, [
            ballLeft.points[0],
            ballRight.points[1],
            ballLeft.points[4],
            ballRight.points[5],
        ], ballRight.nx, grady=1/gradBorder)

        ballFront.addEdgeArc(4, 5, x=x, y=y, z=z, isOrigin=True)
        ballFront.addEdgeArc(0, 1, x=x, y=y, z=z, isOrigin=True)

        ballFront.addFaceProjection("front", name)

        ballBack = self.addBack(ballCenter, name, [
            ballRight.points[2],
            ballLeft.points[3],
            ballRight.points[6],
            ballLeft.points[7],
        ], ballRight.nx, grady=gradBorder)

        ballBack.addEdgeArc(6, 7, x=x, y=y, z=z, isOrigin=True)
        ballBack.addEdgeArc(3, 2, x=x, y=y, z=z, isOrigin=True)

        ballBack.addFaceProjection("back", name)

        ballBottom = self.addBottom(ballCenter, name, [
            ballLeft.points[0],
            ballRight.points[1],
            ballRight.points[2],
            ballLeft.points[3],
        ], ballRight.nx, gradz=1/gradBorder)

        ballBottom.addFaceProjection("bottom", name)


        if (isAddBoundaryConditions):
            idxFace = len(self.faces)

            wallFaceTop = Face(f"{name}WallTop_{idxFace}", boundaryType='wall')
            wallFaceTop.addSubFace(ballCenter.topFace())
            wallFaceTop.addSubFace(ballBack.topFace())
            wallFaceTop.addSubFace(ballFront.topFace())
            wallFaceTop.addSubFace(ballRight.topFace())
            wallFaceTop.addSubFace(ballLeft.topFace())
            wallFaceBottom = Face(f"{name}WallBottom_{idxFace}", boundaryType='wall')
            wallFaceBottom.addSubFace(ballBottom.bottomFace())
            wallFaceFront = Face(f"{name}WallFront_{idxFace}", boundaryType='wall')
            wallFaceFront.addSubFace(ballFront.frontFace())
            wallFaceBack = Face(f"{name}WallBack_{idxFace}", boundaryType='wall')
            wallFaceBack.addSubFace(ballBack.backFace())
            wallFaceRight = Face(f"{name}WallRight_{idxFace}", boundaryType='wall')
            wallFaceRight.addSubFace(ballRight.rightFace())
            wallFaceLeft = Face(f"{name}WallLeft_{idxFace}", boundaryType='wall')
            wallFaceLeft.addSubFace(ballLeft.leftFace())

            self.addBoundary(wallFaceTop)
            self.addBoundary(wallFaceBottom)
            self.addBoundary(wallFaceFront)
            self.addBoundary(wallFaceBack)
            self.addBoundary(wallFaceRight)
            self.addBoundary(wallFaceLeft)

        return(ballCenter, ballRight, ballBack, ballLeft, ballFront, ballBottom)


    def createHollowHalfSphere(
            self,
            name: str,
            innerRadius: float,
            outerRadius: float,
            x: float=0,
            y: float=0,
            z: float=0,
            nr: int=1,
            nt: int=1,
            isAddBoundaryConditions: bool=False
        ):
        """
        Return
        ------
            (ballRight, ballBack, ballLeft, ballFront, ballBottom)
        """
        sqrt2 = np.sqrt(2)
        sqrt3 = np.sqrt(3)

        self.addSphere(name+"_inner", innerRadius, (x, y, z))
        self.addSphere(name+"_outer", outerRadius, (x, y, z))

        ballFront = self.createBlock(name, [
            Point(x-outerRadius/sqrt3, y-outerRadius/sqrt3, z-outerRadius/sqrt3),
            Point(x+outerRadius/sqrt3, y-outerRadius/sqrt3, z-outerRadius/sqrt3),
            Point(x+innerRadius/sqrt3, y-innerRadius/sqrt3, z-innerRadius/sqrt3),
            Point(x-innerRadius/sqrt3, y-innerRadius/sqrt3, z-innerRadius/sqrt3),
            Point(x-outerRadius/sqrt2, y-outerRadius/sqrt2, z),
            Point(x+outerRadius/sqrt2, y-outerRadius/sqrt2, z),
            Point(x+innerRadius/sqrt2, y-innerRadius/sqrt2, z),
            Point(x-innerRadius/sqrt2, y-innerRadius/sqrt2, z),
        ], nx=nt, ny=nr, nz=nt)

        ballFront.addEdgeArc(0, 1, x=x, y=y, z=z, isOrigin=True)
        ballFront.addEdgeArc(2, 3, x=x, y=y, z=z, isOrigin=True)
        ballFront.addEdgeArc(4, 5, x=x, y=y, z=z, isOrigin=True)
        ballFront.addEdgeArc(6, 7, x=x, y=y, z=z, isOrigin=True)
        ballFront.addEdgeArc(0, 4, x=x, y=y, z=z, isOrigin=True)
        ballFront.addEdgeArc(1, 5, x=x, y=y, z=z, isOrigin=True)
        ballFront.addEdgeArc(2, 6, x=x, y=y, z=z, isOrigin=True)
        ballFront.addEdgeArc(3, 7, x=x, y=y, z=z, isOrigin=True)

        ballFront.addFaceProjection("back", name+"_inner")
        ballFront.addFaceProjection("front", name+"_outer")


        ballBack = self.createBlock(name, [
            Point(x-innerRadius/sqrt3, y+innerRadius/sqrt3, z-innerRadius/sqrt3),
            Point(x+innerRadius/sqrt3, y+innerRadius/sqrt3, z-innerRadius/sqrt3),
            Point(x+outerRadius/sqrt3, y+outerRadius/sqrt3, z-outerRadius/sqrt3),
            Point(x-outerRadius/sqrt3, y+outerRadius/sqrt3, z-outerRadius/sqrt3),
            Point(x-innerRadius/sqrt2, y+innerRadius/sqrt2, z),
            Point(x+innerRadius/sqrt2, y+innerRadius/sqrt2, z),
            Point(x+outerRadius/sqrt2, y+outerRadius/sqrt2, z),
            Point(x-outerRadius/sqrt2, y+outerRadius/sqrt2, z),
        ], nx=nt, ny=nr, nz=nt)

        ballBack.addEdgeArc(0, 1, x=x, y=y, z=z, isOrigin=True)
        ballBack.addEdgeArc(2, 3, x=x, y=y, z=z, isOrigin=True)
        ballBack.addEdgeArc(4, 5, x=x, y=y, z=z, isOrigin=True)
        ballBack.addEdgeArc(6, 7, x=x, y=y, z=z, isOrigin=True)
        ballBack.addEdgeArc(0, 4, x=x, y=y, z=z, isOrigin=True)
        ballBack.addEdgeArc(1, 5, x=x, y=y, z=z, isOrigin=True)
        ballBack.addEdgeArc(2, 6, x=x, y=y, z=z, isOrigin=True)
        ballBack.addEdgeArc(3, 7, x=x, y=y, z=z, isOrigin=True)

        ballBack.addFaceProjection("front", name+"_inner")
        ballBack.addFaceProjection("back", name+"_outer")


        ballRight = self.createBlock(name, [
            ballFront.points[2],
            ballFront.points[1],
            ballBack.points[2],
            ballBack.points[1],
            ballFront.points[6],
            ballFront.points[5],
            ballBack.points[6],
            ballBack.points[5],
        ], nx=nr, ny=nt, nz=nt)

        ballRight.addEdgeArc(0, 3, x=x, y=y, z=z, isOrigin=True)
        ballRight.addEdgeArc(1, 2, x=x, y=y, z=z, isOrigin=True)
        ballRight.addEdgeArc(4, 7, x=x, y=y, z=z, isOrigin=True)
        ballRight.addEdgeArc(5, 6, x=x, y=y, z=z, isOrigin=True)

        ballRight.addFaceProjection("left", name+"_inner")
        ballRight.addFaceProjection("right", name+"_outer")

        ballLeft = self.createBlock(name, [
            ballFront.points[0],
            ballFront.points[3],
            ballBack.points[0],
            ballBack.points[3],
            ballFront.points[4],
            ballFront.points[7],
            ballBack.points[4],
            ballBack.points[7],
        ], nx=nr, ny=nt, nz=nt)

        ballLeft.addEdgeArc(0, 3, x=x, y=y, z=z, isOrigin=True)
        ballLeft.addEdgeArc(1, 2, x=x, y=y, z=z, isOrigin=True)
        ballLeft.addEdgeArc(4, 7, x=x, y=y, z=z, isOrigin=True)
        ballLeft.addEdgeArc(5, 6, x=x, y=y, z=z, isOrigin=True)

        ballLeft.addFaceProjection("right", name+"_inner")
        ballLeft.addFaceProjection("left", name+"_outer")

        ballBottom = self.createBlock(name, [
            ballFront.points[0],
            ballFront.points[1],
            ballBack.points[2],
            ballBack.points[3],
            ballFront.points[3],
            ballFront.points[2],
            ballBack.points[1],
            ballBack.points[0],
        ], nx=nt, ny=nt, nz=nr)

        ballBottom.addFaceProjection("top", name+"_inner")
        ballBottom.addFaceProjection("bottom", name+"_outer")


        if (isAddBoundaryConditions):
            idxFace = len(self.faces)

            wallFaceTop = Face(f"{name}WallTop_{idxFace}", boundaryType='wall')
            wallFaceTop.addSubFace(ballBack.topFace())
            wallFaceTop.addSubFace(ballFront.topFace())
            wallFaceTop.addSubFace(ballRight.topFace())
            wallFaceTop.addSubFace(ballLeft.topFace())
            outerWallBack = Face(f"{name}OuterWallBack_{idxFace}", boundaryType='wall')
            outerWallBack.addSubFace(ballBack.backFace())
            outerWallFront = Face(f"{name}OuterWallFront_{idxFace}", boundaryType='wall')
            outerWallFront.addSubFace(ballFront.frontFace())
            outerWallLeft = Face(f"{name}OuterWallLeft_{idxFace}", boundaryType='wall')
            outerWallLeft.addSubFace(ballLeft.leftFace())
            outerWallRight = Face(f"{name}OuterWallRight_{idxFace}", boundaryType='wall')
            outerWallRight.addSubFace(ballRight.rightFace())
            outerWallBottom = Face(f"{name}OuterWallBottom_{idxFace}", boundaryType='wall')
            outerWallBottom.addSubFace(ballBottom.bottomFace())
            innerWallBack = Face(f"{name}InnerWallBack_{idxFace}", boundaryType='wall')
            innerWallBack.addSubFace(ballFront.backFace())
            innerWallFront = Face(f"{name}InnerWallFront_{idxFace}", boundaryType='wall')
            innerWallFront.addSubFace(ballBack.frontFace())
            innerWallLeft = Face(f"{name}InnerWallLeft_{idxFace}", boundaryType='wall')
            innerWallLeft.addSubFace(ballRight.leftFace())
            innerWallRight = Face(f"{name}InnerWallRight_{idxFace}", boundaryType='wall')
            innerWallRight.addSubFace(ballLeft.rightFace())
            innerWallBottom = Face(f"{name}InnerWallBottom_{idxFace}", boundaryType='wall')
            innerWallBottom.addSubFace(ballBottom.topFace())

            self.addBoundary(wallFaceTop)
            self.addBoundary(outerWallBack)
            self.addBoundary(outerWallFront)
            self.addBoundary(outerWallLeft)
            self.addBoundary(outerWallRight)
            self.addBoundary(outerWallBottom)
            self.addBoundary(innerWallBack)
            self.addBoundary(innerWallFront)
            self.addBoundary(innerWallLeft)
            self.addBoundary(innerWallRight)
            self.addBoundary(innerWallBottom)

        return(ballRight, ballBack, ballLeft, ballFront, ballBottom)


    def createSphere1D(
            self,
            name: str,
            innerRadius: float,
            outerRadius: float,
            opening: float,
            x: float=0,
            y: float=0,
            z: float=0,
            nr: int=1,
            gradr: float=1,
            isAddBoundaryConditions: bool=False
        ):
        """
        Create a 1D sphere along the X-axis.

        Parameters
        ----------
        name : str
            Name of the block.
        innerRadius : float
            Inner radius of the sphere.
        outerRadius : float
            Outer radius of the sphere.
        opening : float
            Opening angle of the sphere in deg.
        x : float
            X-position (default `0`).
        y : float=0
            Y-position (default `0`).
        z : float=0
            Z-position (default `0`).
        nr : int
            Number of radial cells (default `1`).
        gradr : float
            Cells radial grading (default `1`).
        isAddBoundaryConditions : bool
            If True, create new boundary conditions (default `False`).

        Return
        ------
            block
        """
        deg = np.pi/180

        tant = np.tan(opening/2 * deg)

        riX = innerRadius / np.sqrt(1 + 2*tant**2)
        roX = outerRadius / np.sqrt(1 + 2*tant**2)
        riYZ = riX * tant
        roYZ = roX * tant

        block = self.createBlock(
            name=name,
            points=[
                Point(x + riX, y - riYZ, z - riYZ),
                Point(x + roX, y - roYZ, z - roYZ),
                Point(x + roX, y + roYZ, z - roYZ),
                Point(x + riX, y + riYZ, z - riYZ),
                Point(x + riX, y - riYZ, z + riYZ),
                Point(x + roX, y - roYZ, z + roYZ),
                Point(x + roX, y + roYZ, z + roYZ),
                Point(x + riX, y + riYZ, z + riYZ),
            ],
            nx=nr, ny=1, nz=1,
            gradx=gradr, grady=1, gradz=1
        )

        if (isAddBoundaryConditions):
            idxFace = len(self.faces)

            wedgeFaceFront = Face(f"{name}WedgeFront_{idxFace}", boundaryType='wedge')
            wedgeFaceFront.addSubFace(block.frontFace())
            wedgeFaceBack = Face(f"{name}WedgeBack_{idxFace}", boundaryType='wedge')
            wedgeFaceBack.addSubFace(block.backFace())
            wedgeFaceTop = Face(f"{name}WedgeTop_{idxFace}", boundaryType='wedge')
            wedgeFaceTop.addSubFace(block.topFace())
            wedgeFaceBottom = Face(f"{name}WedgeBottom_{idxFace}", boundaryType='wedge')
            wedgeFaceBottom.addSubFace(block.bottomFace())
            wallFaceOuter = Face(f"{name}WallOuter_{idxFace}", boundaryType='patch')
            wallFaceOuter.addSubFace(block.rightFace())
            wallFaceInner = Face(f"{name}WallInner_{idxFace}", boundaryType='patch')
            wallFaceInner.addSubFace(block.leftFace())

            self.addBoundary(wedgeFaceFront)
            self.addBoundary(wedgeFaceBack)
            self.addBoundary(wedgeFaceTop)
            self.addBoundary(wedgeFaceBottom)
            self.addBoundary(wallFaceOuter)
            self.addBoundary(wallFaceInner)

        return(block)


    def createRingAlongZ(
            self,
            name: str,
            innerRadius: float, outerRadius: float, lowZ: float, highZ: float,
            nr: int, nt: int, nz: int,
            x: float=0, y: float=0,
            isAddBoundaryConditions: bool=False
        ):
        """
        Return
        ------
            (frontBlock, rightBlock, backBlock, leftBlock)
        """
        sqrt2 = np.sqrt(2)

        frontBlock = self.createBlock(name, [
            Point(x-outerRadius/sqrt2, y-outerRadius/sqrt2, lowZ),
            Point(x+outerRadius/sqrt2, y-outerRadius/sqrt2, lowZ),
            Point(x+innerRadius/sqrt2, y-innerRadius/sqrt2, lowZ),
            Point(x-innerRadius/sqrt2, y-innerRadius/sqrt2, lowZ),
            Point(x-outerRadius/sqrt2, y-outerRadius/sqrt2, highZ),
            Point(x+outerRadius/sqrt2, y-outerRadius/sqrt2, highZ),
            Point(x+innerRadius/sqrt2, y-innerRadius/sqrt2, highZ),
            Point(x-innerRadius/sqrt2, y-innerRadius/sqrt2, highZ)
        ], nt, nr, nz)

        rightBlock = self.createBlock(name, [
            frontBlock.points[2],
            frontBlock.points[1],
            Point(x+outerRadius/sqrt2, y+outerRadius/sqrt2, lowZ),
            Point(x+innerRadius/sqrt2, y+innerRadius/sqrt2, lowZ),
            frontBlock.points[6],
            frontBlock.points[5],
            Point(x+outerRadius/sqrt2, y+outerRadius/sqrt2, highZ),
            Point(x+innerRadius/sqrt2, y+innerRadius/sqrt2, highZ)
        ], nr, nt, nz)

        leftBlock = self.createBlock(name, [
            frontBlock.points[0],
            frontBlock.points[3],
            Point(x-innerRadius/sqrt2, y+innerRadius/sqrt2, lowZ),
            Point(x-outerRadius/sqrt2, y+outerRadius/sqrt2, lowZ),
            frontBlock.points[4],
            frontBlock.points[7],
            Point(x-innerRadius/sqrt2, y+innerRadius/sqrt2, highZ),
            Point(x-outerRadius/sqrt2, y+outerRadius/sqrt2, highZ)
        ], nr, nt, nz)

        backBlock = self.createBlock(name, [
            leftBlock.points[2],
            rightBlock.points[3],
            rightBlock.points[2],
            leftBlock.points[3],
            leftBlock.points[6],
            rightBlock.points[7],
            rightBlock.points[6],
            leftBlock.points[7],
        ], nt, nr, nz)

        frontBlock.addEdgeArc(0, 1, x=x, y=y, isOrigin=True)
        frontBlock.addEdgeArc(2, 3, x=x, y=y, isOrigin=True)
        frontBlock.addEdgeArc(4, 5, x=x, y=y, isOrigin=True)
        frontBlock.addEdgeArc(6, 7, x=x, y=y, isOrigin=True)
        rightBlock.addEdgeArc(1, 2, x=x, y=y, isOrigin=True)
        rightBlock.addEdgeArc(0, 3, x=x, y=y, isOrigin=True)
        rightBlock.addEdgeArc(5, 6, x=x, y=y, isOrigin=True)
        rightBlock.addEdgeArc(4, 7, x=x, y=y, isOrigin=True)
        leftBlock.addEdgeArc(1, 2, x=x, y=y, isOrigin=True)
        leftBlock.addEdgeArc(3, 0, x=x, y=y, isOrigin=True)
        leftBlock.addEdgeArc(5, 6, x=x, y=y, isOrigin=True)
        leftBlock.addEdgeArc(7, 4, x=x, y=y, isOrigin=True)
        backBlock.addEdgeArc(0, 1, x=x, y=y, isOrigin=True)
        backBlock.addEdgeArc(2, 3, x=x, y=y, isOrigin=True)
        backBlock.addEdgeArc(4, 5, x=x, y=y, isOrigin=True)
        backBlock.addEdgeArc(6, 7, x=x, y=y, isOrigin=True)

        if (isAddBoundaryConditions):
            idxFace = len(self.faces)

            topFace = Face(f"{name}Top_{idxFace}")
            botFace = Face(f"{name}Bottom_{idxFace}")

            for block in [frontBlock, rightBlock, backBlock, leftBlock]:
                topFace.addSubFace(block.topFace())
                botFace.addSubFace(block.bottomFace())

            innerWallFaceFront = Face(f"{name}InnerWallFront_{idxFace}", boundaryType='wall')
            innerWallFaceFront.addSubFace(frontBlock.backFace())
            innerWallFaceLeft = Face(f"{name}InnerWallLeft_{idxFace}", boundaryType='wall')
            innerWallFaceLeft.addSubFace(leftBlock.rightFace())
            innerWallFaceRight = Face(f"{name}InnerWallRight_{idxFace}", boundaryType='wall')
            innerWallFaceRight.addSubFace(rightBlock.leftFace())
            innerWallFaceBack = Face(f"{name}InnerWallBack_{idxFace}", boundaryType='wall')
            innerWallFaceBack.addSubFace(backBlock.frontFace())

            outerWallFaceFront = Face(f"{name}OuterWallFront_{idxFace}", boundaryType='wall')
            outerWallFaceFront.addSubFace(frontBlock.frontFace())
            outerWallFaceLeft = Face(f"{name}OuterWallLeft_{idxFace}", boundaryType='wall')
            outerWallFaceLeft.addSubFace(leftBlock.leftFace())
            outerWallFaceRight = Face(f"{name}OuterWallRight_{idxFace}", boundaryType='wall')
            outerWallFaceRight.addSubFace(rightBlock.rightFace())
            outerWallFaceBack = Face(f"{name}OuterWallBack_{idxFace}", boundaryType='wall')
            outerWallFaceBack.addSubFace(backBlock.backFace())

            self.addBoundary(topFace)
            self.addBoundary(botFace)
            self.addBoundary(innerWallFaceFront)
            self.addBoundary(innerWallFaceLeft)
            self.addBoundary(innerWallFaceRight)
            self.addBoundary(innerWallFaceBack)
            self.addBoundary(outerWallFaceFront)
            self.addBoundary(outerWallFaceLeft)
            self.addBoundary(outerWallFaceRight)
            self.addBoundary(outerWallFaceBack)

        return(frontBlock, rightBlock, backBlock, leftBlock)


    def createRingSectorAlongZ(
            self,
            name: str,
            innerRadius: float,
            outerRadius: float,
            angleStart: float,
            angleArc: float,
            lowZ: float,
            highZ: float,
            x: float=0, y: float=0,
            nr: int=1,
            nt: int=1,
            nz: int=1,
            isAddBoundaryConditions: bool=False
        ) -> Block:
        """
        Face orientation:
        - `front`: large arc face
        - `back`: short arc face

        Parameters
        ----------
        name : str
            Name of the block
        innerRadius : float
            Inner radius
        outerRadius : float
            Outer radius
        angleStart : float
            Starting angle in degree
        angleArc : float
            Spanning angle of the arc drawn by the ring in degree
        lowZ : float
            Minimum Z-position
        highZ : float
            Maximum Z-position
        x : float
            Offset position in along the X-axis (default `0`)
        y : float
            Offset position in along the Y-axis (default `0`)
        nr : int
            Number of radial discretization (default `1`)
        nt : int
            Number of azimutal discretization (default `1`)
        nz : int
            Number of axial discretization (default `1`)

        Return
        ------
            (newBlock)
        """
        # Convert from deg to rad
        deg = np.pi/180
        angleStart = angleStart * deg
        angleArc = angleArc * deg

        cost, sint = np.cos(angleStart), np.sin(angleStart)
        costt, sintt = np.cos(angleStart+angleArc), np.sin(angleStart+angleArc)

        newBlock = self.createBlock(name, [
            Point(x+outerRadius*cost, y+outerRadius*sint, lowZ),
            Point(x+outerRadius*costt, y+outerRadius*sintt, lowZ),
            Point(x+innerRadius*costt, y+innerRadius*sintt, lowZ),
            Point(x+innerRadius*cost, y+innerRadius*sint, lowZ),
            Point(x+outerRadius*cost, y+outerRadius*sint, highZ),
            Point(x+outerRadius*costt, y+outerRadius*sintt, highZ),
            Point(x+innerRadius*costt, y+innerRadius*sintt, highZ),
            Point(x+innerRadius*cost, y+innerRadius*sint, highZ),
        ], nx=nt, ny=nr, nz=nz)

        newBlock.addEdgeArc(0, 1, x, y, lowZ, isOrigin=True)
        newBlock.addEdgeArc(2, 3, x, y, lowZ, isOrigin=True)
        newBlock.addEdgeArc(4, 5, x, y, highZ, isOrigin=True)
        newBlock.addEdgeArc(6, 7, x, y, highZ, isOrigin=True)

        if (isAddBoundaryConditions):
            idxFace = len(self.faces)

            topFace = Face(f"{name}Top_{idxFace}")
            botFace = Face(f"{name}Bottom_{idxFace}")

            topFace.addSubFace(newBlock.topFace())
            botFace.addSubFace(newBlock.bottomFace())

            innerWallFace = Face(f"{name}InnerWall_{idxFace}", boundaryType='wall')
            innerWallFace.addSubFace(newBlock.backFace())

            outerWallFace = Face(f"{name}OuterWall_{idxFace}", boundaryType='wall')
            outerWallFace.addSubFace(newBlock.frontFace())

            self.addBoundary(topFace)
            self.addBoundary(botFace)
            self.addBoundary(innerWallFace)
            self.addBoundary(outerWallFace)

        return(newBlock)


    def createCylinderAlongZ(
            self,
            name: str,
            radius: float,
            lowZ: float,
            highZ: float,
            nx: int, ny: int, nz: int,
            x: float=0, y: float=0,
            isAddBoundaryConditions: bool=False
        ) -> Block:
        """
        Return
        ------
            (centerBlock, frontBlock, rightBlock, backBlock, leftBlock)
        """
        sqrt2 = np.sqrt(2)

        centerBlock = self.createBlock(name, [
            Point(x-radius/3, y-radius/3, lowZ),
            Point(x+radius/3, y-radius/3, lowZ),
            Point(x+radius/3, y+radius/3, lowZ),
            Point(x-radius/3, y+radius/3, lowZ),
            Point(x-radius/3, y-radius/3, highZ),
            Point(x+radius/3, y-radius/3, highZ),
            Point(x+radius/3, y+radius/3, highZ),
            Point(x-radius/3, y+radius/3, highZ)
        ], nx, ny, nz)

        frontBlock = self.addFront(centerBlock, name, [
            Point(x-radius/sqrt2, y-radius/sqrt2, lowZ),
            Point(x+radius/sqrt2, y-radius/sqrt2, lowZ),
            Point(x-radius/sqrt2, y-radius/sqrt2, highZ),
            Point(x+radius/sqrt2, y-radius/sqrt2, highZ),
        ], ny=ny)

        rightBlock = self.addRight(centerBlock, name, [
            frontBlock.points[1],
            Point(x+radius/sqrt2, y+radius/sqrt2, lowZ),
            frontBlock.points[5],
            Point(x+radius/sqrt2, y+radius/sqrt2, highZ),
        ], nx=frontBlock.ny)

        backBlock = self.addBack(centerBlock, name, [
            rightBlock.points[2],
            Point(x-radius/sqrt2, y+radius/sqrt2, lowZ),
            rightBlock.points[6],
            Point(x-radius/sqrt2, y+radius/sqrt2, highZ),
        ], ny=frontBlock.ny)

        leftBlock = self.addLeft(centerBlock, name, [
            frontBlock.points[0],
            backBlock.points[3],
            frontBlock.points[4],
            backBlock.points[7],
        ], nx=frontBlock.ny)

        frontBlock.addEdgeArc(0, 1, x=x, y=y, isOrigin=True)
        frontBlock.addEdgeArc(4, 5, x=x, y=y, isOrigin=True)
        rightBlock.addEdgeArc(1, 2, x=x, y=y, isOrigin=True)
        rightBlock.addEdgeArc(5, 6, x=x, y=y, isOrigin=True)
        backBlock.addEdgeArc(2, 3, x=x, y=y, isOrigin=True)
        backBlock.addEdgeArc(6, 7, x=x, y=y, isOrigin=True)
        leftBlock.addEdgeArc(3, 0, x=x, y=y, isOrigin=True)
        leftBlock.addEdgeArc(7, 4, x=x, y=y, isOrigin=True)


        if (isAddBoundaryConditions):
            idxFace = len(self.faces)

            topFace = Face(f"{name}Top_{idxFace}")
            botFace = Face(f"{name}Bottom_{idxFace}")

            for block in [centerBlock, frontBlock, rightBlock, backBlock, leftBlock]:
                topFace.addSubFace(block.topFace())
                botFace.addSubFace(block.bottomFace())

            wallFaceFront = Face(f"{name}WallFront_{idxFace}", boundaryType='wall')
            wallFaceFront.addSubFace(frontBlock.frontFace())
            wallFaceLeft = Face(f"{name}WallLeft_{idxFace}", boundaryType='wall')
            wallFaceLeft.addSubFace(leftBlock.leftFace())
            wallFaceRight = Face(f"{name}WallRight_{idxFace}", boundaryType='wall')
            wallFaceRight.addSubFace(rightBlock.rightFace())
            wallFaceBack = Face(f"{name}WallBack_{idxFace}", boundaryType='wall')
            wallFaceBack.addSubFace(backBlock.backFace())

            self.addBoundary(topFace)
            self.addBoundary(botFace)
            self.addBoundary(wallFaceFront)
            self.addBoundary(wallFaceLeft)
            self.addBoundary(wallFaceRight)
            self.addBoundary(wallFaceBack)

        return(centerBlock, frontBlock, rightBlock, backBlock, leftBlock)


    def createCubeWithHoleAlongZ(
            self,
            name: str,
            lowX: float, lowY: float, lowZ: float,
            highX: float, highY: float, highZ: float,
            radius: float,
            nx: int=1, ny: int=1, nz: int=1, nt: int=1,
            isHoleCylinder: bool=True,
            squareEdgeToHoleCenter=None,
            isAddBoundaryConditions: bool=False
        ):
        """
        Parameters
        ----------
        squareEdgeToHoleCenter : float
            Distance from the corner cube edge to the hole center along one axis.
            (default `radius/sqrt(2)`)

        Return
        ------
            frontBlock, backRightBlock, frontRightBlock, backBlock,
            backLeftBlock, frontLeftBlock, leftBlock, rightBlock
        """
        sqrt2 = np.sqrt(2)
        centerX = (highX+lowX)/2
        centerY = (highY+lowY)/2

        sqrEdge = squareEdgeToHoleCenter
        if (sqrEdge is None):
            sqrEdge = radius/sqrt2

        frontLeftBlock = self.createBlock(name, [
            Point(lowX,            lowY,            lowZ),
            Point(centerX-sqrEdge, lowY,            lowZ),
            Point(centerX-sqrEdge, centerY-sqrEdge, lowZ),
            Point(lowX,            centerY-sqrEdge, lowZ),
            Point(lowX,            lowY,            highZ),
            Point(centerX-sqrEdge, lowY,            highZ),
            Point(centerX-sqrEdge, centerY-sqrEdge, highZ),
            Point(lowX,            centerY-sqrEdge, highZ),
        ], nx, ny, nz)

        frontBlock = self.extrudeRight([frontLeftBlock], name, dx=2*sqrEdge, nx=nt)

        frontRightBlock = self.extrudeRight([frontBlock], name, dx=highX-centerX-sqrEdge, nx=nx)

        (leftBlock, rightBlock) = self.extrudeBack(
            [frontLeftBlock, frontRightBlock], name, dy=2*sqrEdge, ny=nt
        )

        (backLeftBlock, backRightBlock) = self.extrudeBack(
            [leftBlock, rightBlock], name, dy=highY-centerY-sqrEdge, ny=ny
        )

        backBlock = self.addRight(backLeftBlock, name, [
            backRightBlock.points[0],
            backRightBlock.points[3],
            backRightBlock.points[4],
            backRightBlock.points[7],
        ], nx=nt)

        frontLeftBlock.points[2].x = centerX-radius/sqrt2
        frontLeftBlock.points[2].y = centerY-radius/sqrt2
        frontLeftBlock.points[6].x = centerX-radius/sqrt2
        frontLeftBlock.points[6].y = centerY-radius/sqrt2
        frontRightBlock.points[3].x = centerX+radius/sqrt2
        frontRightBlock.points[3].y = centerY-radius/sqrt2
        frontRightBlock.points[7].x = centerX+radius/sqrt2
        frontRightBlock.points[7].y = centerY-radius/sqrt2
        backRightBlock.points[0].x = centerX+radius/sqrt2
        backRightBlock.points[0].y = centerY+radius/sqrt2
        backRightBlock.points[4].x = centerX+radius/sqrt2
        backRightBlock.points[4].y = centerY+radius/sqrt2
        backLeftBlock.points[1].x = centerX-radius/sqrt2
        backLeftBlock.points[1].y = centerY+radius/sqrt2
        backLeftBlock.points[5].x = centerX-radius/sqrt2
        backLeftBlock.points[5].y = centerY+radius/sqrt2

        if (isHoleCylinder):
            frontBlock.addEdgeArc(2, 3, x=centerX, y=centerY, isOrigin=True)
            frontBlock.addEdgeArc(6, 7, x=centerX, y=centerY, isOrigin=True)
            leftBlock.addEdgeArc(1, 2, x=centerX, y=centerY, isOrigin=True)
            leftBlock.addEdgeArc(5, 6, x=centerX, y=centerY, isOrigin=True)
            rightBlock.addEdgeArc(0, 3, x=centerX, y=centerY, isOrigin=True)
            rightBlock.addEdgeArc(4, 7, x=centerX, y=centerY, isOrigin=True)
            backBlock.addEdgeArc(0, 1, x=centerX, y=centerY, isOrigin=True)
            backBlock.addEdgeArc(4, 5, x=centerX, y=centerY, isOrigin=True)

        if (isAddBoundaryConditions):
            idxFace = len(self.faces)

            topFace = Face(f"{name}Top_{idxFace}")
            botFace = Face(f"{name}Bottom_{idxFace}")

            for block in [
                frontBlock, backRightBlock, frontRightBlock, backBlock,
                backLeftBlock, frontLeftBlock, leftBlock, rightBlock
            ]:
                topFace.addSubFace(block.topFace())
                botFace.addSubFace(block.bottomFace())

            wallFaceFront = Face(f"{name}WallFront_{idxFace}", boundaryType='wall')
            wallFaceFront.addSubFace(frontBlock.frontFace())
            wallFaceFront.addSubFace(frontLeftBlock.frontFace())
            wallFaceFront.addSubFace(frontRightBlock.frontFace())
            wallFaceBack = Face(f"{name}WallBack_{idxFace}", boundaryType='wall')
            wallFaceBack.addSubFace(backBlock.backFace())
            wallFaceBack.addSubFace(backRightBlock.backFace())
            wallFaceBack.addSubFace(backLeftBlock.backFace())
            wallFaceLeft = Face(f"{name}WallLeft_{idxFace}", boundaryType='wall')
            wallFaceLeft.addSubFace(frontLeftBlock.leftFace())
            wallFaceLeft.addSubFace(leftBlock.leftFace())
            wallFaceLeft.addSubFace(backLeftBlock.leftFace())
            wallFaceRight = Face(f"{name}WallRight_{idxFace}", boundaryType='wall')
            wallFaceRight.addSubFace(frontRightBlock.rightFace())
            wallFaceRight.addSubFace(rightBlock.rightFace())
            wallFaceRight.addSubFace(backRightBlock.rightFace())
            wallFaceHoleFront = Face(f"{name}WallHoleFront_{idxFace}", boundaryType='wall')
            wallFaceHoleFront.addSubFace(frontBlock.backFace())
            wallFaceHoleBack = Face(f"{name}WallHoleBack_{idxFace}", boundaryType='wall')
            wallFaceHoleBack.addSubFace(backBlock.frontFace())
            wallFaceHoleLeft = Face(f"{name}WallHoleLeft_{idxFace}", boundaryType='wall')
            wallFaceHoleLeft.addSubFace(rightBlock.leftFace())
            wallFaceHoleRight = Face(f"{name}WallHoleRight_{idxFace}", boundaryType='wall')
            wallFaceHoleRight.addSubFace(leftBlock.rightFace())

            self.addBoundary(topFace)
            self.addBoundary(botFace)
            self.addBoundary(wallFaceFront)
            self.addBoundary(wallFaceBack)
            self.addBoundary(wallFaceLeft)
            self.addBoundary(wallFaceRight)
            self.addBoundary(wallFaceHoleFront)
            self.addBoundary(wallFaceHoleBack)
            self.addBoundary(wallFaceHoleLeft)
            self.addBoundary(wallFaceHoleRight)

        return(
            frontBlock, backRightBlock, frontRightBlock, backBlock,
            backLeftBlock, frontLeftBlock, leftBlock, rightBlock
        )


    def createHexagonPrismAlongZ(
            self,
            name: str,
            zmin: float, zmax: float,
            pitch: float,
            x: float=0, y: float=0,
            nr: int=1,
            nt: int=1,
            nz: int=1,
            isAddBoundaryConditions: bool=False
        ) -> Block:
        """
        Block mesh instruction to create an hexagonal prism along the Z-axis.
        It is composed of 6 triagular prisms forming the 6 edges of the hexagon.

        Parameters
        ----------
        name : str
            Name of the hexagon prism
        zmin : float
            Bottom face Z position
        zmax : float
            Top face Z position
        picth : float
            Flat to flat distance
        x : float
            X-position of the prism center (default 0)
        y : float
            Y-position of the prism center (default 0)
        nr : int
            Number of cell layer in the radial direction (default 1)
        nt : int
            Number of cell per triangle in the azimutal direction (default 1)
        nz : int
            Number of cell along the Z direction (default 1)
        isAddBoundaryConditions : bool
            Flag to add automatically boundary faces on the prism (default False).
            If use lattice placement, it is recommended to set to 'true' if the
            internal faces need to be merged.

        Return
        ------
            (frontBlock, frontLeftBlock, frontRightBlock, backBlock, backRightBlock, backLeftBlock)
        """

        side = pitch/np.sqrt(3)

        frontBlock = self.createBlock(name, [
            Point(x-side/2, y-pitch/2, zmin),
            Point(x+side/2, y-pitch/2, zmin),
            Point(x, y, zmin),
            Point(x, y, zmin),
            Point(x-side/2, y-pitch/2, zmax),
            Point(x+side/2, y-pitch/2, zmax),
            Point(x, y, zmax),
            Point(x, y, zmax)
        ], nt, nr, nz)

        frontRightBlock = self.addRight(frontBlock, name, [
            Point(x+side, y, zmin),
            Point(x, y, zmin),
            Point(x+side, y, zmax),
            Point(x, y, zmax),
        ], nx=nt)

        frontLeftBlock = self.addLeft(frontBlock, name, [
            Point(x-side, y, zmin),
            Point(x, y, zmin),
            Point(x-side, y, zmax),
            Point(x, y, zmax),
        ], nx=nt)

        backBlock = self.addBack(frontBlock, name, [
            Point(x+side/2, y+pitch/2, zmin),
            Point(x-side/2, y+pitch/2, zmin),
            Point(x+side/2, y+pitch/2, zmax),
            Point(x-side/2, y+pitch/2, zmax),
        ], ny=nr)

        backLeftBlock = self.addLeft(backBlock, name, [
            frontLeftBlock.points[3],
            frontLeftBlock.points[0],
            frontLeftBlock.points[7],
            frontLeftBlock.points[4],
        ], nx=nt)

        backRightBlock = self.addRight(backBlock, name, [
            frontRightBlock.points[2],
            frontRightBlock.points[1],
            frontRightBlock.points[6],
            frontRightBlock.points[5],
        ], nx=nt)

        if (isAddBoundaryConditions):
            idxFace = len(self.faces)

            topFace = Face(f"{name}Top_{idxFace}")
            botFace = Face(f"{name}Bottom_{idxFace}")

            for block in [frontBlock, backRightBlock, frontRightBlock, backBlock, backLeftBlock, frontLeftBlock]:
                topFace.addSubFace(block.topFace())
                botFace.addSubFace(block.bottomFace())

            wallFaceFront = Face(f"{name}WallFront_{idxFace}", boundaryType='wall')
            wallFaceFront.addSubFace(frontBlock.frontFace())
            wallFaceFrontLeft = Face(f"{name}WallFrontLeft_{idxFace}", boundaryType='wall')
            wallFaceFrontLeft.addSubFace(frontLeftBlock.frontFace())
            wallFaceFrontRight = Face(f"{name}WallFrontRight_{idxFace}", boundaryType='wall')
            wallFaceFrontRight.addSubFace(frontRightBlock.frontFace())
            wallFaceBack = Face(f"{name}WallBack_{idxFace}", boundaryType='wall')
            wallFaceBack.addSubFace(backBlock.backFace())
            wallFaceBackLeft = Face(f"{name}WallBackLeft_{idxFace}", boundaryType='wall')
            wallFaceBackLeft.addSubFace(backLeftBlock.backFace())
            wallFaceBackRight = Face(f"{name}WallBackRight_{idxFace}", boundaryType='wall')
            wallFaceBackRight.addSubFace(backRightBlock.backFace())

            self.addBoundary(topFace)
            self.addBoundary(botFace)
            self.addBoundary(wallFaceFront)
            self.addBoundary(wallFaceFrontLeft)
            self.addBoundary(wallFaceFrontRight)
            self.addBoundary(wallFaceBack)
            self.addBoundary(wallFaceBackLeft)
            self.addBoundary(wallFaceBackRight)

        return(frontBlock, frontLeftBlock, frontRightBlock, backBlock, backRightBlock, backLeftBlock)


    def createHexagonPrismWithHoleAlongZ(
            self,
            name: str,
            zmin: float, zmax: float,
            pitch: float,
            radius: float,
            x: float=0, y: float=0,
            nr: int=1,
            nt: int=1,
            nz: int=1,
            isHoleCylinder: bool=True,
            isAddBoundaryConditions: bool=False
        ):
        """
        Block mesh instruction to create an hexagonal prism along the Z-axis.
        It is composed of 6 triagular prisms forming the 6 edges of the hexagon.

        Parameters
        ----------
        name : str
            Name of the hexagon prism
        zmin : float
            Bottom face Z position
        zmax : float
            Top face Z position
        picth : float
            Flat to flat distance
        radius : float
            Radius of the hole
        x : float
            X-position of the prism center (default 0)
        y : float
            Y-position of the prism center (default 0)
        nr : int
            Number of cell layer in the radial direction (default 1)
        nt : int
            Number of cell per triangle in the azimutal direction (default 1)
        nz : int
            Number of cell along the Z direction (default 1)
        isHoleCylinder : bool
            Create a cylindrical hole if True. If False, create an hexagonal
            hole. (default `True`)
        isAddBoundaryConditions : bool
            Flag to add automatically boundary faces on the prism (default False).
            If use lattice placement, it is recommended to set to 'true' if the
            internal faces need to be merged.

        Return
        ------
            (frontBlock, frontLeftBlock, frontRightBlock, backBlock, backRightBlock, backLeftBlock)
        """

        side = pitch/np.sqrt(3)
        halfSqrt3 = np.sqrt(3)/2

        frontBlock = self.createBlock(name, [
            Point(x-side/2, y-pitch/2, zmin),
            Point(x+side/2, y-pitch/2, zmin),
            Point(x+radius/2, y-radius*halfSqrt3, zmin),
            Point(x-radius/2, y-radius*halfSqrt3, zmin),
            Point(x-side/2, y-pitch/2, zmax),
            Point(x+side/2, y-pitch/2, zmax),
            Point(x+radius/2, y-radius*halfSqrt3, zmax),
            Point(x-radius/2, y-radius*halfSqrt3, zmax)
        ], nt, nr, nz)

        frontRightBlock = self.addRight(frontBlock, name, [
            Point(x+side, y, zmin),
            Point(x+radius, y, zmin),
            Point(x+side, y, zmax),
            Point(x+radius, y, zmax),
        ], nx=nt)

        frontLeftBlock = self.addLeft(frontBlock, name, [
            Point(x-side, y, zmin),
            Point(x-radius, y, zmin),
            Point(x-side, y, zmax),
            Point(x-radius, y, zmax),
        ], nx=nt)

        backBlock = self.createBlock(name, [
            Point(x-radius/2, y+radius*halfSqrt3, zmin),
            Point(x+radius/2, y+radius*halfSqrt3, zmin),
            Point(x+side/2, y+pitch/2, zmin),
            Point(x-side/2, y+pitch/2, zmin),
            Point(x-radius/2, y+radius*halfSqrt3, zmax),
            Point(x+radius/2, y+radius*halfSqrt3, zmax),
            Point(x+side/2, y+pitch/2, zmax),
            Point(x-side/2, y+pitch/2, zmax),
        ], nt, nr, nz)

        backLeftBlock = self.addLeft(backBlock, name, [
            frontLeftBlock.points[3],
            frontLeftBlock.points[0],
            frontLeftBlock.points[7],
            frontLeftBlock.points[4],
        ], nx=nt)

        backRightBlock = self.addRight(backBlock, name, [
            frontRightBlock.points[2],
            frontRightBlock.points[1],
            frontRightBlock.points[6],
            frontRightBlock.points[5],
        ], nx=nt)

        if (isHoleCylinder):
            frontBlock.addEdgeArc(2, 3, x=x, y=y, isOrigin=True)
            frontBlock.addEdgeArc(6, 7, x=x, y=y, isOrigin=True)
            backBlock.addEdgeArc(0, 1, x=x, y=y, isOrigin=True)
            backBlock.addEdgeArc(4, 5, x=x, y=y, isOrigin=True)
            frontRightBlock.addEdgeArc(2, 3, x=x, y=y, isOrigin=True)
            frontRightBlock.addEdgeArc(6, 7, x=x, y=y, isOrigin=True)
            frontLeftBlock.addEdgeArc(2, 3, x=x, y=y, isOrigin=True)
            frontLeftBlock.addEdgeArc(6, 7, x=x, y=y, isOrigin=True)
            backLeftBlock.addEdgeArc(0, 1, x=x, y=y, isOrigin=True)
            backLeftBlock.addEdgeArc(4, 5, x=x, y=y, isOrigin=True)
            backRightBlock.addEdgeArc(0, 1, x=x, y=y, isOrigin=True)
            backRightBlock.addEdgeArc(4, 5, x=x, y=y, isOrigin=True)

        if (isAddBoundaryConditions):
            idxFace = len(self.faces)

            topFace = Face(f"{name}Top_{idxFace}")
            botFace = Face(f"{name}Bottom_{idxFace}")

            for block in [frontBlock, backRightBlock, frontRightBlock, backBlock, backLeftBlock, frontLeftBlock]:
                topFace.addSubFace(block.topFace())
                botFace.addSubFace(block.bottomFace())

            wallFaceFront = Face(f"{name}WallFront_{idxFace}", boundaryType='wall')
            wallFaceFront.addSubFace(frontBlock.frontFace())
            wallFaceFrontLeft = Face(f"{name}WallFrontLeft_{idxFace}", boundaryType='wall')
            wallFaceFrontLeft.addSubFace(frontLeftBlock.frontFace())
            wallFaceFrontRight = Face(f"{name}WallFrontRight_{idxFace}", boundaryType='wall')
            wallFaceFrontRight.addSubFace(frontRightBlock.frontFace())
            wallFaceBack = Face(f"{name}WallBack_{idxFace}", boundaryType='wall')
            wallFaceBack.addSubFace(backBlock.backFace())
            wallFaceBackLeft = Face(f"{name}WallBackLeft_{idxFace}", boundaryType='wall')
            wallFaceBackLeft.addSubFace(backLeftBlock.backFace())
            wallFaceBackRight = Face(f"{name}WallBackRight_{idxFace}", boundaryType='wall')
            wallFaceBackRight.addSubFace(backRightBlock.backFace())
            wallFaceHole = Face(f"{name}WallHole_{idxFace}", boundaryType='wall')
            wallFaceHole.addSubFace(frontBlock.backFace())
            wallFaceHole.addSubFace(frontRightBlock.backFace())
            wallFaceHole.addSubFace(frontLeftBlock.backFace())
            wallFaceHole.addSubFace(backBlock.frontFace())
            wallFaceHole.addSubFace(backRightBlock.frontFace())
            wallFaceHole.addSubFace(backLeftBlock.frontFace())

            self.addBoundary(topFace)
            self.addBoundary(botFace)
            self.addBoundary(wallFaceFront)
            self.addBoundary(wallFaceFrontLeft)
            self.addBoundary(wallFaceFrontRight)
            self.addBoundary(wallFaceBack)
            self.addBoundary(wallFaceBackLeft)
            self.addBoundary(wallFaceBackRight)
            self.addBoundary(wallFaceHole)

        return(frontBlock, frontLeftBlock, frontRightBlock, backBlock, backRightBlock, backLeftBlock)


    def createEdgeHexagonPrismAlongZ(
            self,
            name: str,
            zmin: float, zmax: float,
            pitch: float,
            distanceCenterToEdge: float,
            edgeFaceOrientation: float=0,
            x: float=0, y: float=0,
            nr: int=1,
            nt: int=1,
            nz: int=1,
            isAddBoundaryConditions: bool=False
        ):
        """
        Block mesh instruction to create an hexagonal prism along the Z-axis.
        It is cutted on one side to fit at the edge of a fuel assembly with a
        wrapper.

        Parameters
        ----------
        name : str
            Name of the hexagon prism
        zmin : float
            Bottom face Z position
        zmax : float
            Top face Z position
        picth : float
            Flat to flat distance
        distanceCenterToEdge : float
            Distance from the hexagon center to the edge (e.g wrapper wall)
        edgeFaceOrientation : float
            Edge face orientation in rad
        x : float
            X-position of the prism center (default 0)
        y : float
            Y-position of the prism center (default 0)
        nr : int
            Number of cell layer in the radial direction (default 1)
        nt : int
            Number of cell per triangle in the azimutal direction (default 1)
        nz : int
            Number of cell along the Z direction (default 1)
        isAddBoundaryConditions : bool
            Flag to add automatically boundary faces on the prism (default False).
            If use lattice placement, it is recommended to set to 'true' if the
            internal faces need to be merged.

        Return
        ------
            (frontBlock, backBlock)
        """

        side = pitch/np.sqrt(3)
        halfSqrt3 = np.sqrt(3)/2

        cost = np.cos(edgeFaceOrientation)
        sint = np.sin(edgeFaceOrientation)

        frontBlock = self.createBlock(name, [
            Point(x-side/2, y-pitch/2, zmin),
            Point(x+distanceCenterToEdge, y-pitch/2, zmin),
            Point(x+distanceCenterToEdge, y, zmin),
            Point(x-side, y, zmin),
            Point(x-side/2, y-pitch/2, zmax),
            Point(x+distanceCenterToEdge, y-pitch/2, zmax),
            Point(x+distanceCenterToEdge, y, zmax),
            Point(x-side, y, zmax),
        ], nr, nt, nz)

        frontBlock.translate(dx=-x, dy=-y)
        frontBlock.rotateZ(theta=edgeFaceOrientation)
        frontBlock.translate(dx=x, dy=y)

        x2 = x + distanceCenterToEdge*cost - pitch/2*sint
        y2 = y + pitch/2*cost + distanceCenterToEdge*sint
        x3 = x - side/2*cost - pitch/2*sint
        y3 = y + pitch/2*cost - side/2*sint
        backBlock = self.addBack(frontBlock, name, [
            Point(x2, y2, zmin),
            Point(x3, y3, zmin),
            Point(x2, y2, zmax),
            Point(x3, y3, zmax),
        ], ny=nt)


        if (isAddBoundaryConditions):
            idxFace = len(self.faces)

            topFace = Face(f"{name}Top_{idxFace}")
            botFace = Face(f"{name}Bottom_{idxFace}")

            for block in [frontBlock, backBlock]:
                topFace.addSubFace(block.topFace())
                botFace.addSubFace(block.bottomFace())

            wallFaceFront = Face(f"{name}WallFront_{idxFace}", boundaryType='wall')
            wallFaceFront.addSubFace(frontBlock.frontFace())
            wallFaceFrontLeft = Face(f"{name}WallFrontLeft_{idxFace}", boundaryType='wall')
            wallFaceFrontLeft.addSubFace(frontBlock.leftFace())
            wallFaceBack = Face(f"{name}WallBack_{idxFace}", boundaryType='wall')
            wallFaceBack.addSubFace(backBlock.backFace())
            wallFaceBackLeft = Face(f"{name}WallBackLeft_{idxFace}", boundaryType='wall')
            wallFaceBackLeft.addSubFace(backBlock.leftFace())
            wallFaceEdge = Face(f"{name}WallEdge_{idxFace}", boundaryType='wall')
            wallFaceEdge.addSubFace(backBlock.rightFace())
            wallFaceEdge.addSubFace(frontBlock.rightFace())

            self.addBoundary(topFace)
            self.addBoundary(botFace)
            self.addBoundary(wallFaceFront)
            self.addBoundary(wallFaceFrontLeft)
            self.addBoundary(wallFaceBack)
            self.addBoundary(wallFaceBackLeft)
            self.addBoundary(wallFaceEdge)

        return(frontBlock, backBlock)


    def createCornerHexagonPrismAlongZ(
            self,
            name: str,
            zmin: float, zmax: float,
            pitch: float,
            distanceCenterToEdge: float,
            edgeFaceOrientation: float=0,
            x: float=0, y: float=0,
            nr: int=1,
            nt: int=1,
            nz: int=1,
            isAddBoundaryConditions: bool=False
        ):
        """
        Block mesh instruction to create an hexagonal prism along the Z-axis.
        It is cutted on two sides to fit at the corner of a fuel assembly with a
        wrapper.

        Parameters
        ----------
        name : str
            Name of the hexagon prism
        zmin : float
            Bottom face Z position
        zmax : float
            Top face Z position
        picth : float
            Flat to flat distance
        distanceCenterToEdge : float
            Distance from the hexagon center to the edge (e.g wrapper wall)
        edgeFaceOrientation : float
            Edge face orientation in rad
        x : float
            X-position of the prism center (default 0)
        y : float
            Y-position of the prism center (default 0)
        nr : int
            Number of cell layer in the radial direction (default 1)
        nt : int
            Number of cell per triangle in the azimutal direction (default 1)
        nz : int
            Number of cell along the Z direction (default 1)
        isAddBoundaryConditions : bool
            Flag to add automatically boundary faces on the prism (default False).
            If use lattice placement, it is recommended to set to 'true' if the
            internal faces need to be merged.

        Return
        ------
            (frontBlock, backBlock, wedgeBlock)
        """

        side = pitch/np.sqrt(3)
        halfSqrt3 = np.sqrt(3)/2

        cost = np.cos(edgeFaceOrientation)
        sint = np.sin(edgeFaceOrientation)

        xc = x -side/2
        yc = y -side * (1 - 1/np.sqrt(2))

        frontBlock = self.createBlock(name, [
            Point(x-side/2, y-pitch/2, zmin),
            Point(x+distanceCenterToEdge, y-pitch/2, zmin),
            Point(x+distanceCenterToEdge, y + 2*distanceCenterToEdge * (1 - 1/np.sqrt(2)), zmin),
            Point(xc, yc, zmin),
            Point(x-side/2, y-pitch/2, zmax),
            Point(x+distanceCenterToEdge, y-pitch/2, zmax),
            Point(x+distanceCenterToEdge, y + 2*distanceCenterToEdge * (1 - 1/np.sqrt(2)), zmax),
            Point(xc, yc, zmax),
        ], nr, nt, nz)

        frontBlock.translate(dx=-x, dy=-y)
        frontBlock.rotateZ(theta=edgeFaceOrientation)
        frontBlock.translate(dx=x, dy=y)

        wedgeBlock = self.addLeft(frontBlock, name, [
            Point(x - side*cost, y -side*sint, zmin),
            Point(frontBlock.points[3].x, frontBlock.points[3].y, zmin),
            Point(x - side*cost, y -side*sint, zmax),
            Point(frontBlock.points[3].x, frontBlock.points[3].y, zmax),
        ], nx=nt)

        dx = - side + (side/2 + distanceCenterToEdge)/2
        dy = (side/2 + distanceCenterToEdge) * halfSqrt3
        x2 = x + dx*cost - dy*sint
        y2 = y + dx*sint + dy*cost
        backBlock = self.addLeft(wedgeBlock, name, [
            Point(x2, y2, zmin),
            frontBlock.points[2],
            Point(x2, y2, zmax),
            frontBlock.points[6],
        ], nx=nr)


        if (isAddBoundaryConditions):
            idxFace = len(self.faces)

            topFace = Face(f"{name}Top_{idxFace}")
            botFace = Face(f"{name}Bottom_{idxFace}")

            for block in [frontBlock, backBlock, wedgeBlock]:
                topFace.addSubFace(block.topFace())
                botFace.addSubFace(block.bottomFace())

            wallFaceFront = Face(f"{name}WallFront_{idxFace}", boundaryType='wall')
            wallFaceFront.addSubFace(frontBlock.frontFace())
            wallFaceRight = Face(f"{name}WallRight_{idxFace}", boundaryType='wall')
            wallFaceRight.addSubFace(frontBlock.rightFace())
            wallFaceBack = Face(f"{name}WallBack_{idxFace}", boundaryType='wall')
            wallFaceBack.addSubFace(backBlock.leftFace())
            wallFaceLeft = Face(f"{name}WallLeft_{idxFace}", boundaryType='wall')
            wallFaceLeft.addSubFace(backBlock.frontFace())
            wallFaceEdge = Face(f"{name}WallEdge_{idxFace}", boundaryType='wall')
            wallFaceEdge.addSubFace(wedgeBlock.frontFace())

            internalFace1 = Face(f"{name}Internal1_{idxFace}", boundaryType='int')
            internalFace1.addSubFace(frontBlock.backFace())
            internalFace2 = Face(f"{name}Internal2_{idxFace}", boundaryType='int')
            internalFace2.addSubFace(backBlock.backFace())

            self.mergePatchPairs.append((internalFace1, internalFace2))

            self.addBoundary(topFace)
            self.addBoundary(botFace)
            self.addBoundary(wallFaceFront)
            self.addBoundary(wallFaceRight)
            self.addBoundary(wallFaceBack)
            self.addBoundary(wallFaceLeft)
            self.addBoundary(wallFaceEdge)
            self.addBoundary(internalFace1)
            self.addBoundary(internalFace2)

        return(frontBlock, backBlock, wedgeBlock)


    def createTriangularChannel(
            self,
            name: str,
            lowZ: float, highZ: float,
            pitch: float,
            radius: float,
            thetaZ: float=0,
            x: float=0, y: float=0,
            nx: int=1,
            ny: int=1,
            nz: int=1,
            isAddBoundaryConditions: bool=False
        ):
        """
        Parameters
        ----------
        thetaZ : float
            Rotate block along the Z-axis in rad (default 0)

        Return
        ------
            frontRightBlock, frontLeftBlock, leftBlock, rightBlock,
            backLeftBlock, backRightBlock
        """
        sqrt3 = np.sqrt(3)

        h = pitch/(sqrt3*2)
        d = sqrt3*pitch/2

        # Block 1
        frontRightBlock = self.createBlock(name, [
            Point(x, y - h, lowZ),
            Point(x + pitch/2 - radius, y - h, lowZ),
            Point(x + pitch/2 - radius*sqrt3/2, y - h + radius/2, lowZ),
            Point(x, y, lowZ),
            Point(x, y - h, highZ),
            Point(x + pitch/2 - radius, y - h, highZ),
            Point(x + pitch/2 - radius*sqrt3/2, y - h + radius/2, highZ),
            Point(x, y, highZ),
        ], nx=nx, ny=ny, nz=nz)

        # Block 2
        frontLeftBlock = self.addLeft(frontRightBlock, name, [
            Point(x - pitch/2 + radius, y - h, lowZ),
            Point(x - pitch/2 + radius*sqrt3/2, y - h + radius/2, lowZ),
            Point(x - pitch/2 + radius, y - h, highZ),
            Point(x - pitch/2 + radius*sqrt3/2, y - h + radius/2, highZ),
        ], nx=nx)

        # Block 3
        leftBlock = self.addBack(frontLeftBlock, name, [
            Point(x - pitch/2 + pitch/2/2, y - h + pitch/2 * sqrt3/2, lowZ),
            Point(x - pitch/2 + radius/2, y - h + radius*sqrt3/2, lowZ),
            Point(x - pitch/2 + pitch/2/2, y - h + pitch/2 * sqrt3/2, highZ),
            Point(x - pitch/2 + radius/2, y - h + radius*sqrt3/2, highZ),
        ], ny=ny)

        # Block 4
        rightBlock = self.addBack(frontRightBlock, name, [
            Point(x + pitch/2 - radius/2, y - h + radius*sqrt3/2, lowZ),
            Point(x + pitch/2 - pitch/2/2, y - h + pitch/2 * sqrt3/2, lowZ),
            Point(x + pitch/2 - radius/2, y - h + radius*sqrt3/2, highZ),
            Point(x + pitch/2 - pitch/2/2, y - h + pitch/2 * sqrt3/2, highZ),
        ], ny=ny)

        # Block 5
        backLeftBlock = self.addRight(leftBlock, name, [
            Point(x, y + d-h-radius, lowZ),
            Point(x - radius/2, y + d-h-radius*sqrt3/2, lowZ),
            Point(x, y + d-h-radius, highZ),
            Point(x - radius/2, y + d-h-radius*sqrt3/2, highZ),
        ], nx=nx)

        # Block 6
        backRightBlock = self.addLeft(rightBlock, name, [
            Point(x , y + d-h-radius, lowZ),
            Point(x + radius/2, y + d-h-radius*sqrt3/2, lowZ),
            Point(x , y + d-h-radius, highZ),
            Point(x + radius/2, y + d-h-radius*sqrt3/2, highZ),
        ], nx=nx)

        frontRightBlock.addEdgeArc(1, 2, x=x + pitch/2, y=y - h, isOrigin=True)
        frontRightBlock.addEdgeArc(5, 6, x=x + pitch/2, y=y - h, isOrigin=True)
        frontLeftBlock.addEdgeArc(0, 3, x=x - pitch/2, y=y - h, isOrigin=True)
        frontLeftBlock.addEdgeArc(4, 7, x=x - pitch/2, y=y - h, isOrigin=True)
        leftBlock.addEdgeArc(0, 3, x=x - pitch/2, y=y - h, isOrigin=True)
        leftBlock.addEdgeArc(4, 7, x=x - pitch/2, y=y - h, isOrigin=True)
        rightBlock.addEdgeArc(1, 2, x=x + pitch/2, y=y - h, isOrigin=True)
        rightBlock.addEdgeArc(5, 6, x=x + pitch/2, y=y - h, isOrigin=True)
        backLeftBlock.addEdgeArc(1, 2, x=x, y=y+d-h, isOrigin=True)
        backLeftBlock.addEdgeArc(5, 6, x=x, y=y+d-h, isOrigin=True)
        backRightBlock.addEdgeArc(0, 3, x=x, y=y+d-h, isOrigin=True)
        backRightBlock.addEdgeArc(4, 7, x=x, y=y+d-h, isOrigin=True)

        # Rotate all points
        points = self.getUniquePoints([
            frontRightBlock, frontLeftBlock, leftBlock, rightBlock,
            backLeftBlock, backRightBlock
        ], withEdge=True)
        for point in points:
            point.translate(dx=-x, dy=-y)
            point.rotateZ(theta=thetaZ)
            point.translate(dx=x, dy=y)


        if (isAddBoundaryConditions):
            idxFace = len(self.faces)

            topFace = Face(f"{name}Top_{idxFace}")
            botFace = Face(f"{name}Bottom_{idxFace}")

            for block in [
                frontRightBlock, frontLeftBlock, leftBlock, rightBlock,
                backLeftBlock, backRightBlock
            ]:
                topFace.addSubFace(block.topFace())
                botFace.addSubFace(block.bottomFace())


            wallFaceFront = Face(f"{name}WallFluidFront_{idxFace}", boundaryType='wall')
            wallFaceFront.addSubFace(frontLeftBlock.frontFace())
            wallFaceFront.addSubFace(frontRightBlock.frontFace())
            wallFaceRight = Face(f"{name}WallFluidRight_{idxFace}", boundaryType='wall')
            wallFaceRight.addSubFace(rightBlock.backFace())
            wallFaceRight.addSubFace(backRightBlock.backFace())
            wallFaceLeft = Face(f"{name}WallFluidLeft_{idxFace}", boundaryType='wall')
            wallFaceLeft.addSubFace(leftBlock.backFace())
            wallFaceLeft.addSubFace(backLeftBlock.backFace())
            wallFaceHoleLeft = Face(f"{name}WallHoleLeft_{idxFace}", boundaryType='wall')
            wallFaceHoleLeft.addSubFace(leftBlock.leftFace())
            wallFaceHoleLeft.addSubFace(frontLeftBlock.leftFace())
            wallFaceHoleRight = Face(f"{name}WallHoleRight_{idxFace}", boundaryType='wall')
            wallFaceHoleRight.addSubFace(rightBlock.rightFace())
            wallFaceHoleRight.addSubFace(frontRightBlock.rightFace())
            wallFaceHoleBack = Face(f"{name}WallHoleBack_{idxFace}", boundaryType='wall')
            wallFaceHoleBack.addSubFace(backLeftBlock.rightFace())
            wallFaceHoleBack.addSubFace(backRightBlock.leftFace())

            internalFace1 = Face(f"{name}Internal1_{idxFace}", boundaryType='int')
            internalFace1.addSubFace(backRightBlock.frontFace())
            internalFace2 = Face(f"{name}Internal2_{idxFace}", boundaryType='int')
            internalFace2.addSubFace(backLeftBlock.frontFace())

            # self.mergePatchPairs.append((internalFace1, internalFace2))

            self.addBoundary(topFace)
            self.addBoundary(botFace)
            self.addBoundary(wallFaceFront)
            self.addBoundary(wallFaceRight)
            self.addBoundary(wallFaceLeft)
            self.addBoundary(wallFaceHoleLeft)
            self.addBoundary(wallFaceHoleRight)
            self.addBoundary(wallFaceHoleBack)
            self.addBoundary(internalFace1)
            self.addBoundary(internalFace2)

        return(
            frontRightBlock, frontLeftBlock, leftBlock, rightBlock,
            backLeftBlock, backRightBlock
        )


    def addPipe1DFromDirection(
            self,
            name: str,
            originPosition: Vector | Block,
            direction: Vector,
            length: float,
            equivalentHydraulicDiameter: float,
            n: int=1,
            elbowRadius: float=0,
            isAddBoundaryConditions: bool=False,
            originPositionOutletFaceName: str='top'
        ) -> Block:
        """
        Create a 1D square pipe starting from an origin position and extending
        along a direction.

        Parameters
        ----------
        name : str
            Name of the pipe
        originPosition : Vector | Block
            Origin point as a vector or as a Block. If Block, the origin point
            is computed as the barycenter of the top face.
        direction : Vector
            Direction of the pipe, can be unnormalized
        length : float
            Length of the pipe
        equivalentHydraulicDiameter : float
            Equivalent hydraulic diameter of a cylindrical pipe. Size of the
            square computed to conserve cross-section flow area as:

                `side = sqrt(pi) * eqDh`

        n: int
            Number cells along the pipe (default `1`)
        elbowRadius : float
            Offset the pipe knowing the elbow curvature radius and if
            `originPosition` is a `Block` (default `0`).
        isAddBoundaryConditions : bool
            If `True`, set left, right, front and back faces to `empty` BC
            (default `False`). If `originPosition` is a `Block`, connect the
            target block to the new pipe via the `connectPipes` method.
        originPositionOutletFaceName : str
            Facename of the origin block to attach the BC (default `top`).

        Return
        ------
            (block)
        """
        check_type("originPosition", originPosition, (Vector, Block))
        check_type("direction", direction, Vector)
        check_type("length", length, (float, int))
        check_type("equivalentHydraulicDiameter", equivalentHydraulicDiameter, (float, int))
        check_type("n", n, int)
        check_positive("n", n)
        check_type("elbowRadius", elbowRadius, (float, int))
        check_positive("elbowRadius", elbowRadius)
        check_type("originPositionOutletFaceName", originPositionOutletFaceName, str)
        check_value("originPositionOutletFaceName", originPositionOutletFaceName, _FACE_NAME_TYPES)

        connectingBlock = originPosition

        offset = 0
        # Recompute origin position and offset knowing the elbow radius
        if (isinstance(connectingBlock, Block)):
            pipe1outletFace = connectingBlock.getFace(originPositionOutletFaceName)
            # pipe1inletFace = connectingBlock.getOppositeFace(originPositionOutletFaceName)

            originPosition = connectingBlock.getFaceBarycenter(pipe1outletFace)
            # pipeInlet = connectingBlock.getFaceBarycenter(pipe1inletFace)
            # connectingBlockDir = originPosition - pipeInlet
            connectingBlockDir = connectingBlock.getFaceNormal(originPositionOutletFaceName)

            theta = np.arccos(connectingBlockDir.dot(direction) / (connectingBlockDir.norm() * direction.norm()))
            offset = elbowRadius * np.tan(theta/2)

            connectingBlockDir.normalize(offset)

            originPosition += connectingBlockDir

        # If conserve hydraulic diameter between a square and circle
        # cross-section flow
        # halfDh = equivalentHydraulicDiameter/2
        # Conserve cross-section flow area between a square and circle
        halfDh = equivalentHydraulicDiameter/2

        x = originPosition.x
        y = originPosition.y
        z = originPosition.z
        direction.normalize()
        dx = (length + 2*offset) * direction.x
        dy = (length + 2*offset) * direction.y
        dz = (length + 2*offset) * direction.z

        theta = np.arccos(direction.z)
        phi = np.arctan2(direction.y, direction.x)

        # Create the block at the origin
        block = self.createBlock(name, [
            Point(-halfDh, -halfDh, -length/2),
            Point(+halfDh, -halfDh, -length/2),
            Point(+halfDh, +halfDh, -length/2),
            Point(-halfDh, +halfDh, -length/2),
            Point(-halfDh, -halfDh, +length/2),
            Point(+halfDh, -halfDh, +length/2),
            Point(+halfDh, +halfDh, +length/2),
            Point(-halfDh, +halfDh, +length/2),
        ], nx=1, ny=1, nz=n)

        # Rotate all points and place to final position
        points = self.getUniquePoints([block], withEdge=True)
        for point in points:
            point.rotateY(theta=-theta)
            point.rotateZ(theta=phi)
            point.translate(dx=x+dx/2, dy=y+dy/2, dz=z+dz/2)

        if (isAddBoundaryConditions):
            # 'block' is always from bottom to top
            self.pipeWallBC.addSubFace(block.leftFace())
            self.pipeWallBC.addSubFace(block.rightFace())
            self.pipeWallBC.addSubFace(block.frontFace())
            self.pipeWallBC.addSubFace(block.backFace())

            if (isinstance(connectingBlock, Block)):
                self.connectPipes(
                    pipe1=connectingBlock,
                    pipe2=block,
                    elbowRadius=elbowRadius,
                    pipe1outletFaceName=originPositionOutletFaceName
                )

        return(block)


    def addPipe1DFrom2Points(
            self,
            name: str,
            originPosition: Vector | Block,
            finalPosition: Vector | Block,
            equivalentHydraulicDiameter: float,
            n: int=1,
            elbowRadius: float=0,
            isAddBoundaryConditions: bool=False,
            originPositionOutletFaceName: str='top',
            finalPositionInletFaceName: str='bottom',
            tolerance: float=1e-12
        ) -> Block:
        """
        Create a 1D square pipe between two points. The length is automatically
        computed. Works analytically in U-shape curve and numerically in S-shape
        curve.

        Parameters
        ----------
        name : str
            Name of the pipe
        originPosition : Vector | Block
            Origin point as a vector or as a Block. If Block, the origin point
            is computed as the barycenter of the top face.
        finalPosition: Vector | Block
            Final point as a vector or as a Block. If Block, the final point
            is computed as the barycenter of the bottom face.
        equivalentHydraulicDiameter : float
            Equivalent hydraulic diameter of a cylindrical pipe. Size of the
            square computed to conserve cross-section flow area as
            s = sqrt(pi) * eqDh / 2
        n: int
            Number cells along the pipe (default `1`)
        elbowRadius : float
            Offset the pipe knowing the elbow curvature radius (default `0`).
        isAddBoundaryConditions : bool
            If `True`, set left, right, front and back face to `empty` BC.
            If `True` and `finalPosition` is a `Block`, connect the new pipe to
            `finalPosition` block. Can be used to close loop.
        originPositionOutletFaceName : str
            Facename of the origin block to attach the BC (default `top`).
        finalPositionInletFaceName : str
            Facename of the new pipe block to attach the BC (default `bottom`).
        tolerance : float
            Tolerance when solving for connecting two pipes in S-shape
            (default 1e-12).

        Return
        ------
            (block)
        """
        check_type("originPosition", originPosition, (Vector, Block))
        check_type("finalPosition", finalPosition, (Vector, Block))
        check_type("originPositionOutletFaceName", originPositionOutletFaceName, str)
        check_value("originPositionOutletFaceName", originPositionOutletFaceName, _FACE_NAME_TYPES)
        check_type("finalPositionInletFaceName", finalPositionInletFaceName, str)
        check_value("finalPositionInletFaceName", finalPositionInletFaceName, _FACE_NAME_TYPES)


        # Extract origin position as a vector
        if (isinstance(originPosition, Vector)):
            originPositionVec = originPosition
        elif (isinstance(originPosition, Block)):
            pipe1outletFace = originPosition.getFace(originPositionOutletFaceName)
            originPositionVec = originPosition.getFaceBarycenter(pipe1outletFace)

        # Extract final position as a vector
        if (isinstance(finalPosition, Vector)):
            finalPositionVec = finalPosition
        elif (isinstance(finalPosition, Block)):
            pipe2inletFace = finalPosition.getFace(finalPositionInletFaceName)
            finalPositionVec = finalPosition.getFaceBarycenter(pipe2inletFace)

        direction: Vector = finalPositionVec - originPositionVec

        # Recompute direction if both origin and final position come from Block
        # Suppose the same elbow radius
        if (isinstance(originPosition, Block) and isinstance(finalPosition, Block)):
            # pipe1inletFace = originPosition.getOppositeFace(originPositionOutletFaceName)
            # pipe2outletFace = finalPosition.getOppositeFace(finalPositionInletFaceName)
            # pipe1Inlet = originPosition.getFaceBarycenter(pipe1inletFace)
            # pipe2Outlet = finalPosition.getFaceBarycenter(pipe2outletFace)

            # originDir = originPositionVec - pipe1Inlet
            # finalDir = pipe2Outlet - finalPositionVec
            originDir = originPosition.getFaceNormal(originPositionOutletFaceName)
            finalDir = finalPosition.getFaceNormal(finalPositionInletFaceName)

            rotationAxisOrigin = direction.cross(originDir)
            rotationAxisFinal = direction.cross(finalDir)

            isUshape = rotationAxisOrigin.dot(rotationAxisFinal) > 0

            R1 = rotationAxisOrigin.cross(originDir)
            R1.normalize(elbowRadius)
            R2 = rotationAxisFinal.cross(finalDir)
            R2.normalize(elbowRadius)

            if (isUshape):
                direction = direction + R1 + R2
            else:
                R1p = copy(R1)
                for i in range(1000):
                    newdir = direction + R1 + R2 + R1p + R1p

                    costheta = R1p.dot(newdir) / (R1p.norm() * newdir.norm())

                    R1p.rotateY(0.1*costheta)

                    if (abs(costheta) < tolerance):
                        break

                direction = direction + R1 + R2 + R1p + R1p

        length = direction.norm()

        block = self.addPipe1DFromDirection(
            name=name,
            originPosition=originPosition,
            direction=direction,
            length=length,
            equivalentHydraulicDiameter=equivalentHydraulicDiameter,
            n=n,
            elbowRadius=elbowRadius,
            isAddBoundaryConditions=isAddBoundaryConditions,
            originPositionOutletFaceName=originPositionOutletFaceName
        )

        if (isAddBoundaryConditions and isinstance(finalPosition, Block)):
            self.connectPipes(
                pipe1=block,
                pipe2=finalPosition,
                elbowRadius=elbowRadius,
                pipe1outletFaceName='top', # 'block' is always oriented from bottom to top
                pipe2inletFaceName=finalPositionInletFaceName
            )

        return(block)
    
    def addPipe1DFrom2PointsSCurve(
        self,
        name: str,
        originPosition: Block,
        finalPosition: Block,
        equivalentHydraulicDiameter: float,
        n: int = 1,
        elbowRadius: float = 0,
        isAddBoundaryConditions: bool = False,
        originPositionOutletFaceName: str = 'top',
        finalPositionInletFaceName: str = 'bottom',
        intermediateDirection: Vector = Vector(0, 1, 0),  # Default vertical
        intermediateLength: float = 5.0  # Default length of intermediate pipe
        ) -> tuple[Block, Block]:
        """
        Create a 1D pipe between two blocks using an intermediate segment.
        Useful for S-shaped connections.

        Returns
        -------
        Tuple of (intermediatePipe, finalPipe)
        """
        check_type("originPosition", originPosition, Block)
        check_type("finalPosition", finalPosition, Block)
        check_type("intermediateDirection", intermediateDirection, Vector)

        # Get pipes discretization

        originFace = originPosition.getFace(originPositionOutletFaceName)
        originVec = originPosition.getFaceBarycenter(originFace)

        # Manual copy and normalize
        normalizedDirection = Vector(
            intermediateDirection.x,
            intermediateDirection.y,
            intermediateDirection.z
        )
        normalizedDirection.normalize()

        
        intermediateOutletVec = Vector(
            originVec.x + normalizedDirection.x * intermediateLength,
            originVec.y + normalizedDirection.y * intermediateLength,
            originVec.z + normalizedDirection.z * intermediateLength
        )
        

        finalFace = finalPosition.getFace(finalPositionInletFaceName)
        finalVec = finalPosition.getFaceBarycenter(finalFace)

        finalLength = (finalVec - intermediateOutletVec).norm()
        totalLength = intermediateLength + finalLength

        intermediateN = math.ceil(n * intermediateLength / totalLength)
        finalN = n - intermediateN


        # Create intermediate pipe
        intermediatePipe = self.addPipe1DFromDirection(
            name=f"{name}_intermediate",
            originPosition=originPosition,
            direction=intermediateDirection,
            length=intermediateLength,
            equivalentHydraulicDiameter=equivalentHydraulicDiameter,
            n=intermediateN,
            elbowRadius=elbowRadius,
            isAddBoundaryConditions=isAddBoundaryConditions,
            originPositionOutletFaceName=originPositionOutletFaceName
        )

        # Create final connecting pipe
        finalPipe = self.addPipe1DFrom2Points(
            name=f"{name}_final",
            originPosition=intermediatePipe,
            finalPosition=finalPosition,
            equivalentHydraulicDiameter=equivalentHydraulicDiameter,
            n=finalN,
            elbowRadius=elbowRadius,
            isAddBoundaryConditions=isAddBoundaryConditions,
            originPositionOutletFaceName='top',  # intermediate pipe always goes bottom → top
            finalPositionInletFaceName=finalPositionInletFaceName
        )

        return intermediatePipe, finalPipe



    def connectPipes(
            self,
            pipe1: Block,
            pipe2: Block,
            elbowRadius: float = 0,
            pipe1outletFaceName: str = 'top',
            pipe2inletFaceName: str = 'bottom',
        ) -> None:
        """
        Connect two pipes, i.e `pipe1` outlet to `pipe2` inlet. Create a cyclic
        boundary condition.

        Parameters
        ----------
        pipe1 : Block
            Block to connect
        pipe2 : Block
            Block to connect
        elbowRadius : float
            Offset the pipe knowing the elbow curvature radius (default `0`).
        pipe1outletFaceName : str
            Name of the face to connect the `pipe1` (default `'top'`).
            Options: (`top`, `bottom`, `left`, `right`, `front`, `back`).
        pipe2inletFaceName : str
            Name of the face to connect the `pipe2` (default `'bottom'`).
            Options: (`top`, `bottom`, `left`, `right`, `front`, `back`).
        """
        check_type("pipe1", pipe1, Block)
        check_type("pipe2", pipe2, Block)
        check_type("elbowRadius", elbowRadius, (float, int))
        check_type("pipe1outletFaceName", pipe1outletFaceName, str)
        check_type("pipe2inletFaceName", pipe2inletFaceName, str)
        check_value("pipe1outletFaceName", pipe1outletFaceName, _FACE_NAME_TYPES)
        check_value("pipe2inletFaceName", pipe2inletFaceName, _FACE_NAME_TYPES)

        # Extract faces
        pipe1outletFace = pipe1.getFace(pipe1outletFaceName)
        pipe2inletFace = pipe2.getFace(pipe2inletFaceName)

        # Compute barycenters
        pipe1Outlet = pipe1.getFaceBarycenter(pipe1outletFace)
        pipe2Inlet = pipe2.getFaceBarycenter(pipe2inletFace)

        # Compute directions
        dir1 = pipe1.getFaceNormal(pipe1outletFaceName)
        dir2 = -pipe2.getFaceNormal(pipe2inletFaceName)

        # Check alignment using dot product
        dir1.normalize()
        dir2.normalize()
        dot_product = dir1.dot(dir2)
        aligned = abs(dot_product - 1.0) < 1e-6 or abs(dot_product + 1.0) < 1e-6

        if aligned:
            # === TRANSLATIONAL CYCLIC ===
            separationVector = pipe2Inlet - pipe1Outlet

            inlet = FaceCyclic(
                name=f"{pipe1.name}_outletAMI",
                neighbourPatch=f"{pipe2.name}_inletAMI",
                transform="translational",
                separationVector=separationVector
            )
            inlet.addSubFace(pipe1outletFace)

            outlet = FaceCyclic(
                name=f"{pipe2.name}_inletAMI",
                neighbourPatch=f"{pipe1.name}_outletAMI",
                transform="translational",
                separationVector=-separationVector
            )
            outlet.addSubFace(pipe2inletFace)

        else:
            # === ROTATIONAL CYCLIC ===
            rotationAxis1 = dir2.cross(dir1)
            rotationAxis2 = dir1.cross(dir2)
            rotationAxis1.normalize()
            rotationAxis2.normalize()

            if elbowRadius == 0:
                rotationCentre = pipe1Outlet
            else:
                R1 = dir1.cross(rotationAxis1)
                R1.normalize(elbowRadius)
                rotationCentre = pipe1Outlet + R1

            inlet = FaceCyclic(
                name=f"{pipe1.name}_outletAMI",
                neighbourPatch=f"{pipe2.name}_inletAMI",
                transform="rotational",
                rotationAxis=rotationAxis1,
                rotationCentre=rotationCentre
            )
            inlet.addSubFace(pipe1outletFace)

            outlet = FaceCyclic(
                name=f"{pipe2.name}_inletAMI",
                neighbourPatch=f"{pipe1.name}_outletAMI",
                transform="rotational",
                rotationAxis=rotationAxis2,
                rotationCentre=rotationCentre
            )
            outlet.addSubFace(pipe2inletFace)

        self.addBoundary(inlet)
        self.addBoundary(outlet)

    def createPipeBranch(
        self,
        originPipe: Block,
        branches: list[dict],
        isAddBoundaryConditions: bool = False,
        originOutletFaceName: str = 'top'
    ) -> list[Block]:
        """
        Create N branching pipes from a single origin pipe.

        Parameters
        ----------
        originPipe : Block
            The pipe block from which the branches originate.
        branches : list of dict
            Each dict must contain:
                - name: str
                - direction: Vector
                - length: float
                - equivalentHydraulicDiameter: float
                - elbowRadius: float
                - n: int (optional, default 1)
        isAddBoundaryConditions : bool
            Whether to add wall boundary conditions and connect pipes.
        originOutletFaceName : str
            Face name of the origin pipe to branch from.

        Returns
        -------
        list[Block]
            The list of newly created pipe blocks.
        """
        createdPipes = []

        for branch in branches:
            name = branch["name"]
            direction = branch["direction"]
            length = branch["length"]
            dh = branch["equivalentHydraulicDiameter"]  
            elbow = branch["elbowRadius"]             
            n = branch.get("n", 1)

            newPipe = self.addPipe1DFromDirection(
                name=name,
                originPosition=originPipe,
                direction=direction,
                length=length,
                equivalentHydraulicDiameter=dh,
                n=n,
                elbowRadius=elbow,
                isAddBoundaryConditions=False,
                originPositionOutletFaceName=originOutletFaceName
            )

            createdPipes.append(newPipe)

            if isAddBoundaryConditions:
                self.pipeWallBC.addSubFace(newPipe.leftFace())
                self.pipeWallBC.addSubFace(newPipe.rightFace())
                self.pipeWallBC.addSubFace(newPipe.frontFace())
                self.pipeWallBC.addSubFace(newPipe.backFace())

                self.connectPipesCustomNames(
                    pipe1=originPipe,
                    pipe2=newPipe,
                    elbowRadius=elbow,  
                    pipe1outletFaceName=originOutletFaceName,
                    pipe2inletFaceName='bottom',
                    customOutletName=f"{originPipe.name}to{name}_outletAMI",
                    customInletName=f"{name}_inletAMI"
                )

        return createdPipes



    def reconnectBranchedToFinal(
        self,
        remainingBranches: list,
        firstIntermediate: 'Block',
        finalPipe: 'Block',
        intermediatePipeSpecs: list,
        elbowRadius: float = 0.0
    ) -> list[Block]:
        """
        Reconnect branched pipes to the final pipe via new intermediate pipes.
        Returns the list of newly created intermediate pipes.

        Parameters
        ----------
        remainingBranches : list[Block]
            List of branch pipes that still need to reach the final pipe.
        firstIntermediate : Block
            The first intermediate pipe (already created by the user).
        finalPipe : Block
            The final pipe (already created by the user).
        intermediatePipeSpecs : list[dict]
            List of specs for the intermediate pipes to create for each remaining branch:
            [
                {"name": str, "length": float, "direction": Vector, "equivalentHydraulicDiameter": float, "n": int, "elbowRadius": float},
                ...
            ]
        elbowRadius : float
            Default elbow radius.
        """
        createdIntermediates = []

        # Step 1: connect first intermediate → final pipe
        self.connectPipesCustomNames(
            pipe1=firstIntermediate,
            pipe2=finalPipe,
            customOutletName=f"{firstIntermediate.name}_outletAMI",
            customInletName=f"{finalPipe.name}from{firstIntermediate.name}_inletAMI",
            pipe1outletFaceName='top',
            pipe2inletFaceName='bottom',
            elbowRadius=elbowRadius
        )

        # Step 2: add wall BCs to final pipe
        self.pipeWallBC.addSubFace(finalPipe.leftFace())
        self.pipeWallBC.addSubFace(finalPipe.rightFace())
        self.pipeWallBC.addSubFace(finalPipe.frontFace())
        self.pipeWallBC.addSubFace(finalPipe.backFace())

        # Step 3: loop over remaining branches
        for branch, spec in zip(remainingBranches, intermediatePipeSpecs):
            # create intermediate pipe between branch and final
            intermediate = self.addPipe1DFrom2Points(
                name=spec["name"],
                originPosition=branch,
                finalPosition=finalPipe,
                equivalentHydraulicDiameter=spec["equivalentHydraulicDiameter"],
                n=spec["n"],
                elbowRadius=spec.get("elbowRadius", elbowRadius),
                isAddBoundaryConditions=False
            )

            createdIntermediates.append(intermediate)

            # normal connect branch → intermediate
            self.connectPipes(branch, intermediate, elbowRadius=elbowRadius)

            # custom connect intermediate → final pipe
            self.connectPipesCustomNames(
                pipe1=intermediate,
                pipe2=finalPipe,
                customOutletName=f"{intermediate.name}_outletAMI",
                customInletName=f"{finalPipe.name}from{intermediate.name}_inletAMI",
                pipe1outletFaceName='top',
                pipe2inletFaceName='bottom',
                elbowRadius=elbowRadius
            )

            # add wall BCs to intermediate pipe
            self.pipeWallBC.addSubFace(intermediate.leftFace())
            self.pipeWallBC.addSubFace(intermediate.rightFace())
            self.pipeWallBC.addSubFace(intermediate.frontFace())
            self.pipeWallBC.addSubFace(intermediate.backFace())

        return createdIntermediates





    def connectPipesCustomNames(
        self,
        pipe1: Block,
        pipe2: Block,
        elbowRadius: float = 0,
        pipe1outletFaceName: str = 'top',
        pipe2inletFaceName: str = 'bottom',
        customOutletName: str = None,  
        customInletName: str = None   
    ) -> None:
        """
        Connect two pipes with custom cyclic patch names to avoid duplication.

        CHANGES:
        - Accepts both customOutletName and customInletName
        - Applies patch names to the correct owning blocks
        """
        check_type("pipe1", pipe1, Block)
        check_type("pipe2", pipe2, Block)
        check_type("elbowRadius", elbowRadius, (float, int))
        check_type("pipe1outletFaceName", pipe1outletFaceName, str)
        check_type("pipe2inletFaceName", pipe2inletFaceName, str)
        check_value("pipe1outletFaceName", pipe1outletFaceName, _FACE_NAME_TYPES)
        check_value("pipe2inletFaceName", pipe2inletFaceName, _FACE_NAME_TYPES)

        pipe1outletFace = pipe1.getFace(pipe1outletFaceName)
        pipe2inletFace = pipe2.getFace(pipe2inletFaceName)

        pipe1Outlet = pipe1.getFaceBarycenter(pipe1outletFace)
        pipe2Inlet = pipe2.getFaceBarycenter(pipe2inletFace)

        dir1 = pipe1.getFaceNormal(pipe1outletFaceName)
        dir2 = -pipe2.getFaceNormal(pipe2inletFaceName)

        dir1.normalize()
        dir2.normalize()
        dot_product = dir1.dot(dir2)
        aligned = abs(dot_product - 1.0) < 1e-6 or abs(dot_product + 1.0) < 1e-6

        if aligned:
            separationVector = pipe2Inlet - pipe1Outlet

            inlet = FaceCyclic(
                name=customOutletName, 
                neighbourPatch=customInletName,
                transform="translational",
                separationVector=separationVector
            )
            inlet.addSubFace(pipe1outletFace)

            outlet = FaceCyclic(
                name=customInletName,  
                neighbourPatch=customOutletName,
                transform="translational",
                separationVector=-separationVector
            )
            outlet.addSubFace(pipe2inletFace)

        else:
            rotationAxis1 = dir2.cross(dir1)
            rotationAxis2 = dir1.cross(dir2)
            rotationAxis1.normalize()
            rotationAxis2.normalize()

            if elbowRadius == 0:
                rotationCentre = pipe1Outlet
            else:
                R1 = dir1.cross(rotationAxis1)
                R1.normalize(elbowRadius)
                rotationCentre = pipe1Outlet + R1

            inlet = FaceCyclic(
                name=customOutletName, 
                neighbourPatch=customInletName,
                transform="rotational",
                rotationAxis=rotationAxis1,
                rotationCentre=rotationCentre
            )
            inlet.addSubFace(pipe1outletFace)

            outlet = FaceCyclic(
                name=customInletName, 
                neighbourPatch=customOutletName,
                transform="rotational",
                rotationAxis=rotationAxis2,
                rotationCentre=rotationCentre
            )
            outlet.addSubFace(pipe2inletFace)

        self.addBoundary(inlet)
        self.addBoundary(outlet)


    def createPipeCylindricalManifoldAlongZ(
            self,
            name: str,
            nEntries: int,
            innerRadius: float,
            outerRadius: float,
            equivalentHydraulicDiameter: float,
            lowZ: float=None,
            highZ: float=None,
            angleStart: float=0,
            x: float=0,
            y: float=0,
            nr: int=1,
            nt: int=1,
            isAddGap: bool=False
        ) -> list[Block]:
        """
        Create a cylindrical ring with `nEntries` for pipes. The face to attach
        the pipes is the `front` face. The pipe size is computed using the
        `equivalentHydraulicDiameter`. Between each pipe entry is an
        intermediate arc with `nt` azimutal cells discretization.

        Parameters
        ----------
        equivalentHydraulicDiameter : float
            Equivalent hydraulic diameter of a cylindrical pipe. Size of the
            square computed to conserve cross-section flow area as:

                `side = sqrt(pi) * eqDh`

        angleStart : float
            Angle start in degrees (default `0`).
        nr : int
            Number of cells radially (default `1`).
        nt : int
            Number of cells azimutally for the intermediate arcs (default `1`).
        isAddGap : bool
            If `True`, add a gap at regular interval between the pipes entry of
            the equivalent size of the eqDh (default `False`). If `False`, stich
            all the blocks without gaps in between.

        Return
        ------
        Return a list Block with the first element being the first pipe entry
        block directed toward the `+y` direction (or offset by `angleStart`).
        """
        check_type("name", name, str)
        check_type("nEntries", nEntries, int)
        check_type("innerRadius", innerRadius, (float, int))
        check_type("outerRadius", outerRadius, (float, int))
        check_type("equivalentHydraulicDiameter", equivalentHydraulicDiameter, (float, int))
        check_type("lowZ", lowZ, (float, int), none_ok=True)
        check_type("highZ", highZ, (float, int), none_ok=True)
        check_type("angleStart", angleStart, (float, int))
        check_type("x", x, (float, int))
        check_type("y", y, (float, int))
        check_type("nr", nr, int)
        check_type("nt", nt, int)
        if ((lowZ is None and highZ is None) or (lowZ is not None and highZ is not None)):
            msg = "Provide either 'lowZ' or 'highZ'"
            raise ValueError(msg)

        halfDh = np.sqrt(np.pi) * equivalentHydraulicDiameter/2

        anglePipe = 2*np.arcsin(np.sqrt(np.pi)*equivalentHydraulicDiameter / (2*outerRadius)) * 180/np.pi
        angleSection = (360 - nEntries*anglePipe)/nEntries

        if (highZ is None):
            highZ = lowZ + 2*halfDh
        elif (lowZ is None):
            lowZ = highZ - 2*halfDh

        if (isAddGap):
            blocks = []
            for i in range(nEntries):
                arc = self.createRingSectorAlongZ(
                    name=f"{name}_pipe{i}",
                    innerRadius=innerRadius, outerRadius=outerRadius,
                    angleStart=angleStart - anglePipe/2 + i*360/nEntries, angleArc=anglePipe,
                    lowZ=lowZ, highZ=highZ,
                    nr=nr, nt=1, nz=1,
                    x=x, y=y
                )
                blocks.append(arc)

            angleInnerPipe = 2*np.arcsin(np.sqrt(np.pi)*equivalentHydraulicDiameter / (2*innerRadius)) * 180/np.pi
            angleSection = (360/nEntries - anglePipe - angleInnerPipe)/2

            for arc in blocks[:nEntries]:
                arc1 = self.extrudeNormalArc(
                    arc, 'right', name, angleSection, arcCenter=Vector(x, y, 0), nt=nt
                )
                arc2 = self.extrudeNormalArc(
                    arc, 'left', name, -angleSection, arcCenter=Vector(x, y, 0), nt=nt
                )
                blocks.append(arc1)
                blocks.append(arc2)

            return(blocks)

        # No gaps
        firstArc = self.createRingSectorAlongZ(
            name=f"{name}_pipe0",
            innerRadius=innerRadius, outerRadius=outerRadius,
            angleStart=angleStart - anglePipe/2, angleArc=anglePipe,
            lowZ=lowZ, highZ=highZ,
            nr=nr, nt=1, nz=1,
            x=x, y=y
        )
        blocks = [firstArc]
        for i in range(nEntries-1):
            arc1 = self.extrudeNormalArc(
                blocks[-1], 'right', name, angleSection, arcCenter=Vector(x, y, 0), nt=nt
            )
            arc2 = self.extrudeNormalArc(
                arc1, 'right', f"{name}_pipe{i+1}", anglePipe, arcCenter=Vector(x, y, 0), nt=1
            )
            blocks.append(arc1)
            blocks.append(arc2)

        # Last section to close the loop
        lastArc = self.addRight(blocks[-1], name, [
            firstArc.points[0], firstArc.points[3],
            firstArc.points[4], firstArc.points[7]
        ], nx=nt)
        lastArc.addEdgeArc(0, 1, x, y, isOrigin=True)
        lastArc.addEdgeArc(2, 3, x, y, isOrigin=True)
        lastArc.addEdgeArc(4, 5, x, y, isOrigin=True)
        lastArc.addEdgeArc(6, 7, x, y, isOrigin=True)
        blocks.append(lastArc)

        return(blocks)


    def hexagonalLatticeAssembly(
            self,
            pitch: float,
            lattice: str, latticeNXY: int,
            zmin: float, zmax: float,
            wrapperFlatToFlat: float,
            xAssembly: float=0, yAssembly: float=0,
            nrEdge: int=3,
            ntEdge: int=4,
            nzEdge: int=2,
            funcElementGeneratorPin: dict={}
        ):
        """
        Parameters
        ----------
        funcElementGeneratorPin : dict[str, (lambda name, x, y)]
            Dictionnary of generator functions
        """
        for pinKey, funcGen in funcElementGeneratorPin.items():
            self.latticePlacement(
                funcElementGenerator=funcGen,
                lattice=lattice,
                latticeType='hexagon',
                pitch=pitch,
                nx=latticeNXY,
                ny=latticeNXY,
                x=xAssembly, y=yAssembly,
                elementsToPlace=pinKey
            )

        def angle(x_, y_):
            theta = np.arctan2(y_, x_)  # Angle in radians, range [-π, π]
            if theta < 0:
                theta += 2 * np.pi  # Normalize to [0, 2π)

            quantized = round(theta / (np.pi / 3)) * (np.pi / 3)
            return quantized


        distanceCenterToEdge = 0.5 * (wrapperFlatToFlat - (latticeNXY + (latticeNXY-2-1)/2)*pitch/np.sqrt(3))

        self.latticePlacement(
            funcElementGenerator=lambda name, x, y: self.createEdgeHexagonPrismAlongZ(
                "edge",
                zmin=zmin, zmax=zmax,
                pitch=pitch,
                distanceCenterToEdge=distanceCenterToEdge,
                edgeFaceOrientation=angle(x-xAssembly, y-yAssembly),
                x=x, y=y,
                nr=nrEdge, nt=ntEdge, nz=nzEdge,
                isAddBoundaryConditions=True
            ),
            lattice=lattice,
            latticeType='hexagon',
            pitch=pitch,
            nx=latticeNXY,
            ny=latticeNXY,
            x=xAssembly, y=yAssembly,
            elementsToPlace='E'
        )

        self.latticePlacement(
            funcElementGenerator=lambda name, x, y: self.createCornerHexagonPrismAlongZ(
                "corner",
                zmin=zmin, zmax=zmax,
                pitch=pitch,
                distanceCenterToEdge=distanceCenterToEdge,
                edgeFaceOrientation=np.arctan2(y-yAssembly, x-xAssembly)-np.pi/6,
                x=x, y=y,
                nr=nrEdge, nt=ntEdge, nz=nzEdge,
                isAddBoundaryConditions=True
            ),
            lattice=lattice,
            latticeType='hexagon',
            pitch=pitch,
            nx=latticeNXY,
            ny=latticeNXY,
            x=xAssembly, y=yAssembly,
            elementsToPlace='C'
        )


    def getNeiboursByIndex(
            self,
            lattice: str,
            idx: int,
            latticeType: str,
            nx: int,
            ny: int
        ) -> list[int]:
        """
        Used in lattice mesh generation
        """
        check_type("latticeType", latticeType, str)
        check_value("latticeType", latticeType, _LATTICE_TYPES)

        def dist(x0, y0, x1, y1):
            return((x0-x1)**2 + (y0-y1)**2)

        flatLattice = (' '.join(lattice.split("\n"))).split()

        coords = self.computeCoordinates(latticeType, 1, nx, ny)

        res = []
        for i, element, coord in zip(range(len(flatLattice)), flatLattice, coords):
            d = dist(coords[idx][0], coords[idx][1], coord[0], coord[1])
            if (element != "0" and 0 < d and d <= 1+1e-6):
                res.append(i)

        return(res)


    def latticePlacement(
            self,
            funcElementGenerator,
            lattice: str,
            nx: int,
            ny: int,
            pitch: float,
            latticeType: str,
            flatToFlatDirection: str='y',
            x: float=0,
            y: float=0,
            elementsToPlace: list[str] | str=None,
            isMergePatches: bool=False
        ):
        """
        Place the elementary submesh generated by funcElementGenerator using a
        lattice map.

        The lattice map can be a regular MCNP/Serpent format.
        Recommended use the Honeycomb tool
        https://foam-for-nuclear.gitlab.io/honeycomb/ to translate from OpenMC
        or CASMO.

        Usage ::

            # Define the lattice
            lattice = "0 0 W F W
             0 F F F F
              W F F F W
               F F F F 0
                W F W 0 0"

            mesh = ffn.BlockMesh()

            mesh.latticePlacement(
                funcElementGenerator=lambda name, x, y: mesh.createHexagonPrismAlongZ(
                    name=name,
                    zmin=0,
                    zmax=0.5,
                    pitch=latticePitch,
                    x=x, y=y,
                    nr=1,
                    nt=1,
                    nz=3,
                    isAddBoundaryConditions=True
                ),
                lattice=lattice,
                latticeType="hexagon",
                nx=5, ny=5, pitch=latticePitch,
                elementsToPlace=["F", "W"],
                isMergePatches=True
            )

        Parameters
        ----------
        funcElementGenerator : lambda name, x, y: None
            Generator function that create the elementary mesh to be replicated
        lattice : str
            Lattice map in text format, see Usage
        nx : int
            Number of rows in the lattice
        ny : int
            Number of columns in the lattice
        pitch : float
            Pitch of the lattice
        latticeType : str
            Lattice type (`hexagon` or `square`)
        flatToFlatDirection : str
            Flat to flat alignment direction ('x', 'y'), (default `'y'`).
        x : float
            Offset in the X-direction (default `0`).
        y : float
            Offset in the Y-direction (default `0`).
        elementsToPlace : list[str] | str | None
            Elements present in the lattice to be placed. If None (default) is
            found, place all the non-0 elements in the lattice.
        isMergePatches : bool
            Merge all the internal patches that might overlap during the
            placement (default `False`).
        """
        latticeElements = self.getLatticeAsList(
            lattice=lattice,
            nx=nx, ny=ny,
            pitch=pitch,
            latticeType=latticeType,
            flatToFlatDirection=flatToFlatDirection,
            x=x, y=y,
            elementsToPlace=elementsToPlace
        )

        for elementName, coord0, coord1 in latticeElements:
            funcElementGenerator(f"{elementName}", coord0, coord1)

        if (isMergePatches):
            self.addMergePatchPairs()


    def fillLatticeRingGap(
            self,
            name: str,
            ringRadius: float,
            nxLat: int,
            nyLat: int,
            pitch: float,
            latticeType: str,
            nxBlock: int,
            nyBlock: int,
            nzBlock: int,
            zmin: float,
            zmax: float,
            flatToFlatDirection: str='y',
            x: float=0,
            y: float=0,
            isAddBoundaryConditions: bool=False
        ):
        """
        Create cutted blocks that cross an infinite cylinder along the Z-axis.
        This function doesn't need to know the lattice because it assumes that
        the lattice expands the maximum space inside a circle. The blocks
        generated are just here to fill the gap between the lattice and the
        cylinder.

        Parameters
        ----------
        nxLat : int
            Number of rows in the lattice.
        nyLat : int
            Number of columns in the lattice.
        pitch : float
            Pitch of the lattice.
        latticeType : str {"hexagon", "square"}
            Lattice type.
        nxBlock : int
            Number of element in each blocks in the X-direction.
        nyBlock : int
            Number of element in each blocks in the Y-direction.
        nzBlock : int
            Number of element in each blocks in the Z-direction.
        zmin : float
            Minumum Z value.
        zmax : float
            Maximum Z value.
        flatToFlatDirection : str {"x", "y"}
            Flat to flat alignment direction (default `y`).
        x : float
            Offset in the X-direction (default `0`).
        y : float
            Offset in the Y-direction (default `0`).
        """
        check_type("latticeType", latticeType, str)
        check_value("latticeType", latticeType, _LATTICE_TYPES)

        coords = self.computeCoordinates(
            latticeType, pitch, nxLat+2, nyLat+2, xoffset=x, yoffset=y
        )

        isInside = lambda x, y: np.sqrt(x**2 + y**2) <= ringRadius

        if (isAddBoundaryConditions):
            idxFace = len(self.faces)

            topFace = Face(f"{name}Top_{idxFace}")
            botFace = Face(f"{name}Bottom_{idxFace}")
            wallFace = Face(f"{name}Wall_{idxFace}")

            self.addBoundary(topFace)
            self.addBoundary(botFace)
            self.addBoundary(wallFace)

        if (latticeType == 'square'):
            blockCorners = []
            cornerInsideRing = []

            # Create coords list
            for coord in coords:
                if (flatToFlatDirection == 'y'):
                    blockCorners.append((
                        coord[0]-pitch/2, # xmin
                        coord[0]+pitch/2, # xmax
                        coord[1]-pitch/2, # ymin
                        coord[1]+pitch/2, # ymax
                    ))
                elif (flatToFlatDirection == 'x'):
                    blockCorners.append((
                        coord[1]-pitch/2, # xmin
                        coord[1]+pitch/2, # xmax
                        coord[0]-pitch/2, # ymin
                        coord[0]+pitch/2, # ymax
                    ))

            # Check block corner inside the ring
            for xmin, xmax, ymin, ymax in blockCorners:
                cornerInsideRing.append((
                    isInside(xmin, ymin),
                    isInside(xmax, ymin),
                    isInside(xmax, ymax),
                    isInside(xmin, ymax),
                ))

            # Filter to block that have at least one corner inside but not all
            # corners
            chopedBlocks = [
                block for block, test in zip(blockCorners, cornerInsideRing)
                if any(test) and not all(test)
            ]

            for block, test in zip(blockCorners, cornerInsideRing):
                nCornerInside = len([val for val in test if val])

                # Create a wedge
                if (nCornerInside == 1):
                    # Back right
                    if (test == (True, False, False, False)):
                        xmin, xmax, ymin, ymax = block

                        newBlock = self.createBlock(name, [
                            Point(xmin, ymin, zmin),
                            Point(x + np.sqrt(ringRadius**2 - (ymin-x)**2), ymin, zmin),
                            Point(xmin, y + np.sqrt(ringRadius**2 - (xmin-x)**2), zmin),
                            Point(xmin, ymin, zmin),
                            Point(xmin, ymin, zmax),
                            Point(x + np.sqrt(ringRadius**2 - (ymin-x)**2), ymin, zmax),
                            Point(xmin, y + np.sqrt(ringRadius**2 - (xmin-x)**2), zmax),
                            Point(xmin, ymin, zmax),
                        ], nx=nxBlock, ny=nyBlock, nz=nzBlock)

                        newBlock.addEdgeArc(1, 2, x, y, isOrigin=True)
                        newBlock.addEdgeArc(5, 6, x, y, isOrigin=True)

                        if (isAddBoundaryConditions):
                            idxFace = len(self.faces)

                            wallFace.addSubFace(newBlock.rightFace())

                            wall1 = Face(f"{name}WallCorner1_{idxFace}")
                            wall1.addSubFace(newBlock.frontFace())
                            wall2 = Face(f"{name}WallCorner2_{idxFace}")
                            wall2.addSubFace(newBlock.backFace())

                            self.addBoundary(wall1)
                            self.addBoundary(wall2)

                    # Front left
                    elif (test == (False, False, True, False)):
                        xmin, xmax, ymin, ymax = block

                        newBlock = self.createBlock(name, [
                            Point(x - np.sqrt(ringRadius**2 - (ymax-x)**2), ymax, zmin),
                            Point(xmax, y - np.sqrt(ringRadius**2 - (xmax-x)**2), zmin),
                            Point(xmax, ymax, zmin),
                            Point(xmax, ymax, zmin),
                            Point(x - np.sqrt(ringRadius**2 - (ymax-x)**2), ymax, zmax),
                            Point(xmax, y - np.sqrt(ringRadius**2 - (xmax-x)**2), zmax),
                            Point(xmax, ymax, zmax),
                            Point(xmax, ymax, zmax),
                        ], nx=nxBlock, ny=nyBlock, nz=nzBlock)

                        newBlock.addEdgeArc(0, 1, x, y, isOrigin=True)
                        newBlock.addEdgeArc(4, 5, x, y, isOrigin=True)

                        if (isAddBoundaryConditions):
                            idxFace = len(self.faces)

                            wallFace.addSubFace(newBlock.frontFace())

                            wall1 = Face(f"{name}WallCorner1_{idxFace}")
                            wall1.addSubFace(newBlock.leftFace())
                            wall2 = Face(f"{name}WallCorner2_{idxFace}")
                            wall2.addSubFace(newBlock.rightFace())

                            self.addBoundary(wall1)
                            self.addBoundary(wall2)

                    # Front right
                    elif (test == (False, False, False, True)):
                        xmin, xmax, ymin, ymax = block

                        newBlock = self.createBlock(name, [
                            Point(xmin, y - np.sqrt(ringRadius**2 - (xmin-x)**2), zmin),
                            Point(x + np.sqrt(ringRadius**2 - (ymax-x)**2), ymax, zmin),
                            Point(xmin, ymax, zmin),
                            Point(xmin, ymax, zmin),
                            Point(xmin, y - np.sqrt(ringRadius**2 - (xmin-x)**2), zmax),
                            Point(x + np.sqrt(ringRadius**2 - (ymax-x)**2), ymax, zmax),
                            Point(xmin, ymax, zmax),
                            Point(xmin, ymax, zmax),
                        ], nx=nxBlock, ny=nyBlock, nz=nzBlock)

                        newBlock.addEdgeArc(0, 1, x, y, isOrigin=True)
                        newBlock.addEdgeArc(4, 5, x, y, isOrigin=True)

                        if (isAddBoundaryConditions):
                            idxFace = len(self.faces)

                            wallFace.addSubFace(newBlock.frontFace())

                            wall1 = Face(f"{name}WallCorner1_{idxFace}")
                            wall1.addSubFace(newBlock.leftFace())
                            wall2 = Face(f"{name}WallCorner2_{idxFace}")
                            wall2.addSubFace(newBlock.rightFace())

                            self.addBoundary(wall1)
                            self.addBoundary(wall2)

                    # Back left
                    elif (test == (False, True, False, False)):
                        xmin, xmax, ymin, ymax = block

                        newBlock = self.createBlock(name, [
                            Point(x - np.sqrt(ringRadius**2 - (ymin-x)**2), ymin, zmin),
                            Point(xmax, ymin, zmin),
                            Point(xmax, ymin, zmin),
                            Point(xmax, y + np.sqrt(ringRadius**2 - (xmax-x)**2), zmin),
                            Point(x - np.sqrt(ringRadius**2 - (ymin-x)**2), ymin, zmax),
                            Point(xmax, ymin, zmax),
                            Point(xmax, ymin, zmax),
                            Point(xmax, y + np.sqrt(ringRadius**2 - (xmax-x)**2), zmax),
                        ], nx=nxBlock, ny=nyBlock, nz=nzBlock)

                        newBlock.addEdgeArc(0, 3, x, y, isOrigin=True)
                        newBlock.addEdgeArc(4, 7, x, y, isOrigin=True)

                        if (isAddBoundaryConditions):
                            idxFace = len(self.faces)

                            wallFace.addSubFace(newBlock.leftFace())

                            wall1 = Face(f"{name}WallCorner1_{idxFace}")
                            wall1.addSubFace(newBlock.frontFace())
                            wall2 = Face(f"{name}WallCorner2_{idxFace}")
                            wall2.addSubFace(newBlock.backFace())

                            self.addBoundary(wall1)
                            self.addBoundary(wall2)

                    if (isAddBoundaryConditions):
                        topFace.addSubFace(newBlock.topFace())
                        botFace.addSubFace(newBlock.bottomFace())

                # Deform a cube
                elif (nCornerInside == 2):
                    # Back
                    if (test == (True, True, False, False)):
                        xmin, xmax, ymin, ymax = block

                        newBlock = self.createBlock(name, [
                            Point(xmin, ymin, zmin),
                            Point(xmax, ymin, zmin),
                            Point(xmax, y + np.sqrt(ringRadius**2 - (xmax-x)**2), zmin),
                            Point(xmin, y + np.sqrt(ringRadius**2 - (xmin-x)**2), zmin),
                            Point(xmin, ymin, zmax),
                            Point(xmax, ymin, zmax),
                            Point(xmax, y + np.sqrt(ringRadius**2 - (xmax-x)**2), zmax),
                            Point(xmin, y + np.sqrt(ringRadius**2 - (xmin-x)**2), zmax),
                        ], nx=nxBlock, ny=nyBlock, nz=nzBlock)

                        newBlock.addEdgeArc(2, 3, x, y, isOrigin=True)
                        newBlock.addEdgeArc(6, 7, x, y, isOrigin=True)

                        if (isAddBoundaryConditions):
                            idxFace = len(self.faces)

                            wallFace.addSubFace(newBlock.backFace())

                            wall1 = Face(f"{name}WallFront_{idxFace}")
                            wall1.addSubFace(newBlock.frontFace())
                            wall2 = Face(f"{name}WallRight_{idxFace}")
                            wall2.addSubFace(newBlock.rightFace())
                            wall3 = Face(f"{name}WallLeft_{idxFace}")
                            wall3.addSubFace(newBlock.leftFace())

                            self.addBoundary(wall1)
                            self.addBoundary(wall2)
                            self.addBoundary(wall3)

                    # Front
                    elif (test == (False, False, True, True)):
                        xmin, xmax, ymin, ymax = block

                        newBlock = self.createBlock(name, [
                            Point(xmin, y - np.sqrt(ringRadius**2 - (xmin-x)**2), zmin),
                            Point(xmax, y - np.sqrt(ringRadius**2 - (xmax-x)**2), zmin),
                            Point(xmax, ymax, zmin),
                            Point(xmin, ymax, zmin),
                            Point(xmin, y - np.sqrt(ringRadius**2 - (xmin-x)**2), zmax),
                            Point(xmax, y - np.sqrt(ringRadius**2 - (xmax-x)**2), zmax),
                            Point(xmax, ymax, zmax),
                            Point(xmin, ymax, zmax),
                        ], nx=nxBlock, ny=nyBlock, nz=nzBlock)

                        newBlock.addEdgeArc(0, 1, x, y, isOrigin=True)
                        newBlock.addEdgeArc(4, 5, x, y, isOrigin=True)

                        if (isAddBoundaryConditions):
                            idxFace = len(self.faces)

                            wallFace.addSubFace(newBlock.frontFace())

                            wall1 = Face(f"{name}WallBack_{idxFace}")
                            wall1.addSubFace(newBlock.backFace())
                            wall2 = Face(f"{name}WallLeft_{idxFace}")
                            wall2.addSubFace(newBlock.leftFace())
                            wall3 = Face(f"{name}WallRight_{idxFace}")
                            wall3.addSubFace(newBlock.rightFace())

                            self.addBoundary(wall1)
                            self.addBoundary(wall2)
                            self.addBoundary(wall3)

                    # Right
                    elif (test == (True, False, False, True)):
                        xmin, xmax, ymin, ymax = block

                        newBlock = self.createBlock(name, [
                            Point(xmin, ymin, zmin),
                            Point(x + np.sqrt(ringRadius**2 - (ymin-y)**2), ymin, zmin),
                            Point(x + np.sqrt(ringRadius**2 - (ymax-y)**2), ymax, zmin),
                            Point(xmin, ymax, zmin),
                            Point(xmin, ymin, zmax),
                            Point(x + np.sqrt(ringRadius**2 - (ymin-y)**2), ymin, zmax),
                            Point(x + np.sqrt(ringRadius**2 - (ymax-y)**2), ymax, zmax),
                            Point(xmin, ymax, zmax),
                        ], nx=nxBlock, ny=nyBlock, nz=nzBlock)

                        newBlock.addEdgeArc(1, 2, x, y, isOrigin=True)
                        newBlock.addEdgeArc(5, 6, x, y, isOrigin=True)

                        if (isAddBoundaryConditions):
                            idxFace = len(self.faces)

                            wallFace.addSubFace(newBlock.rightFace())

                            wall1 = Face(f"{name}WallFront_{idxFace}")
                            wall1.addSubFace(newBlock.frontFace())
                            wall2 = Face(f"{name}WallLeft_{idxFace}")
                            wall2.addSubFace(newBlock.leftFace())
                            wall3 = Face(f"{name}WallBack_{idxFace}")
                            wall3.addSubFace(newBlock.backFace())

                            self.addBoundary(wall1)
                            self.addBoundary(wall2)
                            self.addBoundary(wall3)

                    # Left
                    elif (test == (False, True, True, False)):
                        xmin, xmax, ymin, ymax = block

                        newBlock = self.createBlock(name, [
                            Point(x - np.sqrt(ringRadius**2 - (ymin-y)**2), ymin, zmin),
                            Point(xmax, ymin, zmin),
                            Point(xmax, ymax, zmin),
                            Point(x - np.sqrt(ringRadius**2 - (ymax-y)**2), ymax, zmin),
                            Point(x - np.sqrt(ringRadius**2 - (ymin-y)**2), ymin, zmax),
                            Point(xmax, ymin, zmax),
                            Point(xmax, ymax, zmax),
                            Point(x - np.sqrt(ringRadius**2 - (ymax-y)**2), ymax, zmax),
                        ], nx=nxBlock, ny=nyBlock, nz=nzBlock)

                        newBlock.addEdgeArc(0, 3, x, y, isOrigin=True)
                        newBlock.addEdgeArc(4, 7, x, y, isOrigin=True)

                        if (isAddBoundaryConditions):
                            idxFace = len(self.faces)

                            wallFace.addSubFace(newBlock.leftFace())

                            wall1 = Face(f"{name}WallFront_{idxFace}")
                            wall1.addSubFace(newBlock.frontFace())
                            wall2 = Face(f"{name}WallRight_{idxFace}")
                            wall2.addSubFace(newBlock.rightFace())
                            wall3 = Face(f"{name}WallBack_{idxFace}")
                            wall3.addSubFace(newBlock.backFace())

                            self.addBoundary(wall1)
                            self.addBoundary(wall2)
                            self.addBoundary(wall3)

                    if (isAddBoundaryConditions):
                        topFace.addSubFace(newBlock.topFace())
                        botFace.addSubFace(newBlock.bottomFace())

                # Double hex
                elif (nCornerInside == 3):
                    # Back Right
                    if (test == (True, True, False, True)):
                        xmin, xmax, ymin, ymax = block

                        x1 = xmax
                        y1 = y + np.sqrt(ringRadius**2 - (xmax-x)**2)
                        x2 = x + np.sqrt(ringRadius**2 - (ymax-y)**2)
                        y2 = ymax
                        xMid = (x1+x2)/2
                        yMid = y + np.sqrt(ringRadius**2 - (xMid-x)**2)

                        newBlock1 = self.createBlock(name, [
                            Point(xmin, ymin, zmin),
                            Point(xmax, ymin, zmin),
                            Point(x1, y1, zmin),
                            Point(xMid, yMid, zmin),
                            Point(xmin, ymin, zmax),
                            Point(xmax, ymin, zmax),
                            Point(x1, y1, zmax),
                            Point(xMid, yMid, zmax),
                        ], nx=nxBlock, ny=nyBlock, nz=nzBlock)

                        newBlock2 = self.createBlock(name, [
                            newBlock1.points[0],
                            newBlock1.points[3],
                            Point(x2, y2, zmin),
                            Point(xmin, ymax, zmin),
                            newBlock1.points[4],
                            newBlock1.points[7],
                            Point(x2, y2, zmax),
                            Point(xmin, ymax, zmax),
                        ], nx=nxBlock, ny=nyBlock, nz=nzBlock)

                        newBlock1.addEdgeArc(2, 3, x, y, isOrigin=True)
                        newBlock1.addEdgeArc(6, 7, x, y, isOrigin=True)
                        newBlock2.addEdgeArc(1, 2, x, y, isOrigin=True)
                        newBlock2.addEdgeArc(5, 6, x, y, isOrigin=True)

                        if (isAddBoundaryConditions):
                            idxFace = len(self.faces)

                            wallFace.addSubFace(newBlock1.backFace())
                            wallFace.addSubFace(newBlock2.rightFace())

                            # wall1 = Face(f"{name}WallFront_{idxFace}")
                            # wall1.addSubFace(newBlock.frontFace())
                            # wall2 = Face(f"{name}WallRight_{idxFace}")
                            # wall2.addSubFace(newBlock.rightFace())
                            # wall3 = Face(f"{name}WallBack_{idxFace}")
                            # wall3.addSubFace(newBlock.backFace())
                            # wall3 = Face(f"{name}WallBack_{idxFace}")
                            # wall3.addSubFace(newBlock.backFace())

                            # self.addBoundary(wall1)
                            # self.addBoundary(wall2)
                            # self.addBoundary(wall3)

                    # Front left
                    if (test == (False, True, True, True)):
                        xmin, xmax, ymin, ymax = block

                        x1 = xmin
                        y1 = y - np.sqrt(ringRadius**2 - (xmin-x)**2)
                        x2 = x - np.sqrt(ringRadius**2 - (ymin-y)**2)
                        y2 = ymin
                        xMid = (x1+x2)/2
                        yMid = y - np.sqrt(ringRadius**2 - (xMid-x)**2)

                        newBlock1 = self.createBlock(name, [
                            Point(x2, y2, zmin),
                            Point(xmax, ymin, zmin),
                            Point(xmax, ymax, zmin),
                            Point(xMid, yMid, zmin),
                            Point(x2, y2, zmax),
                            Point(xmax, ymin, zmax),
                            Point(xmax, ymax, zmax),
                            Point(xMid, yMid, zmax),
                        ], nx=nxBlock, ny=nyBlock, nz=nzBlock)

                        newBlock2 = self.createBlock(name, [
                            Point(x1, y1, zmin),
                            newBlock1.points[3],
                            newBlock1.points[2],
                            Point(xmin, ymax, zmin),
                            Point(x1, y1, zmax),
                            newBlock1.points[7],
                            newBlock1.points[6],
                            Point(xmin, ymax, zmax),
                        ], nx=nxBlock, ny=nyBlock, nz=nzBlock)

                        newBlock1.addEdgeArc(0, 3, x, y, isOrigin=True)
                        newBlock1.addEdgeArc(4, 7, x, y, isOrigin=True)
                        newBlock2.addEdgeArc(0, 1, x, y, isOrigin=True)
                        newBlock2.addEdgeArc(4, 5, x, y, isOrigin=True)

                        if (isAddBoundaryConditions):
                            idxFace = len(self.faces)

                            wallFace.addSubFace(newBlock1.leftFace())
                            wallFace.addSubFace(newBlock2.frontFace())

                    # Front right
                    if (test == (True, False, True, True)):
                        xmin, xmax, ymin, ymax = block

                        x1 = xmax
                        y1 = y - np.sqrt(ringRadius**2 - (xmax-x)**2)
                        x2 = x + np.sqrt(ringRadius**2 - (ymin-y)**2)
                        y2 = ymin
                        xMid = (x1+x2)/2
                        yMid = y - np.sqrt(ringRadius**2 - (xMid-x)**2)

                        newBlock1 = self.createBlock(name, [
                            Point(xmin, ymin, zmin),
                            Point(x2, y2, zmin),
                            Point(xMid, yMid, zmin),
                            Point(xmin, ymax, zmin),
                            Point(xmin, ymin, zmax),
                            Point(x2, y2, zmax),
                            Point(xMid, yMid, zmax),
                            Point(xmin, ymax, zmax),
                        ], nx=nxBlock, ny=nyBlock, nz=nzBlock)

                        newBlock2 = self.createBlock(name, [
                            newBlock1.points[2],
                            Point(x1, y1, zmin),
                            Point(xmax, ymax, zmin),
                            newBlock1.points[3],
                            newBlock1.points[6],
                            Point(x1, y1, zmax),
                            Point(xmax, ymax, zmax),
                            newBlock1.points[7],
                        ], nx=nxBlock, ny=nyBlock, nz=nzBlock)

                        newBlock1.addEdgeArc(1, 2, x, y, isOrigin=True)
                        newBlock1.addEdgeArc(5, 6, x, y, isOrigin=True)
                        newBlock2.addEdgeArc(0, 1, x, y, isOrigin=True)
                        newBlock2.addEdgeArc(4, 5, x, y, isOrigin=True)

                        if (isAddBoundaryConditions):
                            idxFace = len(self.faces)

                            wallFace.addSubFace(newBlock1.rightFace())
                            wallFace.addSubFace(newBlock2.frontFace())

                    # Back left
                    if (test == (True, True, True, False)):
                        xmin, xmax, ymin, ymax = block

                        x1 = xmin
                        y1 = y + np.sqrt(ringRadius**2 - (xmin-x)**2)
                        x2 = x - np.sqrt(ringRadius**2 - (ymax-y)**2)
                        y2 = ymax
                        xMid = (x1+x2)/2
                        yMid = y + np.sqrt(ringRadius**2 - (xMid-x)**2)

                        newBlock1 = self.createBlock(name, [
                            Point(xmin, ymin, zmin),
                            Point(xmax, ymin, zmin),
                            Point(xMid, yMid, zmin),
                            Point(x1, y1, zmin),
                            Point(xmin, ymin, zmax),
                            Point(xmax, ymin, zmax),
                            Point(xMid, yMid, zmax),
                            Point(x1, y1, zmax),
                        ], nx=nxBlock, ny=nyBlock, nz=nzBlock)

                        newBlock2 = self.createBlock(name, [
                            newBlock1.points[2],
                            newBlock1.points[1],
                            Point(xmax, ymax, zmin),
                            Point(x2, y2, zmin),
                            newBlock1.points[6],
                            newBlock1.points[5],
                            Point(xmax, ymax, zmax),
                            Point(x2, y2, zmax),
                        ], nx=nxBlock, ny=nyBlock, nz=nzBlock)

                        newBlock1.addEdgeArc(2, 3, x, y, isOrigin=True)
                        newBlock1.addEdgeArc(6, 7, x, y, isOrigin=True)
                        newBlock2.addEdgeArc(0, 3, x, y, isOrigin=True)
                        newBlock2.addEdgeArc(4, 7, x, y, isOrigin=True)

                        if (isAddBoundaryConditions):
                            idxFace = len(self.faces)

                            wallFace.addSubFace(newBlock1.backFace())
                            wallFace.addSubFace(newBlock2.leftFace())

                    if (isAddBoundaryConditions):
                        topFace.addSubFace(newBlock1.topFace())
                        topFace.addSubFace(newBlock2.topFace())
                        botFace.addSubFace(newBlock1.bottomFace())
                        botFace.addSubFace(newBlock2.bottomFace())

        elif (latticeType == 'hexagon'):
            msg = "Lattice type 'hexagon' not yet implemented in 'fillLatticeRingGap'"
            raise NotImplementedError(msg)

    def getOverlappingFaces(
            self,
            overlappingFaces: list[Face]=[],
            includeFacename: list[str]=[],
            excludeFacename: list[str]=[],
            description: str='overlappingFaces'
        ) -> list[tuple[Face]]:
        """
        Look for overlapping patches i.e exact one-to-one patch face matching.
        Be very careful with this method, it is extremely computational heavy
        both for Python and `blockMesh`.

        Parameters
        ----------
        overlappingFaces : list[Face]
            List of faces to not take into account in the search of overlapping
            faces.
        includeFacename : list[str]
            List of keyword include in face names to include the faces. For
            example (`faces = ['wall1', 'wall2', 'top']`,
            `includeFacename = ['wall']`, `results = ['wall1', 'wall2']`)
        excludeFacename : list[str]
            List of keyword include in face names to exclude the faces. For
            example (`faces = ['wall1', 'wall2', 'top']`,
            `excludeFacename = ['wall']`, `results = ['top']`)
        description : str
            Optional description during progress bar (default `overlappingFaces`).

        Returns
        -------
        List of face pairs
        """
        check_type("overlappingFaces", overlappingFaces, list)
        check_type("includeFacename", includeFacename, list)
        check_type("excludeFacename", excludeFacename, list)

        overlappingFaceNames = [face.name for face in overlappingFaces]

        # Remove already overlap faces
        faces = [
            (
                i,
                face,
                face.name in overlappingFaceNames # True means not look for it
            )
            for i, face in enumerate(self.faces)
        ]

        # Filter remaining faces to include
        if (len(includeFacename) != 0):
            faces = [
                (
                    i,
                    face,
                    isOverlapi if any([regex in face.name for regex in includeFacename]) else True
                )
                for i, face, isOverlapi in faces
            ]

        # Filter remaining faces to exclude
        if (len(excludeFacename) != 0):
            faces = [
                (
                    i,
                    face,
                    True if any([regex in face.name for regex in excludeFacename]) else isOverlapi
                )
                for i, face, isOverlapi in faces
            ]

        # Search
        coupledFaces = []
        for i, facei, isOverlapi in tqdm.tqdm(faces, desc=f"{self.region} {description}"):
            if (isOverlapi):
                continue

            isBreak = False
            for j, facej, isOverlapj in faces[i:]:
                if (isOverlapj or facei == facej):
                    continue

                # Loop over single subface for each patch, need at least one
                # subface to say it overlaps
                for faceii in facei.faces:
                    for facejj in facej.faces:
                        if (isEqualSubFace(faceii, facejj)):
                            coupledFaces.append((facei, facej))
                            # overlappingFaces.append(facei)
                            # overlappingFaces.append(facej)
                            faces[i] = (i, facei, True)
                            faces[j] = (j, facej, True)
                            isBreak = True
                            break
                    if (isBreak):
                        break
                if (isBreak):
                    break

        return(coupledFaces)


    def getStandaloneFaces(
            self,
            includeFacename: list[str]=[],
            excludeFacename: list[str]=[]
        ) -> list[Face]:
        """
        Extract all standalone faces.

        Parameters
        ----------
        includeFacename : list[str]
            List of keyword include in face names to include the faces. For
            example (`faces = ['wall1', 'wall2', 'top']`,
            `includeFacename = ['wall']`, `results = ['wall1', 'wall2']`)
        excludeFacename : list[str]
            List of keyword include in face names to exclude the faces. For
            example (`faces = ['wall1', 'wall2', 'top']`,
            `excludeFacename = ['wall']`, `results = ['top']`)

        Return
        ------
        List of standalone faces
        """
        check_type("includeFacename", includeFacename, list)
        check_type("excludeFacename", excludeFacename, list)

        mergedFaces  = [f1 for f1, f2 in self.mergePatchPairs] + [f2 for f1, f2 in self.mergePatchPairs]
        mergedFaces += [f1 for f1, f2 in self.baffleFaces] + [f2 for f1, f2 in self.baffleFaces]
        coupledFaces = self.getOverlappingFaces(
            overlappingFaces=mergedFaces,
            description='standaloneFaces'
        )

        coupledFaces = [f1.name for f1, f2 in coupledFaces] + [f2.name for f1, f2 in coupledFaces]
        coupledFaces += [face.name for face in mergedFaces]

        standaloneFaces = [face for face in self.faces if face.name not in coupledFaces]

        if (len(includeFacename) != 0):
            standaloneFaces = [face for face in standaloneFaces if any([regex in face.name for regex in includeFacename])]

        if (len(excludeFacename) != 0):
            standaloneFaces = [face for face in standaloneFaces if all([regex not in face.name for regex in excludeFacename])]

        return(standaloneFaces)


    def addBaffles(
            self,
            baffleName: str='baffle',
            includeFacename: list[str]=[],
            excludeFacename: list[str]=[]
        ) -> list[Face]:
        """
        Add faces into a baffle zone to merged in `createPatchDict`. Create 2
        new baffle faces named `<baffleName>0` and `<baffleName>1`

        Parameters
        ----------
        includeFacename : list[str]
            List of keyword include in face names to include the faces. For
            example (`faces = ['wall1', 'wall2', 'top']`,
            `includeFacename = ['wall']`, `results = ['wall1', 'wall2']`)
        excludeFacename : list[str]
            List of keyword include in face names to exclude the faces. For
            example (`faces = ['wall1', 'wall2', 'top']`,
            `excludeFacename = ['wall']`, `results = ['top']`)

        Returns
        -------
        List of baffle faces
        """
        check_type("includeFacename", includeFacename, list)
        check_type("excludeFacename", excludeFacename, list)

        mergedFaces  = [f1 for f1, f2 in self.mergePatchPairs] + [f2 for f1, f2 in self.mergePatchPairs]

        baffleFaces = self.getOverlappingFaces(
            overlappingFaces=mergedFaces,
            includeFacename=includeFacename,
            excludeFacename=excludeFacename,
            description='baffles'
        )
        self.baffleFaces += baffleFaces
        leftFaces  = [face0 for face0, face1 in baffleFaces]
        rightFaces = [face1 for face0, face1 in baffleFaces]

        baffle0 = Face(
            name=f"{baffleName}0",
            boundaryType="mappedWall",
            inGroups=['wall', 'mappedPatch', 'baffleFaces'],
            extraParameters={
                'sampleMode': 'nearestPatchFace',
                'samplePatch': f'{baffleName}1'
            }
        )
        baffle1 = Face(
            name=f"{baffleName}1",
            boundaryType="mappedWall",
            inGroups=['wall', 'mappedPatch', 'baffleFaces'],
            extraParameters={
                'sampleMode': 'nearestPatchFace',
                'samplePatch': f'{baffleName}0'
            }
        )

        for face in leftFaces:
            for face_i in face.faces:
                baffle0.addSubFace(face_i)
        for face in rightFaces:
            for face_i in face.faces:
                baffle1.addSubFace(face_i)

        leftNames = [face.name for face in leftFaces]
        rightNames = [face.name for face in rightFaces]

        self.faces = [face for face in self.faces if (face.name not in leftNames and face.name not in rightNames)]

        self.addBoundary(baffle0)
        self.addBoundary(baffle1)

        return(leftFaces + rightFaces)


    def addMergePatchPairs(
            self,
            includeFacename: list[str]=[],
            excludeFacename: list[str]=[],
        ) -> None:
        """
        Add merge patch pairs in `mergePatchPairs` dict of `blockMeshDict`. Find
        exact one-to-one patch pair matching for all faces.

        Usage ::

            # List of faces: `faces = ['wall1', 'wall2', 'top']`

            newMesh.addMergePatchPairs(
                includeFacename=['wall']
            )
            # >>> ['wall1', 'wall2']

            newMesh.addMergePatchPairs(
                excludeFacename=['wall']
            )
            # >>> ['top']

        If more than 2 faces need to be merged, please use `mergeBoundaryFaces`.

        Be very careful with this method, it is extremely computational heavy
        both for Python and `blockMesh`.

        Parameters
        ----------
        includeFacename : list[str]
            List of keyword include in face names to include the faces. If
            `includeFacename` is kept empty, look for all faces.
        excludeFacename : list[str]
            List of keyword include in face names to exclude the faces. If
            `excludeFacename` is kept empty, look for all faces.
        """
        check_type("includeFacename", includeFacename, list)
        check_type("excludeFacename", excludeFacename, list)

        mergedFaces  = [f1 for f1, f2 in self.mergePatchPairs] + [f2 for f1, f2 in self.mergePatchPairs]
        mergedFaces += [f1 for f1, f2 in self.baffleFaces] + [f2 for f1, f2 in self.baffleFaces]

        mergeableFaces = self.getOverlappingFaces(
            mergedFaces,
            includeFacename=includeFacename,
            excludeFacename=excludeFacename,
            description='mergePatchPairs'
        )
        self.mergePatchPairs += mergeableFaces


    def mergePatchesWithName(
            self,
            name: str,
            includeFacename: list[str]=[],
            excludeFacename: list[str]=[],
            patchType: str='patch',
            inGroups: list[str]=[],
            sampleMode: str=None,
            samplePatch: str=None
        ) -> Face:
        """
        Merge patches with name and create a new patch using `name`. Doesn't
        rely on `createPatch`.

        Parameters
        ----------
        name : str
            Name of the new patch
        patchType : str
            Patch type (e.g `patch`, `wall`, ...)
        """
        check_type("includeFacename", includeFacename, list)
        check_type("excludeFacename", excludeFacename, list)
        check_type("patchType", patchType, str)
        check_type("inGroups", inGroups, list)
        check_type("sampleMode", sampleMode, str, none_ok=True)
        check_type("samplePatch", samplePatch, str, none_ok=True)

        extraParameters = {}
        if (sampleMode is not None):
            extraParameters['sampleMode'] = sampleMode
        if (samplePatch is not None):
            extraParameters['samplePatch'] = samplePatch

        newFace = Face(
            name=name,
            boundaryType=patchType,
            inGroups=inGroups,
            extraParameters=extraParameters
        )

        strictExcludeFacenames = [face_i.name for face_i, face_j in self.baffleFaces]
        strictExcludeFacenames += [face_j.name for face_i, face_j in self.baffleFaces]
        strictExcludeFacenames += [face_i.name for face_i, face_j in self.mergePatchPairs]
        strictExcludeFacenames += [face_j.name for face_i, face_j in self.mergePatchPairs]

        facesToInclude = [
            face for face in self.faces
            if
                any([faceName in face.name for faceName in includeFacename])
                and
                all([faceName not in face.name for faceName in excludeFacename])
                and
                all([faceName != face.name for faceName in strictExcludeFacenames])
        ]

        for face in facesToInclude:
            for face_i in face.faces:
                newFace.addSubFace(face_i)

        faceNameToInclude = [face.name for face in facesToInclude]
        self.faces = [face for face in self.faces if (face.name not in faceNameToInclude)]

        self.addBoundary(newFace)

        return(newFace)


    def mergeBoundaryFaces(
            self,
            includeFacename: list[str],
            newBoundaryName: str,
            newBoundaryType: str='patch'
        ) -> Face:
        """
        Merge/aggregate multiple boundary face already defined in `boundary`
        dict of `blockMeshDict`. To be used before `addMergePatchPairs`. This
        method is directly modifying the `blockMeshDict` and doesn't require to
        use `mergePatchesWithName`.

        Parameters
        ----------
        includeFacename : list[str]
            List of keyword include in face names to include the faces.
        newBoundaryName : str
            New boundary condition name
        newBoundaryType : str
            New boundary condition type (default `patch`)
        """
        check_type("includeFacename", includeFacename, list)
        check_type("newBoundaryName", newBoundaryName, str)
        check_type("newBoundaryType", newBoundaryType, str)

        # Create a filtered list with only the faces that matches the
        # includeFacename list
        facesToInclude: list[Face] = [
            face for face in self.faces
            if len([regex for regex in includeFacename if regex in face.name]) > 0
        ]

        # Filter the self.faces list with only the faces that never matches the
        # includeFacename list
        self.faces = [
            face for face in self.faces
            if len([regex for regex in includeFacename if regex in face.name]) == 0
        ]

        # Create new Face and append to the list of boundary conditions
        newFace = Face(
            name=newBoundaryName,
            boundaryType=newBoundaryType
        )
        for face in facesToInclude:
            for facei in face.faces:
                newFace.addSubFace(facei)

        self.addBoundary(newFace)

        return(newFace)


    def mergePatchPairsByName(
            self,
            facename1: str,
            facename2: str
        ) -> None:
        """
        Merge already defined boundary faces using their names. Append to
        `mergePatchPairs` the merge of `facename1` and `facename2`. Doesn't
        create any new boundary.
        """
        check_type("facename1", facename1, str)
        check_type("facename2", facename2, str)

        self.mergePatchPairs.append((
            [face for face in self.faces if face.name == facename1][0],
            [face for face in self.faces if face.name == facename2][0]
        ))


    def duplicateBlock(
            self,
            targetBlock: Block,
            newName: str,
            xoffset: float=0,
            yoffset: float=0,
            zoffset: float=0
        ):
        """
        To be improved
        """
        newPoints = [
            Point(point.x+xoffset, point.y+yoffset, point.z+zoffset, isIndexed=False)
            for point in targetBlock.points
        ]

        newBlock = Block(
            newName,
            newPoints,
            nx=targetBlock.nx,
            ny=targetBlock.ny,
            nz=targetBlock.nz,
        )

        # for i, newPoint in enumerate(newPoints):
        #     isFound = False
        #     for block in self.blocks:
        #         for point in block.leftFace():
        #             if (newPoint == point and point.id not in [p.id for p in newPoints]):
        #                 newPoints[i] = point
        #                 isFound = True
        #                 break
        #         if (isFound):
        #             break
        #     if (not isFound):
        #         newPoints[i].setIndex()

        # for block in self.blocks:
        #     if (isEqualSubFace(newBlock.rightFace(), block.leftFace())):
        #         for point in newBlock.points:
        #             for point_r in newBlock.rightFace():

        #         break

        for i in range(len(newBlock.points)):
            if (newBlock.points[i].id == -1):
                newBlock.points[i].setIndex()
        # print([p.id for p in newPoints])

        # newBlock = self.createBlock(
        #     newName,
        #     newPoints,
        #     nx=targetBlock.nx,
        #     ny=targetBlock.ny,
        #     nz=targetBlock.nz,
        # )

        self.blocks.append(newBlock)

        return(newBlock)


    def getUniquePoints(self, blocks: list[Block], withEdge: bool=False) -> list[Point]:
        points = []
        for block in blocks:
            for point in block.points:
                if (id(point) not in [id(p) for p in points]):
                    points.append(point)
            if (withEdge):
                for edge in block.edges:
                    if (id(edge.midPoint) not in [id(p) for p in points]):
                        points.append(edge.midPoint)
        return(points)


    def addBoundary(self, face: Face) -> None:
        """
        Add a boundary to the face collection, can be used later for boundary
        condition definition.
        """
        check_type("face", face, Face)
        self.faces.append(face)


    def addTop(self, targetBlock: Block, name: str, points: list[Point], nz: int=1, gradz: float=1) -> Block:
        for block in self.blocks:
            if (block == targetBlock):
                newBlock = Block(
                    name, block.points[4:]+points, block.nx, block.ny, nz,
                    gradx=block.gradx, grady=block.grady, gradz=gradz
                )
                self.blocks.append(newBlock)
                return(newBlock)

    def addBottom(self, targetBlock: Block, name: str, points: list[Point], nz: int=1, gradz: float=1) -> Block:
        for block in self.blocks:
            if (block == targetBlock):
                newBlock = Block(
                    name, points+block.points[:4], block.nx, block.ny, nz,
                    gradx=block.gradx, grady=block.grady, gradz=gradz
                )
                self.blocks.append(newBlock)
                return(newBlock)

    def addFront(self, targetBlock: Block, name: str, points: list[Point], ny: int=1, grady: float=1) -> Block:
        for block in self.blocks:
            if (block == targetBlock):
                newBlock = Block(
                    name, [
                        points[0], points[1], block.points[1], block.points[0],
                        points[2], points[3], block.points[5], block.points[4]
                    ],
                    block.nx, ny, block.nz,
                    gradx=block.gradx, grady=grady, gradz=block.gradz
                )
                self.blocks.append(newBlock)
                return(newBlock)

    def addBack(self, targetBlock: Block, name: str, points: list[Point], ny: int=1, grady: float=1) -> Block:
        for block in self.blocks:
            if (block == targetBlock):
                newBlock = Block(
                    name, [
                        block.points[3], block.points[2], points[0], points[1],
                        block.points[7], block.points[6], points[2], points[3]
                    ],
                    block.nx, ny, block.nz,
                    gradx=block.gradx, grady=grady, gradz=block.gradz
                )
                self.blocks.append(newBlock)
                return(newBlock)

    def addLeft(self, targetBlock: Block, name: str, points: list[Point], nx: int=1, gradx: float=1) -> Block:
        for block in self.blocks:
            if (block == targetBlock):
                newBlock = Block(
                    name, [
                        points[0], block.points[0], block.points[3], points[1],
                        points[2], block.points[4], block.points[7], points[3]
                    ],
                    nx, block.ny, block.nz,
                    gradx=gradx, grady=block.grady, gradz=block.gradz
                )
                self.blocks.append(newBlock)
                return(newBlock)

    def addRight(self, targetBlock: Block, name: str, points: list[Point], nx: int=1, gradx: float=1) -> Block:
        for block in self.blocks:
            if (block == targetBlock):
                newBlock = Block(
                    name, [
                        block.points[1], points[0], points[1], block.points[2],
                        block.points[5], points[2], points[3], block.points[6]
                    ],
                    nx, block.ny, block.nz,
                    gradx=gradx, grady=block.grady, gradz=block.gradz
                )
                self.blocks.append(newBlock)
                return(newBlock)

    def getPointsForBlock(
            self,
            originalPoints: list[Point],
            existingPoints: dict[Point],
            fx=lambda x: x,
            fy=lambda y: y,
            fz=lambda z: z
        ) -> list[Point]:
        """
        Generate a list of Point objects that are unique in coordinates, i.e
        that are not overlapping.

        Parameters
        ----------
        originalPoints : list[Point]
            List of base points, e.g reference points for extrusion
        existingPoints : dict[Point]
            List of existing points, is filled by this function if necessary
        fx : lambda
            Displacement along the X-axis
        fy : lambda
            Displacement along the Y-axis
        fz : lambda
            Displacement along the Z-axis
        """
        newPointsForCurrBlock: list[Point] = []
        for point in originalPoints:
            newPointCoord = (fx(point.x), fy(point.y), fz(point.z))
            if (newPointCoord not in existingPoints.keys()):
                newPoint = Point(newPointCoord[0], newPointCoord[1], newPointCoord[2])
                existingPoints[newPointCoord] = newPoint
                newPointsForCurrBlock.append(newPoint)
            else:
                newPointsForCurrBlock.append(existingPoints[newPointCoord])
        return(newPointsForCurrBlock)

    def extrudeTop(self, targetBlocks: list[Block], name: str, dz: float, nz: int=1, gradz: float=1) -> tuple[Block] | Block:
        """
        Extrude a block using the top faces of the targetBlocks. The face of the
        new blocks is the same as the targetBlocks top faces. nx and ny are also
        preserved.
        """
        newBlocks = []
        newPoints: dict[Point] = {}
        for targetBlock in targetBlocks:
            newPointsForCurrBlock = self.getPointsForBlock(
                targetBlock.points[4:],
                newPoints,
                fz=lambda z: z + dz
            )

            newBlocks.append(self.addTop(
                targetBlock, name, newPointsForCurrBlock, nz=nz, gradz=gradz
            ))

        if (len(newBlocks) > 1):
            return tuple(newBlocks)
        return(newBlocks[0])

    def extrudeBottom(self, targetBlocks: list[Block], name: str, dz: float, nz: int=1, gradz: float=1) -> tuple[Block] | Block:
        """
        Extrude a block using the bottom faces of the targetBlocks. The face of the
        new blocks is the same as the targetBlocks bottom faces. nx and ny are also
        preserved.
        """
        newBlocks = []
        newPoints: dict[Point] = {}
        for targetBlock in targetBlocks:
            newPointsForCurrBlock = self.getPointsForBlock(
                targetBlock.points[:4],
                newPoints,
                fz=lambda z: z - dz
            )

            newBlocks.append(self.addBottom(
                targetBlock, name, newPointsForCurrBlock, nz=nz, gradz=gradz
            ))

        if (len(newBlocks) > 1):
            return tuple(newBlocks)
        return(newBlocks[0])

    def extrudeFront(self, targetBlocks: list[Block], name: str, dy: float, ny: int=1, grady: float=1) -> tuple[Block] | Block:
        """
        Extrude a block using the front faces of the targetBlocks. The face of the
        new blocks is the same as the targetBlocks front faces. nx and nz are also
        preserved.
        """
        newBlocks = []
        newPoints: dict[Point] = {}
        for targetBlock in targetBlocks:
            newPointsForCurrBlock = self.getPointsForBlock(
                [targetBlock.points[0], targetBlock.points[1], targetBlock.points[4], targetBlock.points[5]],
                newPoints,
                fy=lambda y: y - dy
            )

            newBlocks.append(self.addFront(
                targetBlock, name, newPointsForCurrBlock, ny=ny, grady=grady
            ))

        if (len(newBlocks) > 1):
            return tuple(newBlocks)
        return(newBlocks[0])

    def extrudeBack(self, targetBlocks: list[Block], name: str, dy: float, ny: int=1, grady: float=1) -> tuple[Block] | Block:
        """
        Extrude a block using the back faces of the targetBlocks. The face of the
        new blocks is the same as the targetBlocks back faces. nx and nz are also
        preserved.
        """
        newBlocks = []
        newPoints: dict[Point] = {}
        for targetBlock in targetBlocks:
            newPointsForCurrBlock = self.getPointsForBlock(
                [targetBlock.points[2], targetBlock.points[3], targetBlock.points[6], targetBlock.points[7]],
                newPoints,
                fy=lambda y: y + dy
            )

            newBlocks.append(self.addBack(
                targetBlock, name, newPointsForCurrBlock, ny=ny, grady=grady
            ))

        if (len(newBlocks) > 1):
            return tuple(newBlocks)
        return(newBlocks[0])

    def extrudeLeft(self, targetBlocks: list[Block], name: str, dx: float, nx: int=1, gradx: float=1) -> tuple[Block] | Block:
        """
        Extrude a block using the back left of the targetBlocks. The face of the
        new blocks is the same as the targetBlocks back left. ny and nz are also
        preserved.
        """
        newBlocks = []
        newPoints: dict[Point] = {}
        for targetBlock in targetBlocks:
            newPointsForCurrBlock = self.getPointsForBlock(
                [targetBlock.points[0], targetBlock.points[3], targetBlock.points[4], targetBlock.points[7]],
                newPoints,
                fx=lambda x: x - dx
            )

            newBlocks.append(self.addLeft(
                targetBlock, name, newPointsForCurrBlock, nx=nx, gradx=gradx
            ))

        if (len(newBlocks) > 1):
            return tuple(newBlocks)
        return(newBlocks[0])

    def extrudeRight(self, targetBlocks: list[Block], name: str, dx: float, nx: int=1, gradx: float=1) -> tuple[Block] | Block:
        """
        Extrude a block using the back right of the targetBlocks. The face of the
        new blocks is the same as the targetBlocks back right. ny and nz are also
        preserved.
        """
        newBlocks = []
        newPoints: dict[Point] = {}
        for targetBlock in targetBlocks:
            newPointsForCurrBlock = self.getPointsForBlock(
                [targetBlock.points[1], targetBlock.points[2], targetBlock.points[5], targetBlock.points[6]],
                newPoints,
                fx=lambda x: x + dx
            )

            newBlocks.append(self.addRight(
                targetBlock, name, newPointsForCurrBlock, nx=nx, gradx=gradx
            ))

        if (len(newBlocks) > 1):
            return tuple(newBlocks)
        return(newBlocks[0])

    def extrudeNormal(
            self,
            targetBlock: Block,
            facename: str,
            name: str,
            length: float,
            n: int=1,
            grad: float=1
        ) -> Block:
        """
        Extrude a block from the normal face of the target block.

        Equivalence between `facename` and new block generation function:
        - `top`: `addTop(...)`
        - `bottom`: `addBottom(...)`
        - `left`: `addLeft(...)`
        - `right`: `addRight(...)`
        - `front`: `addFront(...)`
        - `back`: `addBack(...)`

        Parameters
        ----------
        targetBlock : Block
            Target block to extrude.
        facename : str
            Name of the face to extrude from the `targetBlock`. Options (`top`,
            `bottom`, `left`, `right`, `front`, `back`).
        name : str
            Name of the extruded block.
        length : float
            Length of the extrusion.
        n : int
            Number of cells along the extrusion (default `1`).
        grad : float
            Grading along the extrusion (default `1`).
        """
        check_type("targetBlock", targetBlock, Block)
        check_type("facename", facename, str)
        check_value("facename", facename, _FACE_NAME_TYPES)
        check_type("name", name, str)
        check_type("length", length, (float, int))
        check_type("n", n, int)
        check_type("grad", grad, (float, int))

        normalVector: Vector = targetBlock.getFaceNormal(facename)
        normalVector.normalize(length)

        facePoints = targetBlock.getFace(facename)
        newPoints = []
        for point in facePoints:
            newPos = Vector(point.x, point.y, point.z) + normalVector
            newPoints.append(Point(newPos.x, newPos.y, newPos.z))

        if (facename == 'top'):
            return(self.addTop(targetBlock=targetBlock, name=name, points=newPoints, nz=n, gradz=grad))
        if (facename == 'bottom'):
            newPoints[1], newPoints[3] = newPoints[3], newPoints[1]
            return(self.addBottom(targetBlock=targetBlock, name=name, points=newPoints, nz=n, gradz=grad))
        if (facename == 'left'):
            newPoints[1], newPoints[2], newPoints[3] = newPoints[3], newPoints[1], newPoints[2]
            return(self.addLeft(targetBlock=targetBlock, name=name, points=newPoints, nx=n, gradx=grad))
        if (facename == 'right'):
            newPoints[2], newPoints[3] = newPoints[3], newPoints[2]
            return(self.addRight(targetBlock=targetBlock, name=name, points=newPoints, nx=n, gradx=grad))
        if (facename == 'front'):
            newPoints[2], newPoints[3] = newPoints[3], newPoints[2]
            return(self.addFront(targetBlock=targetBlock, name=name, points=newPoints, ny=n, grady=grad))
        if (facename == 'back'):
            newPoints[0], newPoints[1], newPoints[3] = newPoints[3], newPoints[0], newPoints[1]
            return(self.addBack(targetBlock=targetBlock, name=name, points=newPoints, ny=n, grady=grad))


    def extrudeNormalArc(
            self,
            targetBlock: Block,
            facename: str,
            name: str,
            angleSpan: float,
            arcCenter: Vector,
            nt: int=1,
            grad: float=1,
            rotationAxis: str='z'
        ) -> Block:
        """
        Extrude arc block, curved along the `rotationAxis`, from the normal face
        of the target block. Fully rely on `self.extrudeNormal(...)`

        Parameters
        ----------
        targetBlock : Block
            Target block to extrude.
        facename : str
            Name of the face to extrude from the `targetBlock`. Options (`top`,
            `bottom`, `left`, `right`, `front`, back`).
        name : str
            Name of the extruded block.
        angleSpan : float
            Arc angle in degrees. Create an arc in the anti-clockwise
            orientation if `angleSpan` is positive.
        arcCenter : Vector
            Origin position of the arc.
        nt : int
            Number of cells along the extrusion (default `1`).
        grad : float
            Grading along the extrusion (default `1`).
        rotationAxis : str
            Arc rotation axis around a cylinder. Options (`x`, `y`, `z`).
        """
        check_type("targetBlock", targetBlock, Block)
        check_type("facename", facename, str)
        check_value("facename", facename, _FACE_NAME_TYPES)
        check_type("name", name, str)
        check_type("angleSpan", angleSpan, (float, int))
        check_type("arcCenter", arcCenter, Vector)
        check_type("nt", nt, int)
        check_type("rotationAxis", rotationAxis, str)
        check_value("rotationAxis", rotationAxis, {'x', 'y', 'z'})

        newBlock = self.extrudeNormal(targetBlock, facename, name, length=0, n=nt, grad=grad)

        if (facename == "right"):
            points = newBlock.rightFace()
            arcPoints = [(0, 1), (2, 3), (4, 5), (6, 7)]
        elif (facename == "left"):
            points = newBlock.leftFace()
            arcPoints = [(0, 1), (2, 3), (4, 5), (6, 7)]
        elif (facename == "front"):
            points = newBlock.frontFace()
            arcPoints = [(0, 3), (2, 1), (4, 7), (6, 5)]
        elif (facename == "back"):
            points = newBlock.backFace()
            arcPoints = [(0, 3), (2, 1), (4, 7), (6, 5)]
        elif (facename == "top"):
            points = newBlock.topFace()
            arcPoints = [(0, 4), (1, 5), (2, 6), (3, 7)]
        elif (facename == "bottom"):
            points = newBlock.bottomFace()
            arcPoints = [(0, 4), (1, 5), (2, 6), (3, 7)]

        #  Add arcs
        centerX = arcCenter.x if rotationAxis != 'x' else None
        centerY = arcCenter.y if rotationAxis != 'y' else None
        centerZ = arcCenter.z if rotationAxis != 'z' else None
        for pointA, pointB in arcPoints:
            newBlock.addEdgeArc(pointA, pointB, centerX, centerY, centerZ, isOrigin=True)

        # Compute displacement
        dx = arcCenter.x if rotationAxis != 'x' else 0
        dy = arcCenter.y if rotationAxis != 'y' else 0
        dz = arcCenter.z if rotationAxis != 'z' else 0

        # Update position
        for point in points:
            # Move to origin
            point.translate(-dx, -dy, -dz)

            # Rotate
            if (rotationAxis == 'x'):
                point.rotateX(angleSpan * np.pi/180)
            elif (rotationAxis == 'y'):
                point.rotateY(angleSpan * np.pi/180)
            elif (rotationAxis == 'z'):
                point.rotateZ(angleSpan * np.pi/180)

            # Back to initial location
            point.translate(dx, dy, dz)

        return(newBlock)


    def addBySymmetryX(
            self,
            name: str,
            faceToClone: list[Point],
            targetFaceToAttach: list[Point],
            nx: int=1
        ) -> None:
        self.createBlock(name, [
            Point(-targetFaceToAttach[0].x, targetFaceToAttach[0].y, targetFaceToAttach[0].z),
            faceToClone[0],
            faceToClone[1],
            Point(-targetFaceToAttach[1].x, targetFaceToAttach[1].y, targetFaceToAttach[1].z),
            Point(-targetFaceToAttach[2].x, targetFaceToAttach[2].y, targetFaceToAttach[2].z),
            faceToClone[2],
            faceToClone[3],
            Point(-targetFaceToAttach[3].x, targetFaceToAttach[3].y, targetFaceToAttach[3].z),
        ], nx=nx)


    def add_pipes_from_json(
            self,
            filename: str,
            originPosition: Vector | Block=Vector(0, 0, 0),
            firstBlockName: str=None,
            elbowRadius: float=0,
            isAddBoundaryConditions: bool=False
        ):
        """
        Create a piping circuit from a JSON file (`filename`) generated by
        `foamForNuclear.mesh.PipingGUI()`. This piping GUI can be called via:

            pipingGui = foamForNuclear.mesh.PipingGUI()
            pipingGui.run()

        The pipes are created using `addPipe1DFromDirection` for each block and
        `addPipe1DFrom2Points` for the last block to close the loop if any. One
        can specify the first starting block using the `firstBlockName`
        parameter.

        Note:
        - Currently tested for a single loop.
        - The `elbowRadius` is the same between every pipe section.

        Parameters
        ----------
        filename : str
            Name of the JSON circuitry file.
        originPosition : Vector | Block
            Origin position of the first block.
        firstBlockName : str
            Name of the first block to place from the circuitry file.
        elbowRadius : float
            Elbow radius applied to all the piping links (default 0).
        isAddBoundaryConditions : bool
            (default False)
        """
        with open(filename, 'r', encoding='utf-8') as f:
            data = json.load(f)

        # Extract pipe links
        linksTargetFromSource = {}
        for edge in data.get('edges', []):
            source = edge.get('source')
            target = edge.get('target')
            sourceHandle = edge.get('sourceHandle')
            # targetHandle = edge.get('targetHandle')
            if (sourceHandle == 'r'):
                linksTargetFromSource[target] = source
            else:
                linksTargetFromSource[source] = target

        # Place first block
        blocks = {}
        for node in data.get('nodes', []):
            node_data = node.get('data')
            name = node_data.get('label')
            if (name != firstBlockName):
                continue

            nid = node.get('id')
            rotation = -node_data.get('rotation', 0) * np.pi/180
            length = node_data.get('length_m')
            width = node_data.get('width_m')
            nCells = int(node_data.get('cells', 1))

            newBlock = self.addPipe1DFromDirection(
                name=f"{name}_{nid}",
                originPosition=originPosition,
                direction=Vector(np.cos(rotation), 0, np.sin(rotation)),
                length=length,
                equivalentHydraulicDiameter=width,
                n=nCells,
                elbowRadius=elbowRadius,
                isAddBoundaryConditions=isAddBoundaryConditions
            )

            blocks[nid] = newBlock
            break

        # Extract blocks and continue the loop
        for node in data.get('nodes', []):
            nid = node.get('id')
            if (nid in blocks.keys()):
                continue

            node_data = node.get('data')
            name = node_data.get('label')
            rotation = -node_data.get('rotation', 0) * np.pi/180
            isNextBlockPlaced = False

            if (len(blocks.keys()) > 0):
                originBlock = linksTargetFromSource[nid]
                originPosition = blocks[originBlock]

                nextBlock = [key for key, nid_ in linksTargetFromSource.items() if nid == nid_][0]
                isNextBlockPlaced = nextBlock in blocks.keys()

            length = node_data.get('length_m')
            width = node_data.get('width_m')
            nCells = int(node_data.get('cells', 1))

            if (isNextBlockPlaced):
                newBlock = self.addPipe1DFrom2Points(
                    name=f"{name}_{nid}",
                    originPosition=originPosition,
                    finalPosition=blocks[nextBlock],
                    equivalentHydraulicDiameter=width,
                    n=nCells,
                    elbowRadius=elbowRadius,
                    isAddBoundaryConditions=isAddBoundaryConditions
                )
            else:
                newBlock = self.addPipe1DFromDirection(
                    name=f"{name}_{nid}",
                    originPosition=originPosition,
                    direction=Vector(np.cos(rotation), 0, np.sin(rotation)),
                    length=length,
                    equivalentHydraulicDiameter=width,
                    n=nCells,
                    elbowRadius=elbowRadius,
                    isAddBoundaryConditions=isAddBoundaryConditions
                )

            blocks[nid] = newBlock

        return(blocks)


    def printBlocks(self) -> str:
        text = "blocks\n"
        text += "(\n"
        for block in self.blocks:
            if (self.isReducedCells):
                block.nx = 3
                block.ny = 3
                block.nz = 3
            if (block.isPrint):
                text += f"{tab}{block}\n"
        text += ");\n"
        return(text)

    def printPoints(self) -> str:
        minId = 1e32
        allPoints = []
        # Find min point id
        for block in self.blocks:
            for point in block.points:
                if (point.id not in [p.id for p in allPoints]):
                    allPoints.append(point)
                if (point.id < minId):
                    minId = point.id

        # Reset point id if necessary
        for point in allPoints:
            point.id = point.id - minId

        pointsPlaced = []
        for block in self.blocks:
            for point in block.points:
                if (point.id not in [point.id for point in pointsPlaced]):
                    pointsPlaced.append(point)

        text = "vertices\n"
        text += "(\n"
        for point in sorted(pointsPlaced, key=lambda p: p.id):
            text += f"{tab}{point}\n"
        text += ");\n"
        return(text)

    def printEdges(self) -> str:
        text = "edges\n"
        text += "(\n"
        for block in self.blocks:
            for edge in block.edges:
                text += f"{tab}{edge}\n"

            for edgeProjection in block.edgeProjection:
                text += f"{tab}projectCurve {edgeProjection['verticeIdx1']} {edgeProjection['verticeIdx2']} ({edgeProjection['geometryNames']})\n"
        text += ");\n"
        return(text)

    def printFaces(self) -> str:
        text = "faces\n(\n"
        for i, block in enumerate(self.blocks):
            for faceProjection in block.faceProjection:
                text += f"{tab}project ({i} {faceProjection['faceId']}) {faceProjection['geometryName']}\n"
        text += ");\n"
        return(text)

    def printBoundaries(self) -> str:
        text = "boundary\n"
        text += "(\n"
        for face in self.faces:
            if (face.isPrint and not face.is_empty):
                text += f"{face}\n"
        text += ");\n"
        return(text)

    def printMergePatchPairs(self) -> str:
        text = "mergePatchPairs\n(\n"
        for facei, facej in self.mergePatchPairs:
            if (isinstance(facei, str) and isinstance(facej, str)):
                text += f"{tab}({facei} {facej})\n"
            elif (isinstance(facei, Face) and isinstance(facej, Face)):
                text += f"{tab}({facei.name} {facej.name})\n"
            else:
                raise ValueError("mergePatchPairs must contains str or Face")
        text += ");\n"
        return(text)

    def printGeometries(self) -> str:
        """
        Print `geometry` section in the blockMeshDict. This section usually
        groups the geometries shape used for projection of edge and faces.
        """
        text = "geometry\n{\n"

        # Spheres
        for sphere in self.spheres:
            text += f"{tab}{sphere['name']}\n"
            text += f"{tab}"+"{\n"
            text += f"{2*tab}type   sphere;\n"
            text += f"{2*tab}origin ({' '.join([str(e) for e in list(sphere['origin'])])});\n"
            text += f"{2*tab}radius {sphere['radius']};\n"
            text += f"{tab}"+"}\n"

        # Cylinders
        for cylinder in self.cylinders:
            text += f"{tab}{cylinder['name']}\n"
            text += f"{tab}"+"{\n"
            text += f"{2*tab}type   cylinder;\n"
            text += f"{2*tab}point1 ({' '.join([str(e) for e in list(cylinder['point1'])])});\n"
            text += f"{2*tab}point2 ({' '.join([str(e) for e in list(cylinder['point2'])])});\n"
            text += f"{2*tab}radius {cylinder['radius']};\n"
            text += f"{tab}"+"}\n"

        # Cones
        for cone in self.cones:
            text += f"{tab}{cone['name']}\n"
            text += f"{tab}"+"{\n"
            text += f"{2*tab}type       cone;\n"
            text += f"{2*tab}point1     ({' '.join([str(e) for e in list(cone['point1'])])});\n"
            text += f"{2*tab}point2     ({' '.join([str(e) for e in list(cone['point2'])])});\n"
            text += f"{2*tab}radius1    {cone['radius1']};\n"
            text += f"{2*tab}radius2    {cone['radius2']};\n"
            text += f"{tab}"+"}\n"

        text += "}\n"
        return(text)

    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        """
        Export the mesh into OpenFOAM's blockMeshDict.
        """
        Mesh.export_to_openfoam(self)

        text = ""

        text += addParameter("scale", self.scale, isAddExtraLine=True)

        text += self.printGeometries()+"\n"

        text += self.printPoints()+"\n"

        text += self.printBlocks()+"\n"

        text += self.printEdges()+"\n"

        text += self.printFaces()+"\n"

        text += self.printBoundaries()+"\n"

        text += self.printMergePatchPairs()+"\n"

        return(text)


class BlockMeshWedge(BlockMesh):
    """
    Special blockMesh builder for wedge geometry.
    The wedge faces are called "front" and "back".

    Parameters
    ----------
    appertureAngle : float
        Wedge angle in degree

    Attributes
    ----------
    faceFrontWedge : Face
        Front face of the wedge
    faceBackWedge : Face
        Back face of the wedge
    """
    def __init__(
            self,
            appertureAngle: float,
            region: str="",
            scale: float=1
        ):
        super().__init__(region, scale)

        self.appertureAngle = appertureAngle

        self.faceFrontWedge = Face("front", boundaryType='wedge')
        self.faceBackWedge  = Face("back", boundaryType='wedge')

        self.addBoundary(self.faceFrontWedge)
        self.addBoundary(self.faceBackWedge)


    @property
    def appertureAngle(self):
        return self._appertureAngle

    @appertureAngle.setter
    def appertureAngle(self, appertureAngle) -> None:
        check_type("appertureAngle", appertureAngle, (float, int))
        self._appertureAngle = appertureAngle * np.pi/180


    def createWedge(
            self,
            name: str,
            innerRadius: float,
            outerRadius: float,
            lowZ: float,
            highZ: float,
            nr: int, nz: int
        ) -> Block:
        """
        Redefined `BlockMesh.createWedge` by adding wedge patches
        """
        block = super().createWedge(
            name, innerRadius, outerRadius, lowZ, highZ, self.appertureAngle*180/np.pi, nr, nz
        )

        self.faceFrontWedge.addSubFace(block.frontFace())
        self.faceBackWedge.addSubFace(block.backFace())

        return(block)


    def createWedgeConical(
            self,
            name: str,
            innerRadiusBottom: float,
            outerRadiusBottom: float,
            innerRadiusTop: float,
            outerRadiusTop: float,
            lowZ: float,
            highZ: float,
            nr: int,
            nz: int
        ) -> Block:
        """
        Redefined `BlockMesh.createWedgeConical` by adding wedge patches
        """
        block = super().createWedgeConical(
            name,
            innerRadiusBottom, outerRadiusBottom,
            innerRadiusTop, outerRadiusTop,
            lowZ, highZ,
            self.appertureAngle*180/np.pi,
            nr, nz
        )

        self.faceFrontWedge.addSubFace(block.frontFace())
        self.faceBackWedge.addSubFace(block.backFace())

        return(block)


    def extrudeTop(self, targetBlocks: list[Block], name: str, dz: float, nz: int=1) -> tuple[Block]:
        """
        Redefined `BlockMesh.extrudeTop` by adding wedge patches
        """
        blocks = super().extrudeTop(targetBlocks, name, dz, nz)

        if (isinstance(blocks, Block)):
            self.faceFrontWedge.addSubFace(blocks.frontFace())
            self.faceBackWedge.addSubFace(blocks.backFace())
        else:
            for block in blocks:
                self.faceFrontWedge.addSubFace(block.frontFace())
                self.faceBackWedge.addSubFace(block.backFace())

        return blocks


    def extrudeBottom(self, targetBlocks: list[Block], name: str, dz: float, nz: int=1) -> tuple[Block]:
        """
        Redefined `BlockMesh.extrudeBottom` by adding wedge patches
        """
        blocks = super().extrudeBottom(targetBlocks, name, dz, nz)

        if (isinstance(blocks, Block)):
            self.faceFrontWedge.addSubFace(blocks.frontFace())
            self.faceBackWedge.addSubFace(blocks.backFace())
        else:
            for block in blocks:
                self.faceFrontWedge.addSubFace(block.frontFace())
                self.faceBackWedge.addSubFace(block.backFace())

        return blocks


    def extrudeLeft(self, targetBlocks: list[Block], name: str, dr: float, nr: int=1) -> tuple[Block]:
        """
        Redefined `BlockMesh.extrudeLeft` by adding wedge patches and keeping
        the wedge angle hwile extruding
        """
        newBlocks = []
        newPoints: dict[Point] = {}
        for targetBlock in targetBlocks:
            newPointsForCurrBlock = self.getPointsForBlock(
                [targetBlock.points[0], targetBlock.points[3], targetBlock.points[4], targetBlock.points[7]],
                newPoints,
                fx=lambda x: x - dr*np.cos(self.appertureAngle/2),
                fy=lambda y: y - dr*np.sin(self.appertureAngle/2) if y > 0 else y + dr*np.sin(self.appertureAngle/2),
            )

            newBlock = self.addLeft(
                targetBlock, name, newPointsForCurrBlock, nx=nr
            )

            newBlocks.append(newBlock)
            self.faceFrontWedge.addSubFace(newBlock.frontFace())
            self.faceBackWedge.addSubFace(newBlock.backFace())

        if (len(newBlocks) > 1):
            return tuple(newBlocks)
        return(newBlocks[0])


    def extrudeRight(self, targetBlocks: list[Block], name: str, dr: float, nr: int=1) -> tuple[Block]:
        """
        Redefined `BlockMesh.extrudeRight` by adding wedge patches and keeping
        the wedge angle hwile extruding
        """
        newBlocks = []
        newPoints: dict[Point] = {}
        for targetBlock in targetBlocks:
            newPointsForCurrBlock = self.getPointsForBlock(
                [targetBlock.points[1], targetBlock.points[2], targetBlock.points[5], targetBlock.points[6]],
                newPoints,
                fx=lambda x: x + dr*np.cos(self.appertureAngle/2),
                fy=lambda y: y - dr*np.sin(self.appertureAngle/2) if y < 0 else y + dr*np.sin(self.appertureAngle/2),
            )

            newBlock = self.addRight(
                targetBlock, name, newPointsForCurrBlock, nx=nr
            )

            newBlocks.append(newBlock)
            self.faceFrontWedge.addSubFace(newBlock.frontFace())
            self.faceBackWedge.addSubFace(newBlock.backFace())

        if (len(newBlocks) > 1):
            return tuple(newBlocks)
        return(newBlocks[0])


    def addRightFaceEdgePolyLine(self, block: Block, rzCoords: list[tuple]) -> None:
        """
        Deforme the right face of the mesh following the rzCoords list in (r, z).

        Parameters
        ----------
        block : Block
            Block to attach the face deformation
        rzCoords : list[tuple]
            List of (r, z) tuple
        """
        pointsFront, pointsBack = [], []
        for r, z in rzCoords:
            pointsFront.append(Vector(
                r*np.cos(-self.appertureAngle/2),
                r*np.sin(-self.appertureAngle/2),
                z
            ))
            pointsBack.append(Vector(
                r*np.cos(self.appertureAngle/2),
                r*np.sin(self.appertureAngle/2),
                z
            ))

        block.addEdgePolyLine(
            pointIdx1=1,
            pointIdx2=5,
            points=pointsFront
        )
        block.addEdgePolyLine(
            pointIdx1=2,
            pointIdx2=6,
            points=pointsBack
        )
