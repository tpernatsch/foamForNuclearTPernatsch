from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import *
from foamForNuclear.openfoamFile import OpenFOAMFile


_SET_TYPE_TYPES = {
    'pointSet', 'pointZoneSet', 'faceSet', 'faceZoneSet', 'cellSet',
    'cellZoneSet'
}
_ACTION_TYPE_TYPES = {
    'add', 'subtract', 'new', 'subset', 'invert', 'clear', 'remove', 'list',
    'ignore', 'delete'
}
_SOURCE_TYPE_TYPES = {
    # pointSet
    'boxToPoint', # : select points inside bound box(es)
    'cellToPoint', # : select points from cell set(s)
    'clipPlaneToPoint', # : select points above a plane
    'cylinderToPoint', # : select points in a cylinder or annulus
    'faceToPoint', # : select points in face set(s)
    'labelToPoint', # : select points from a list of point indices
    'nearestToPoint', # : select points closest to a set of points
    'pointToPoint', # : select points in point set(s)
    'searchableSurfaceToPoint', # : select points enclosed by a searchableSurface
    'setToPointZone', # : convert a pointSet to a pointZone
    'sphereToPoint', # : select points within a bounding sphere
    'surfaceToPoint', # : select points based on relation (distance, inside/outside) to a surface
    'zoneToPoint', # : convert a pointZone to a pointSet
    # faceSet
    'boxToFace', # : select faces inside bound box(es)
    'cellToFace', # : select faces from cell set(s)
    'clipPlaneToFace', # : select faces above a plane
    'cylinderAnnulusToFace', # : select faces in a cylinder annulus
    'cylinderToFace', # : select faces in a cylinder
    'faceToFace', # : select faces in face set(s)
    'holeToFace', # : select a set of faces that closes a hole
    'labelToFace', # : select faces from a list of face indices
    'normalToFace', # : select faces with normal aligned to a specified direction
    'patchToFace', # : select faces associated with given patch(es)
    'pointToFace', # : select faces connected to any points in point set(s)
    'regionToFace', # : select faces in a mesh region
    'searchableSurfaceToFace', # : select faces enclosed by a searchableSurface
    'sphereToFace', # : select faces within a bounding sphere
    'zoneToFace', # : convert a faceZone to a faceSet
    # faceZoneSet
    'cellToFaceZone', # : select faces with only 1 neighbour in cell set(s)
    'faceZoneToFaceZone', # : select faces in a faceZone
    'planeToFaceZone', # : select faces based on the adjacent cell centres spanning a given plane
    'searchableSurfaceToFaceZone', # : select faces whose cell-cell centre vector intersects a given searchableSurface
    'setAndNormalToFaceZone', # : select faces from a faceSet where the face normal is aligned to a specified direction
    'setToFaceZone', # : convert a faceSet to a faceZone
    'setsToFaceZone', # : convert a list of faceSets to a faceZone
    # cellSet
    'boxToCell', # : select cells inside bound box(es)
    'cellToCell', # : select cells from in cell set(s)
    'clipPlaneToCell', # : select cells with centres above a given plane
    'cylinderAnnulusToCell', # : select cells in a cylinder annulus
    'cylinderToCell', # : select cells in a cylinder
    'faceToCell', # : select cells with a face in the given face set(s)
    'faceZoneToCell', # : select cells with a face in the given face zone(s)
    'fieldToCell', # : select cells based on field values
    'haloToCell', # : select cells connected to the outside a cellSet
    'labelToCell', # : select cells from a list of cell indices
    'nbrToCell', # : select cells based on each cell’s number of neighbour cells
    'nearestToCell', # : select cells with centres nearest to a set of points
    'patchToCell', # : select cells adjacent to a lost of patches
    'pointToCell', # : select cells with with any points in the given point set(s)
    'regionToCell', # : select cells in a mesh region
    'rotatedBoxToCell', # : select sells in a rotated box
    'searchableSurfaceToCell', # : select cells enclosed by a searchableSurface
    'shapeToCell', # : select cells based on their shape
    'sphereToCell', # : select cells within a bounding sphere
    'surfaceToCell', # : select cells based on relation (distance, inside/outside) to a surface
    'targetVolumeToCell', # : selects cells based on a target volume obtained by sweeping a surface through the mesh
    'zoneToCell', # : converts a cellZone to a cellSet
    # cellZoneSet
    'setToCellZone', # : converts a cellSet to a cellZone
}


class TopoSetAction(OpenFOAMDict):
    """
    Action used exclusively in topoSetDict.
    """
    def __init__(
            self,
            name: str,
            setType: str,
            actionType: str,
            sourceType: str,
            extraParameters: dict={}
        ):
        super().__init__()
        self.actionName = name
        self.setType = setType
        self.actionType = actionType
        self.sourceType = sourceType
        self.extraParameters = extraParameters


    def __repr__(self):
        for key, param in self.extraParameters.items():
            self.__setitem__(key, param)
        return(super().__repr__(depth=1))

    @property
    def actionName(self):
        return self._actionName

    @actionName.setter
    def actionName(self, actionName):
        check_type('actionName', actionName, str)
        self._actionName = actionName
        self.__setitem__("name", actionName)

    @property
    def setType(self):
        return self._setType

    @setType.setter
    def setType(self, setType):
        check_type('setType', setType, str)
        check_value('setType', setType, _SET_TYPE_TYPES)
        self._setType = setType
        self.__setitem__("type", setType)

    @property
    def actionType(self):
        return self._actionType

    @actionType.setter
    def actionType(self, actionType):
        check_type('actionType', actionType, str)
        check_value('actionType', actionType, _ACTION_TYPE_TYPES)
        self._actionType = actionType
        self.__setitem__("action", actionType)

    @property
    def sourceType(self):
        return self._sourceType

    @sourceType.setter
    def sourceType(self, sourceType):
        check_type('sourceType', sourceType, str)
        check_value('sourceType', sourceType, _SOURCE_TYPE_TYPES)
        self._sourceType = sourceType
        self.__setitem__("source", sourceType)



class TopoSetDict(OpenFOAMFile):
    """
    Operates on cellSets/faceSets/pointSets through a dictionary, normally
    `system/topoSetDict`.
    """
    def __init__(self, region: str=""):
        super().__init__("topoSetDict", folder="system", region=region)

        self.actions: list[TopoSetAction] = []


    def append(self, item):
        """
        Parameters
        ----------
        item : TopoSetAction
        """
        check_type("item", item, TopoSetAction)
        self.actions.append(item)


    def add_action(
            self,
            name: str,
            setType: str,
            actionType: str,
            sourceType: str,
            extraParameters: dict={}
        ):
        self.append(TopoSetAction(
            name=name,
            setType=setType,
            actionType=actionType,
            sourceType=sourceType,
            extraParameters=extraParameters
        ))

    @property
    def is_empty(self):
        return(len(self.actions) == 0)


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""

        text += "actions\n(\n"

        for action in self.actions:
            text += f"{action!r}"

        text += ");\n"

        return(text)
