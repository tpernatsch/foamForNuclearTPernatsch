from matplotlib.patches import Patch
from foamForNuclear.boundaryConditions.zeroGradient import ZeroGradient
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import *
from foamForNuclear.openfoamFile import OpenFOAMFile

_FACE_TYPE_TYPES = {"faceZone", "searchableSurface"}
_SAMPLE_MODE_TYPES = {"nearestCell", "nearestPatchFace", "nearestPatchFaceAMI", "nearestPatchPoint", "nearestFace", "nearestOnlyCell"}


class BaffleDict:
    """
    Baffle dict used exclusively in createBafflesDict.
    """
    def __init__(
            self,
            name: str,
            faceType: str,
            zoneName: str,
            masterType: str="wall",
            slaveType: str="wall",
            neighbourPatch: str=None,
            sampleMode: str=None,
            sampleRegion: str=None,
            offsetMode: str=None,
            transform: str=None,
            offset: Vector=None,
            patchFields: OpenFOAMDict=None
        ):
        self.name = name
        self.faceType = faceType
        self.zoneName = zoneName
        self.masterType = masterType
        self.slaveType = slaveType
        self.neighbourPatch = neighbourPatch
        self.sampleMode = sampleMode
        self.sampleRegion = sampleRegion
        self.offsetMode = offsetMode
        self.transform = transform
        self.offset = offset
        self.patchFields = patchFields

        self.add_field_boundary_condition('".*"', ZeroGradient())


    def __repr__(self):
        text = ""
        text += f"{tab}{self.name}\n"
        text += tab + "{\n"
        text += addParameter("type", self.faceType, indent=2)
        text += addParameter("zoneName", self.zoneName, indent=2)

        text += f"{2*tab}patches\n"
        text += 2*tab + "{\n"

        text += f"{3*tab}master\n"
        text += 3*tab + "{\n"
        text += addParameter("name", f"{self.name}_master", indent=4)
        text += addParameter("type", self.masterType, indent=4)
        text += addParameter("neighbourPatch", f"{self.name}_slave", indent=4)

        if (self.masterType == "mappedWall"):
            text += addParameter("sampleMode", self.sampleMode, indent=4)
            text += addParameter("samplePatch", f"{self.name}_slave", indent=4)
            text += addParameter("sampleRegion", self.sampleRegion, indent=4)
            text += addParameter("offsetMode", self.offsetMode, indent=4)
            text += addParameter("transform", self.transform, indent=4)
            text += addParameter("offset", self.offset, indent=4)

        text += f"{4*tab}patchFields{self.patchFields.__repr__(depth=4)}"

        text += 3*tab + "}\n"

        text += f"{3*tab}slave\n"
        text += 3*tab + "{\n"
        text += addParameter("name", f"{self.name}_slave", indent=4)
        text += addParameter("type", self.slaveType, indent=4)
        text += addParameter("neighbourPatch", f"{self.name}_master", indent=4)

        if (self.masterType == "mappedWall"):
            text += addParameter("sampleMode", self.sampleMode, indent=4)
            text += addParameter("samplePatch", f"{self.name}_master", indent=4)
            text += addParameter("sampleRegion", self.sampleRegion, indent=4)
            text += addParameter("offsetMode", self.offsetMode, indent=4)
            text += addParameter("transform", self.transform, indent=4)
            text += addParameter("offset", self.offset, indent=4)

        text += f"{4*tab}patchFields\n"
        text += 4*tab + "{\n"
        text += 5*tab + "${...master.patchFields}\n"
        text += 4*tab + "}\n"

        text += 3*tab + "}\n"

        text += 2*tab + "}\n"

        text += tab + "}\n"
        return(text)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        check_type('name', name, str)
        self._name = name

    @property
    def faceType(self):
        return self._faceType

    @faceType.setter
    def faceType(self, faceType):
        check_type('faceType', faceType, str)
        check_value('faceType', faceType, _FACE_TYPE_TYPES)
        self._faceType = faceType

    @property
    def zoneName(self):
        return self._zoneName

    @zoneName.setter
    def zoneName(self, zoneName):
        check_type('zoneName', zoneName, str)
        self._zoneName = zoneName

    @property
    def masterType(self):
        return self._masterType

    @masterType.setter
    def masterType(self, masterType):
        check_type('masterType', masterType, str)
        self._masterType = masterType

    @property
    def slaveType(self):
        return self._slaveType

    @slaveType.setter
    def slaveType(self, slaveType):
        check_type('slaveType', slaveType, str)
        self._slaveType = slaveType

    @property
    def sampleMode(self):
        return self._sampleMode

    @sampleMode.setter
    def sampleMode(self, sampleMode):
        check_type('sampleMode', sampleMode, str, none_ok=True)
        if (sampleMode is not None):
            check_value('sampleMode', sampleMode, _SAMPLE_MODE_TYPES)
        self._sampleMode = sampleMode

    @property
    def patchFields(self):
        return self._patchFields

    @patchFields.setter
    def patchFields(self, patchFields):
        check_type('patchFields', patchFields, OpenFOAMDict, none_ok=True)
        if (patchFields is not None):
            self._patchFields = patchFields
        else:
            self._patchFields = OpenFOAMDict({})

    def add_field_boundary_condition(
            self,
            fieldName: str,
            boundaryCondition: Patch
        ):
        check_type("fieldName", fieldName, str)
        # check_type("boundaryCondition", boundaryCondition, Patch)

        self.patchFields[fieldName] = boundaryCondition


class CreateBafflesDict(OpenFOAMFile):
    """
    Makes internal faces into boundary faces. Does not duplicate points.

    Parameters
    ----------
    region : str
        Name of the region.
    internalFacesOnly : bool
        Whether to convert internal faces only (so leave boundary faces intact).
        This is only relevant if your face selection type can pick up boundary
        faces (default `True`).
    noFields : bool
        Optionally do not read/convert/write any fields (default `None`).
    """
    def __init__(
            self,
            region: str="",
            internalFacesOnly: bool=True,
            noFields: bool=None,
        ):
        super().__init__("createBafflesDict", folder="system", region=region)

        self.internalFacesOnly = internalFacesOnly
        self.noFields = noFields

        self.newBaffles: list[BaffleDict] = []

    @property
    def internalFacesOnly(self):
        return self._internalFacesOnly

    @internalFacesOnly.setter
    def internalFacesOnly(self, internalFacesOnly):
        check_type('internalFacesOnly', internalFacesOnly, bool)
        self._internalFacesOnly = internalFacesOnly

    @property
    def noFields(self):
        return self._noFields

    @noFields.setter
    def noFields(self, noFields):
        check_type('noFields', noFields, bool, none_ok=True)
        self._noFields = noFields


    def append(self, item):
        """
        Parameters
        ----------
        item : BaffleDict
        """
        check_type("item", item, BaffleDict)
        self.newBaffles.append(item)


    def add_baffle_mapped_wall(
            self,
            name: str,
            zoneName: str,
            faceType: str="faceZone",
        ):
        self.append(BaffleDict(
            name=name,
            faceType=faceType,
            zoneName=zoneName,
            masterType="mappedWall",
            slaveType="mappedWall",
            sampleMode="nearestPatchFace",
            sampleRegion=self.region,
            offsetMode="uniform",
            transform="coincidentFullMatch",
            offset=Vector(0, 0, 0)
        ))

    @property
    def is_empty(self):
        return(len(self.newBaffles) == 0)


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""

        text += addParameter("internalFacesOnly", self.internalFacesOnly, isAddExtraLine=True)

        if (self.noFields is not None):
            text += addParameter("noFields", self.noFields, isAddExtraLine=True)

        text += "baffles\n{\n"

        for newBaffle in self.newBaffles:
            text += f"{newBaffle!r}"

        text += "}\n"

        return(text)
