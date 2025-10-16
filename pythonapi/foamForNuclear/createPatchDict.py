from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import *
from foamForNuclear.openfoamFile import OpenFOAMFile


_SAMPLE_MODE_TYPES = {'nearestPatchFace'}
_CONSTRUCT_FROM_TYPES = {'patches', 'set'}


class PatchInfo:
    """
    Patch info used exclusively in createPatchDict.
    """
    def __init__(
            self,
            name: str,
            regex: str,
            patchType: str,
            inGroups: list[str]=[],
            sampleMode: str=None,
            samplePatch: str=None,
            sampleRegion: str="region0",
            constructFrom: str="patches"
        ):
        self.name = name
        self.regex = regex
        self.patchType = patchType
        self.inGroups = inGroups
        self.sampleMode = sampleMode
        self.samplePatch = samplePatch
        self.sampleRegion = sampleRegion
        self.constructFrom = constructFrom


    def __repr__(self):
        text = ""
        text += tab + "{\n"

        text += addParameter("name", self.name, indent=2)
        text += f"{2*tab}patchInfo\n"
        text += f"{2*tab}"+"{\n"

        text += addParameter("type", self.patchType, indent=3)
        if (len(self.inGroups) > 0):
            text += f"{3*tab}inGroups ({' '.join(self.inGroups)});\n"
        if (self.sampleMode is not None and self.samplePatch is not None):
            text += addParameter("sampleMode", self.sampleMode, indent=3)
            text += addParameter("sampleRegion", self.sampleRegion, indent=3)
            text += addParameter("samplePatch", self.samplePatch, indent=3)

        text += f"{2*tab}"+"}\n"

        text += addParameter("constructFrom", self.constructFrom, indent=2)
        if (self.constructFrom == "patches"):
            text += f'{2*tab}patches ("{self.regex}");\n'
        elif (self.constructFrom == "set"):
            text += addParameter("set", self.regex, indent=2)

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
    def regex(self):
        return self._regex

    @regex.setter
    def regex(self, regex):
        check_type('regex', regex, str)
        self._regex = regex

    @property
    def patchType(self):
        return self._patchType

    @patchType.setter
    def patchType(self, patchType):
        check_type('patchType', patchType, str)
        self._patchType = patchType

    @property
    def inGroups(self):
        return self._inGroups

    @inGroups.setter
    def inGroups(self, inGroups):
        check_type('inGroups', inGroups, list)
        self._inGroups = inGroups

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
    def samplePatch(self):
        return self._samplePatch

    @samplePatch.setter
    def samplePatch(self, samplePatch):
        check_type('samplePatch', samplePatch, str, none_ok=True)
        self._samplePatch = samplePatch

    @property
    def sampleRegion(self):
        return self._sampleRegion

    @sampleRegion.setter
    def sampleRegion(self, sampleRegion):
        check_type('sampleRegion', sampleRegion, str, none_ok=True)
        self._sampleRegion = sampleRegion

    @property
    def constructFrom(self):
        return self._constructFrom

    @constructFrom.setter
    def constructFrom(self, constructFrom):
        check_type('constructFrom', constructFrom, str)
        check_value('constructFrom', constructFrom, _CONSTRUCT_FROM_TYPES)
        self._constructFrom = constructFrom


class CreatePatchDict(OpenFOAMFile):
    """
    Create patches out of selected boundary faces, which are either from
    existing patches or from a faceSet.
    """
    def __init__(self, region: str=""):
        super().__init__("createPatchDict", folder="system", region=region)

        self.newPatches: list[PatchInfo] = []


    def append(self, item):
        """
        Parameters
        ----------
        item : PatchInfo
        """
        check_type("item", item, PatchInfo)
        self.newPatches.append(item)


    def add_patch(
            self,
            name: str,
            regex: str,
            patchType: str,
            inGroups: list[str]=[],
            sampleMode: str=None,
            samplePatch: str=None,
            constructFrom: str="patches"
        ):
        self.append(PatchInfo(
            name, regex, patchType, inGroups, sampleMode, samplePatch,
            sampleRegion=self.region,
            constructFrom=constructFrom
        ))

    @property
    def is_empty(self):
        return(len(self.newPatches) == 0)


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""

        text += addParameter("pointSync", 'true', isAddExtraLine=True)

        text += "patches\n(\n"

        for newPatch in self.newPatches:
            text += f"{newPatch!r}"

        text += ");\n"

        return(text)
