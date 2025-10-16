import numpy as np
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.createBafflesDict import CreateBafflesDict
from foamForNuclear.createPatchDict import CreatePatchDict
from foamForNuclear.topoSetDict import TopoSetDict

_LATTICE_TYPES =  {"hexagon", "square"}


class Mesh:
    """
    Base class of Mesh.
    """
    def __init__(self, region: str=""):
        self.region = region

        self.createPatchDict = CreatePatchDict(region=region)
        self.createBafflesDict = CreateBafflesDict(region=region)
        self.topoSetDict = TopoSetDict(region=region)


    def computeCoordinates(
            self,
            latticeType: str,
            pitch: float,
            nx: int,
            ny: int,
            xoffset: float=0,
            yoffset: float=0
        ) -> list[list[float]]:
        check_type("latticeType", latticeType, str)
        check_value("latticeType", latticeType, _LATTICE_TYPES)

        coord = []
        nRows = ny
        nCols = nx

        if (latticeType == "hexagon"):
            summitToSummit = pitch * 2/np.sqrt(3)
            for i in range(nRows):
                for j in range(nCols):
                    x = summitToSummit * 0.75 * (i - (nRows-1)/2) + xoffset
                    y = pitch * (j - (nCols-1)/2) + yoffset
                    coord.append([x, y + (i - (nRows-1)/2)*pitch/2])

        elif (latticeType == "square"):
            for i in range(nRows):
                for j in range(nCols):
                    x = pitch * (i - (nRows-1)/2) + xoffset
                    y = pitch * (j - (nCols-1)/2) + yoffset
                    coord.append([x, y])

        return(coord)


    def getLatticeAsList(
            self,
            lattice: str,
            nx: int,
            ny: int,
            pitch: float,
            latticeType: str,
            flatToFlatDirection: str='y',
            x: float=0,
            y: float=0,
            elementsToPlace: list[str] | str=None,
        ):
        """
        Parameters
        ----------
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
        """
        check_type("latticeType", latticeType, str)
        check_value("latticeType", latticeType, _LATTICE_TYPES)
        check_type("elementsToPlace", elementsToPlace, (list, str), none_ok=True)

        if (isinstance(elementsToPlace, str)):
            elementsToPlace = [elementsToPlace]

        coords = self.computeCoordinates(
            latticeType, pitch, nx, ny, xoffset=x, yoffset=y
        )

        flatLattice = (' '.join(lattice.split("\n"))).split()

        result = []
        for i, coord in enumerate(coords):
            elementName = flatLattice[i]

            if (elementsToPlace is None or elementName in elementsToPlace):
                if (flatToFlatDirection == 'y'):
                    result.append((elementName, coord[0], coord[1]))
                elif (flatToFlatDirection == 'x'):
                    result.append((elementName, coord[1], coord[0]))

        return(result)


    def mergePatchesWithName(
            self,
            name: str,
            regex: str,
            patchType: str='patch',
            inGroups: list[str]=[],
            sampleMode: str=None,
            samplePatch: str=None
        ):
        """
        Merge patches with name using `createPatch`

        Parameters
        ----------
        name : str
            Name of the new patch
        regex : str
            Regex to merge the patch in OpenFOAM style. Example:

                `".*Wall"` = all patches which ends with `Wall`
                `"core.*"` = all patches which starts with `core`
                `".*Temp.*"` = all patches which contains `Temp`

        patchType : str
            Patch type (e.g `patch`, `wall`, ...)
        """
        self.createPatchDict.add_patch(name, regex, patchType, inGroups, sampleMode, samplePatch)


    def export_to_openfoam(self):
        if (not self.createPatchDict.is_empty):
            self.createPatchDict.region = self.region
            self.createPatchDict.export_to_openfoam()

        if (not self.createBafflesDict.is_empty):
            self.createBafflesDict.region = self.region
            self.createBafflesDict.export_to_openfoam()

        if (not self.topoSetDict.is_empty):
            self.topoSetDict.region = self.region
            self.topoSetDict.export_to_openfoam()
