from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.openfoamFile import OpenFOAMFile
from foamForNuclear.common import *


_DECOMPOSE_METHOD_TYPES = {"simple", "scotch", "hierarchical", "manual"}


class DecomposeParDict(OpenFOAMFile):
    """
    Decompose a mesh and fields of a case for parallel execution.

    Parameters
    ----------
    region : str
        Name of the region to apply the decomposition.
    numberOfSubdomains : int
        Default to `1`
    method : str {"simple", "scotch", "hierarchical", "manual"}
        Method of mesh composition (default `simple`).

    Attributes
    ----------
    simpleCoeffs : OpenFOAMDict
    hierarchicalCoeffs : OpenFOAMDict
    scotchCoeffs : OpenFOAMDict
    manualCoeffs : OpenFOAMDict
    """

    def __init__(
            self,
            region: str="",
            numberOfSubdomains: int=1,
            method: str="simple"
        ):
        super().__init__(name="decomposeParDict", folder="system", region=region)

        self.numberOfSubdomains = numberOfSubdomains
        self.method = method
        self.simpleCoeffs = OpenFOAMDict({
            "n": Vector(1, 1, 1),
            "delta": 0.001
        })
        self.hierarchicalCoeffs = OpenFOAMDict({
            "n": Vector(1, 1, 1),
            "delta": 0.001,
            "order": "xyz"
        })
        self.scotchCoeffs = OpenFOAMDict({})
        self.manualCoeffs = OpenFOAMDict({
            "dataFile": None,
        })


    @property
    def numberOfSubdomains(self):
        return self._numberOfSubdomains

    @numberOfSubdomains.setter
    def numberOfSubdomains(self, numberOfSubdomains) -> None:
        check_type("numberOfSubdomains", numberOfSubdomains, int)
        if (numberOfSubdomains >= 1):
            self._numberOfSubdomains = numberOfSubdomains
        else:
            msg = f"decomposeParDict in region {self.region} must be greater or equal to 1, found {numberOfSubdomains}"
            raise ValueError(msg)

    @property
    def method(self):
        return self._method

    @method.setter
    def method(self, method) -> None:
        check_type("method", method, str)
        check_value("method", method, _DECOMPOSE_METHOD_TYPES)
        self._method = method


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""
        text += addParameter("numberOfSubdomains", self.numberOfSubdomains, isAddExtraLine=True)
        text += addParameter("method", self.method, isAddExtraLine=True)

        if (self.method == "simple"):
            text += f"simpleCoeffs{self.simpleCoeffs!r}\n"
        elif (self.method == "scotch"):
            text += f"scotchCoeffs{self.scotchCoeffs!r}\n"
        elif (self.method == "hierarchical"):
            text += f"hierarchicalCoeffs{self.hierarchicalCoeffs!r}\n"
        elif (self.method == "manual"):
            text += f"manualCoeffs{self.manualCoeffs!r}\n"

        return(text)


    def extract_info_from_log(self, caseFolder: str="./"):
        info = {}

        filename = f"{caseFolder}/log.decomposePar"
        with open(filename, 'r') as f:
            for line in f:
                if ("Max number of cells" in line):
                    info["nMaxCells"] = int(line.split()[5])
                elif ("Max number of processor patches" in line):
                    info["nMaxProcessorPatches"] = int(line.split()[6])
                elif ("Max number of faces between processors" in line):
                    info["nMaxProcessorPatches"] = int(line.split()[7])

        return(info)
