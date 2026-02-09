from __future__ import annotations

import attrs as attr
from attrs import define, field, validators as v
from foamForNuclear._attrs_tools import auto_type_validator, call_method_on_change, _propagate_region_to

from abc import abstractmethod
import os

from foamForNuclear.boundaryConditions.zeroGradient import ZeroGradient
from foamForNuclear.common import *
from foamForNuclear.checkvalue import check_type, CheckedList
from foamForNuclear.mesh.dicts import DecomposeParDict
from foamForNuclear.mesh.dicts import DynamicMeshDict
from foamForNuclear.fields import Field
from foamForNuclear.numerics import fvSchemes, fvSolution
from foamForNuclear.mesh.mesh import Mesh
from foamForNuclear.preprocessing import SetFieldRegion, SetFieldsDict
from foamForNuclear.timeFolder import TimeFolder

@define(
    slots=True,
    on_setattr=[attr.setters.validate, call_method_on_change("propagate_region", "region")],
    field_transformer=auto_type_validator
)
class Solver:
    """
    Base class for solver.

    Parameters
    ----------
    region : str
        Name of the region
    solver : str
        Name of the solver to be used in the region (e.g `diffusionNeutronics`, `onePhase`, ...)
    removeBaffles : bool
        Flag to remove baffles. Create a ghost region without baffles. Can help
        if solving using thermal-hydraulics with porous media. (default `False`)
    timeFolder : TimeFolder
        Object representing the time folder
    mesh : Mesh
        Mesh object
    isMeshDeformation : bool
        Flag to allow mesh deformation based on a displacement field (default `False`)
    displacementFieldName : str
        Name of the displacement field (default `disp`)

    Attributes
    ----------
    region : str
        Name of the region
    solver : str
        Name of the solver to be used in the region (e.g `diffusionNeutronics`, `onePhase`, ...)
    removeBaffles : bool
        Flag to remove baffles. Create a ghost region without baffles. Can help
        if solving using thermal-hydraulics with porous media.
    timeFolder : TimeFolder
        Object representing the time folder
    mesh : Mesh
        Mesh object
    fvSchemes : fvSchemes
        Select schemes and solving strategy (see https://www.openfoam.com/documentation/user-guide/6-solving/6.2-numerical-schemes)
    fvSolution : fvSolution
        Select solution, residual and tolerances (see https://www.openfoam.com/documentation/user-guide/6-solving/6.3-solution-and-algorithm-control)
    decomposeParDict : DecomposeParDict
        Mesh decomposition dict for parallel computing.
    isMeshDeformation : bool
        Flag to allow mesh deformation based on a displacement field (default `False`)
    displacementFieldName : str
        Name of the displacement field (default `disp`).
    setFieldsDict : SetFieldsDict
        setFieldsDict
    dynamicMeshDict : DynamicMeshDict
        Dynamic mesh motion dict.
    """

    region: str = field(default="")
    solver: str = field(default="none")
    removeBaffles: bool = False
    timeFolder: TimeFolder | None = None
    mesh: Mesh | None = None
    isMeshDeformation: bool = False
    displacementFieldName: str = "disp"

    fvSchemes: fvSchemes = field(factory=fvSchemes)
    fvSolution: fvSolution = field(factory=fvSolution)
    decomposeParDict: DecomposeParDict = field(factory=DecomposeParDict)
    setFieldsDict: SetFieldsDict = field(factory=SetFieldsDict)
    dynamicMeshDict: DynamicMeshDict = field(
        factory=lambda: DynamicMeshDict(dynamicFvMesh="dynamicMotionSolverFvMesh"))

    def propagate_region(self, new_value: str | None = None) -> None:
        if new_value is None:
            new_value = self.region
        _propagate_region_to(self, new_value, "fvSchemes", "fvSolution")

    def __attrs_post_init__(self):
        _propagate_region_to(self, self.region, "fvSchemes", "fvSolution")

    # @region.on_setattr
    # def _on_region_change(self, _attr, _value):
    #     propagate_region_to(self, "fvSchemes", "fvSolution")

    def __repr__(self, depth = 0):
        text = ""
        text += f"{depth * tab}{type(self).__name__}\n"
        text += f"{(depth+1)*tab}Solver: {self.solver}\n"
        text += f"{(depth+1)*tab}Region: {self.region}\n"

        requiredFields = ', '.join(self.get_required_fields())

        if (len(self.get_required_fields()) > 0):
            text += f"{(depth+1)*tab}Required fields: {requiredFields}\n"

        text += f"{(depth+1)*tab}Mesh generation: "
        if (self.mesh is not None):
            text += f"{type(self.mesh).__name__}\n"
        else:
            text += "unknown\n"

        if (self.removeBaffles):
            text += f"{(depth+1)*tab}Remove baffles\n"
        else:
            text += f"{(depth+1)*tab}Keep baffles\n"

        if (self.isMeshDeformation):
            text += f"{(depth+1)*tab}Mesh deformation from '{self.displacementFieldName}' field\n"

        if (self.decomposeParDict.numberOfSubdomains > 1):
            text += f"{(depth+1)*tab}Parallel computing over {self.decomposeParDict.numberOfSubdomains} subdomains\n"
            text += f"{(depth+1)*tab}Mesh decomposition using '{self.decomposeParDict.method}' method\n"

        if (not self.dynamicMeshDict.is_empty):
            text += f"{(depth+1)*tab}Dynamic motion activated\n"
            text += f"{self.dynamicMeshDict.__repr__(depth=depth+2)}"

        return(text)


    # @property
    # def region(self):
    #     return self._region

    # @region.setter
    # def region(self, region) -> None:
    #     if region is not None:
    #         check_type("region", region, str)
    #         self._region = region
    #     else:
    #         self._region = ''

    # @property
    # def solver(self):
    #     return self._solver

    # @solver.setter
    # def solver(self, solver) -> None:
    #     if solver is not None:
    #         check_type("solver", solver, str)
    #         self._solver = solver
    #     else:
    #         self._solver = ''

    # @fvSchemes.setter
    # def fvSchemes(self, fvSchemes_) -> None:
    #     check_type("fvSchemes", fvSchemes_, fvSchemes)
    #     self._fvSchemes = fvSchemes_
    #     self._fvSchemes.region = self.region

    # @fvSolution.setter
    # def fvSolution(self, fvSolution_) -> None:
    #     check_type("fvSolution", fvSolution_, fvSolution)
    #     self._fvSolution = fvSolution_
    #     self._fvSolution.region = self.region


    def create_folders(self) -> None:
        if (self.region is not None and self.region != ""):
            if (not os.path.exists(f"constant/{self.region}")):
                os.mkdir(f"constant/{self.region}")
            if (not os.path.exists(f"system/{self.region}")):
                os.mkdir(f"system/{self.region}")


    def add_relaxation_on_equation(self, equationName: str, value: float) -> None:
        self.fvSolution.add_relaxation_on_equation(equationName, value)


    def add_relaxation_on_field(self, fieldName: str, value: float) -> None:
        self.fvSolution.add_relaxation_on_field(fieldName, value)


    def add_initial_values(self, item: SetFieldRegion):
        check_type("item: SetFieldRegion", item, SetFieldRegion)
        self.setFieldsDict.append(item)


    @abstractmethod
    def get_required_fields(self) -> list[str]:
        """
        Return a list of required fields depending on the solver selected.
        They are labeled "MUST_READ"
        """
        return([])

    @abstractmethod
    def get_default_fields(self) -> list[str]:
        """
        Return a list of default fields that can be used as input depending on
        the solver selected. They are labeled "READ_IF_PRESENT"
        """
        return([])


    def create_zone_field(
            self,
            cellZones: list[str]
        ) -> Field:
        """
        Create a dummy field where each cellZone is assign a unique id. This
        field is never used in a calculation, it is purely for visualisation
        purpose. By default, each cell is set to index 0. The created field
        needs to be manually attached to the :class:`TimeFolder`.

        Parameters
        ----------
        cellZones : list[str]
            List of the zone names.

        Return
        ------
        Return a field named `zoneIdx` with zeroGradient BC everywhere.
        """
        idxField = Field("zoneIdx", internalField=0, region=self.region)
        idxField.set_boundary_condition('".*"', ZeroGradient())

        for idx, cellZone in enumerate(cellZones):
            self.setFieldsDict.add_zone_to_cell([(idxField, idx+1)], cellZone)

        return(idxField)


    def export_to_openfoam(self):
        self.create_folders()

        self.fvSchemes.region = self.region
        self.fvSolution.region = self.region
        self.decomposeParDict.region = self.region

        if (self.mesh is not None): # and isinstance(self.mesh, (BlockMesh, PolyMesh))):
            self.mesh.region = self.region
            self.mesh.export_to_openfoam()

        if (not self.setFieldsDict.is_empty):
            self.setFieldsDict.region = self.region
            self.setFieldsDict.export_to_openfoam()

        self.fvSchemes.export_to_openfoam()
        self.fvSolution.export_to_openfoam()

        if (self.decomposeParDict.numberOfSubdomains > 1):
            self.decomposeParDict.export_to_openfoam()

        if (self.dynamicMeshDict is not None and not self.dynamicMeshDict.is_empty):
            self.dynamicMeshDict.region = self.region
            self.dynamicMeshDict.export_to_openfoam()



class Solvers(CheckedList):
    """
    Main solver collection object.
    """
    def __init__(self, solvers: list[Solver]=None):
        super().__init__(Solver, 'solver collection', solvers)


    def set_time_folder(self, timeFolder: TimeFolder):
        for solver in self:
            solver.timeFolder = timeFolder


    def set_solver_in_parallel_settings(self):
        maxSubdomains = self.get_max_subdomains()

        if (maxSubdomains > 1):
            for solver in self:
                if (solver.decomposeParDict.numberOfSubdomains < maxSubdomains):
                    print(f"{Keyword.WARNING}: Set {solver.region} solver to {maxSubdomains} subdomains")
                    print(f"{tab}Found {solver.decomposeParDict.numberOfSubdomains} subdomains instead of {maxSubdomains}")
                    print(f"{tab}Set decompose method to 'scotch'")
                    solver.decomposeParDict.numberOfSubdomains = maxSubdomains
                    solver.decomposeParDict.method = 'scotch'


    def get_max_subdomains(self) -> int:
        if (len(self) > 0):
            return(max([
                solver.decomposeParDict.numberOfSubdomains
                for solver in self
            ]))

        return(1)


    def is_one_solver_parallel(self) -> bool:
        return(any([
            solver.decomposeParDict.numberOfSubdomains > 1
            for solver in self
        ]))


    def write_to_controlDict(self):
        text = "removeBaffles\n"
        text += "{\n"
        for solver in self:
            if (isinstance(solver, Solver) and solver.region != ""):
                text += f"{tab}{solver.region:15} {1 if solver.removeBaffles else 0};\n"
        text += "}\n\n"
        return(text)


    def export_region_properties_to_openfoam(self):
        if (not os.path.exists("constant")):
            os.mkdir("constant")

        with open("constant/regionProperties", 'w') as f:
            f.write(openfoamHeader)
            f.write(openfoamFileHeader("regionProperties"))

            f.write("regions\n")
            f.write("(\n")
            for solver in self:
                f.write(f"{tab}{solver.region:23} ({solver.region})\n")
            f.write(");\n\n")

            f.write(openfoamFooterLine)


    def export_to_openfoam(self):
        self.set_solver_in_parallel_settings()

        for solver in self:
            solver.export_to_openfoam()

        self.export_region_properties_to_openfoam()
