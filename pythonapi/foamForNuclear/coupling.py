from foamForNuclear.common import *
from foamForNuclear.checkvalue import check_positive, check_type, check_value
from foamForNuclear.solver import Solver


_MULTIPHYSICS_LOOP_TYPES = {
    "picardLoop", "picardLoopNoFluid", "FSILoop", "CHTLoop", "multiScaleLoop"
}


class FieldTransfer:
    """
    Field transfer object that stores src/destination regions and mapped fields.

    Parameters
    ----------
    sourceRegion : Solver
        Source solver to map from.
    destinationRegion : Solver
        Destination solver to map to.
    sourceField : str
        Source field name to map data from.
    destinationField : str
        Destination field name to map data to.
    """

    def __init__(
            self,
            sourceRegion: Solver,
            destinationRegion: Solver,
            sourceField: str,
            destinationField: str
        ):
        self.sourceRegion: Solver = sourceRegion
        self.destinationRegion: Solver = destinationRegion
        self.sourceField = sourceField
        self.destinationField = destinationField


    def __repr__(self):
        return(f"{self.sourceField} ({self.sourceRegion.region}) -> {self.destinationField} ({self.destinationRegion.region})")


    @property
    def sourceRegion(self):
        return self._sourceRegion

    @sourceRegion.setter
    def sourceRegion(self, sourceRegion) -> None:
        check_type("sourceRegion", sourceRegion, Solver, none_ok=True)
        self._sourceRegion = sourceRegion

    @property
    def destinationRegion(self):
        return self._destinationRegion

    @destinationRegion.setter
    def destinationRegion(self, destinationRegion) -> None:
        check_type("destinationRegion", destinationRegion, Solver, none_ok=True)
        self._destinationRegion = destinationRegion

    @property
    def sourceField(self):
        return self._sourceField

    @sourceField.setter
    def sourceField(self, sourceField) -> None:
        check_type("sourceField", sourceField, str, none_ok=True)
        self._sourceField = sourceField

    @property
    def destinationField(self):
        return self._destinationField

    @destinationField.setter
    def destinationField(self, destinationField) -> None:
        check_type("destinationField", destinationField, str, none_ok=True)
        self._destinationField = destinationField



class MultiPhysicsLoop(CheckedList):
    """
    Multi-physics loop applying multiple iterations over multiple coupled
    solvers:

    `picardLoop`: This class creates a multiphysics tightly coupled loop, made
    of an arbitrary number of physics. The class member-functions are
    inherited from solver.H and are the same as those of the other solvers.
    Each member function consists of a loop where the same member function
    is called iteratively for all the solvers included in the `picardLoop`.
    Note that another `picardLoop` can be created within a `picardLoop`,
    hence creating a hierarchical structure which might be helpful for
    multi-scale problems. The `picardLoop` automatically solves all the
    equations of the solvers for each time step. In some cases, it is
    convenient to only solve the energy equation of the thermal-hydraulics
    solver, while the fluid-mechanics is solved only once at the beginning
    of the iteration. This can be done by setting the `picardLoopNoFluid`
    solver instead.

    `picardLoopNoFluid`: Identical to the `picardLoop`, but it does not solve
    the fluid-mechancs equations at each of the iterations. Instead, it
    solves the fluid-mechanics equations only once at the beginning of the
    iteration, and then solves the energy equation of the thermal-hydraulics
    solver iteratively with the other physics.


    Parameters
    ----------
    solver : str {"picardLoop", "picardLoopNoFluid", "FSILoop", "CHTLoop", "multiScaleLoop"}
        Multiphysics loop type.
    region : str
        Name of the multi-physics loop, e.g `Level_1`.
    maxResidual : float
        Minimum residual after which the iteration is ended.
    maxIterations : int
        Maximum number of iterations before the loop is ended.
    couplingStartTime : float
        Default `0`
    solvers : list[Solver | MultiPhysicsLoop]
        List of solver/multi-physics loop to iter on
    """

    def __init__(
            self,
            solver: str,
            region: str,
            maxResidual: float,
            maxIterations: int,
            couplingStartTime: float=None,
            includeFluidMechanicsInLoop: bool=None,
            solvers: list[Solver]=None
        ):
        super().__init__((Solver, MultiPhysicsLoop), "MultiPhysicsLoop", solvers)

        self.region = region
        self.solver = solver
        self.maxResidual = maxResidual
        self.maxIterations = maxIterations
        self.couplingStartTime = couplingStartTime
        self.includeFluidMechanicsInLoop = includeFluidMechanicsInLoop

    def __repr__(self):
        text = f"{tab}{self.region}\n"
        text += tab + "{\n"
        text += f"{2*tab}subSolvers\n"
        text += 2*tab + "{\n"
        for solver in self:
            text += addParameter(solver.region, solver.solver, indent=3)
        text += 2*tab + "}\n"
        text += addParameter('maxResidual', self.maxResidual, indent=2)
        text += addParameter('maxIterations', self.maxIterations, indent=2)
        text += addParameter('couplingStartTime', self.couplingStartTime, indent=2, none_ok=False)
        text += addParameter('includeFluidMechanicsInLoop', self.includeFluidMechanicsInLoop, indent=2, none_ok=False)
        text += tab + "}\n"
        return(text)

    @property
    def region(self):
        return self._region

    @region.setter
    def region(self, region) -> None:
        check_type("region", region, str)
        self._region = region

    @property
    def solver(self):
        return self._solver

    @solver.setter
    def solver(self, solver) -> None:
        check_type("solver", solver, str)
        check_value("solver", solver, _MULTIPHYSICS_LOOP_TYPES)
        self._solver = solver

    @property
    def maxResidual(self):
        return self._maxResidual

    @maxResidual.setter
    def maxResidual(self, maxResidual) -> None:
        check_type("maxResidual", maxResidual, (float, int))
        check_positive("maxResidual", maxResidual)
        self._maxResidual = maxResidual

    @property
    def maxIterations(self):
        return self._maxIterations

    @maxIterations.setter
    def maxIterations(self, maxIterations) -> None:
        check_type("maxIterations", maxIterations, int)
        check_positive("maxIterations", maxIterations)
        self._maxIterations = maxIterations

    @property
    def couplingStartTime(self):
        return self._couplingStartTime

    @couplingStartTime.setter
    def couplingStartTime(self, couplingStartTime) -> None:
        check_type("couplingStartTime", couplingStartTime, (float, int), none_ok=True)
        self._couplingStartTime = couplingStartTime

    @property
    def includeFluidMechanicsInLoop(self):
        return self._includeFluidMechanicsInLoop

    @includeFluidMechanicsInLoop.setter
    def includeFluidMechanicsInLoop(self, includeFluidMechanicsInLoop) -> None:
        check_type("includeFluidMechanicsInLoop", includeFluidMechanicsInLoop, bool, none_ok=True)
        self._includeFluidMechanicsInLoop = includeFluidMechanicsInLoop


class CHTLoop(MultiPhysicsLoop):
    """
    Conjugate Heat Transfer (CHT) multi-physics loop.

    Parameters
    ----------
    maxResidual : float
        Minimum residual after which the iteration is ended.
    maxIterations : int
        Maximum number of iterations before the loop is ended.
    useHTC : bool
        Default `False`
    oneWayCoupling : bool
        Default `False`
    fluidKappa : str
        Default `kappaEff`
    solidKappa : str
        Default `k`
    """
    def __init__(
            self,
            region: str,
            maxResidual: float,
            maxIterations: int,
            fluidRegionName: str,
            solidRegionName: str,
            fluidPatches: list[str],
            solidPatches: list[str],
            useHTC: bool=None,
            oneWayCoupling: bool=None,
            fluidKappa: str=None,
            solidKappa: str=None,
            couplingStartTime: float=None,
            includeFluidMechanicsInLoop: bool=None,
            solvers: list[Solver | MultiPhysicsLoop]=None
        ):
        super().__init__(
            "CHTLoop",
            region,
            maxResidual,
            maxIterations,
            couplingStartTime,
            includeFluidMechanicsInLoop,
            solvers
        )
        self.fluidRegionName = fluidRegionName
        self.solidRegionName = solidRegionName
        self.fluidPatches = fluidPatches
        self.solidPatches = solidPatches
        self.useHTC = useHTC
        self.oneWayCoupling = oneWayCoupling
        self.fluidKappa = fluidKappa
        self.solidKappa = solidKappa

    def __repr__(self):
        text = super().__repr__()
        text = text[:-6] # Remove last "tab}\n"

        text += "\n"
        text += f"{2*tab}CHTLoopCoeffs\n"
        text += 2*tab + "{\n"
        text += addParameter('fluidRegionName', self.fluidRegionName, indent=3)
        text += addParameter('solidRegionName', self.solidRegionName, indent=3)
        text += addParameter('fluidPatches', self.fluidPatches, indent=3)
        text += addParameter('solidPatches', self.solidPatches, indent=3)
        text += addParameter('useHTC', self.useHTC, indent=3, none_ok=False)
        text += addParameter('oneWayCoupling', self.oneWayCoupling, indent=3, none_ok=False)
        text += addParameter('couplingStartTime', self.couplingStartTime, indent=3, none_ok=False)
        text += addParameter('fluidKappa', self.fluidKappa, indent=3, none_ok=False)
        text += addParameter('solidKappa', self.solidKappa, indent=3, none_ok=False)
        text += 2*tab + "}\n"
        text += tab + "}\n"

        return text

    @property
    def fluidRegionName(self):
        return self._fluidRegionName

    @fluidRegionName.setter
    def fluidRegionName(self, fluidRegionName) -> None:
        check_type("fluidRegionName", fluidRegionName, str)
        self._fluidRegionName = fluidRegionName

    @property
    def solidRegionName(self):
        return self._solidRegionName

    @solidRegionName.setter
    def solidRegionName(self, solidRegionName) -> None:
        check_type("solidRegionName", solidRegionName, str)
        self._solidRegionName = solidRegionName

    @property
    def fluidPatches(self):
        return self._fluidPatches

    @fluidPatches.setter
    def fluidPatches(self, fluidPatches) -> None:
        check_type("fluidPatches", fluidPatches, (list, List))
        if (isinstance(fluidPatches, list)):
            self._fluidPatches = List(fluidPatches)
        else:
            self._fluidPatches = fluidPatches

    @property
    def solidPatches(self):
        return self._solidPatches

    @solidPatches.setter
    def solidPatches(self, solidPatches) -> None:
        check_type("solidPatches", solidPatches, (list, List))
        if (isinstance(solidPatches, list)):
            self._solidPatches = List(solidPatches)
        else:
            self._solidPatches = solidPatches

    @property
    def useHTC(self):
        return self._useHTC

    @useHTC.setter
    def useHTC(self, useHTC) -> None:
        check_type("useHTC", useHTC, bool, none_ok=True)
        self._useHTC = useHTC

    @property
    def oneWayCoupling(self):
        return self._oneWayCoupling

    @oneWayCoupling.setter
    def oneWayCoupling(self, oneWayCoupling) -> None:
        check_type("oneWayCoupling", oneWayCoupling, bool, none_ok=True)
        self._oneWayCoupling = oneWayCoupling

    @property
    def fluidKappa(self):
        return self._fluidKappa

    @fluidKappa.setter
    def fluidKappa(self, fluidKappa) -> None:
        check_type("fluidKappa", fluidKappa, str, none_ok=True)
        self._fluidKappa = fluidKappa

    @property
    def solidKappa(self):
        return self._solidKappa

    @solidKappa.setter
    def solidKappa(self, solidKappa) -> None:
        check_type("solidKappa", solidKappa, str, none_ok=True)
        self._solidKappa = solidKappa


class FSILoop(MultiPhysicsLoop):
    """
    Fluid-Structure Interaction (FSI) multi-physics loop.

    Parameters
    ----------
    maxResidual : float
        Minimum residual after which the iteration is ended.
    maxIterations : int
        Maximum number of iterations before the loop is ended.
    fluidRegionName : str
        Name of the fluid region.
    solidRegionName : str
        Name of the solid region.
    fluidSidePatchName : str
        Name of the coupled patch on the solid side.
    solidSidePatchName : str
        Name of the coupled patch on the fluid side.
    underRelaxationFactor : float
        Under-relaxation factor for the coupling.
    couplingStartTime : float
        Time after which the coupling starts (default `0`).
    thermalCoupling : bool
        Default `False`.
    """
    def __init__(
            self,
            region: str,
            maxResidual: float,
            maxIterations: int,
            FSIInterface: bool,
            fluidRegionName: str,
            solidRegionName: str,
            fluidSidePatchName: str,
            solidSidePatchName: str,
            underRelaxationFactor: float,
            thermalCoupling: bool=None,
            useHTC: bool=None,
            fluidKappa: str=None,
            solidKappa: str=None,
            couplingStartTime: float=None,
            includeFluidMechanicsInLoop: bool=None,
            solvers: list[Solver | MultiPhysicsLoop]=None
        ):
        super().__init__(
            "FSILoop",
            region,
            maxResidual,
            maxIterations,
            couplingStartTime,
            includeFluidMechanicsInLoop,
            solvers
        )
        self.FSIInterface = FSIInterface
        self.fluidRegionName = fluidRegionName
        self.solidRegionName = solidRegionName
        self.fluidSidePatchName = fluidSidePatchName
        self.solidSidePatchName = solidSidePatchName
        self.underRelaxationFactor = underRelaxationFactor
        self.thermalCoupling = thermalCoupling
        self.useHTC = useHTC
        self.fluidKappa = fluidKappa
        self.solidKappa = solidKappa

    def __repr__(self):
        text = super().__repr__()
        text = text[:-6] # Remove last "tab}\n"

        text += "\n"
        text += f"{2*tab}FSILoopCoeffs\n"
        text += 2*tab + "{\n"
        text += addParameter('FSIInterface', self.FSIInterface, indent=3)
        text += addParameter('fluidRegionName', self.fluidRegionName, indent=3)
        text += addParameter('solidRegionName', self.solidRegionName, indent=3)
        text += addParameter('fluidSidePatchName', self.fluidSidePatchName, indent=3)
        text += addParameter('solidSidePatchName', self.solidSidePatchName, indent=3)
        text += addParameter('underRelaxationFactor', self.underRelaxationFactor, indent=3)
        text += addParameter('couplingStartTime', self.couplingStartTime, indent=3, none_ok=False)
        text += addParameter('thermalCoupling', self.thermalCoupling, indent=3, none_ok=False)
        text += addParameter('useHTC', self.useHTC, indent=3, none_ok=False)
        text += addParameter('fluidKappa', self.fluidKappa, indent=3, none_ok=False)
        text += addParameter('solidKappa', self.solidKappa, indent=3, none_ok=False)
        text += 2*tab + "}\n"
        text += tab + "}\n"

        return text

    @property
    def FSIInterface(self):
        return self._FSIInterface

    @FSIInterface.setter
    def FSIInterface(self, FSIInterface) -> None:
        check_type("FSIInterface", FSIInterface, bool)
        self._FSIInterface = FSIInterface

    @property
    def fluidRegionName(self):
        return self._fluidRegionName

    @fluidRegionName.setter
    def fluidRegionName(self, fluidRegionName) -> None:
        check_type("fluidRegionName", fluidRegionName, str)
        self._fluidRegionName = fluidRegionName

    @property
    def solidRegionName(self):
        return self._solidRegionName

    @solidRegionName.setter
    def solidRegionName(self, solidRegionName) -> None:
        check_type("solidRegionName", solidRegionName, str)
        self._solidRegionName = solidRegionName

    @property
    def fluidSidePatchName(self):
        return self._fluidSidePatchName

    @fluidSidePatchName.setter
    def fluidSidePatchName(self, fluidSidePatchName) -> None:
        check_type("fluidSidePatchName", fluidSidePatchName, str)
        self._fluidSidePatchName = fluidSidePatchName

    @property
    def solidSidePatchName(self):
        return self._solidSidePatchName

    @solidSidePatchName.setter
    def solidSidePatchName(self, solidSidePatchName) -> None:
        check_type("solidSidePatchName", solidSidePatchName, str)
        self._solidSidePatchName = solidSidePatchName

    @property
    def underRelaxationFactor(self):
        return self._underRelaxationFactor

    @underRelaxationFactor.setter
    def underRelaxationFactor(self, underRelaxationFactor) -> None:
        check_type("underRelaxationFactor", underRelaxationFactor, (int, float))
        self._underRelaxationFactor = underRelaxationFactor

    @property
    def thermalCoupling(self):
        return self._thermalCoupling

    @thermalCoupling.setter
    def thermalCoupling(self, thermalCoupling) -> None:
        check_type("thermalCoupling", thermalCoupling, bool, none_ok=True)
        self._thermalCoupling = thermalCoupling

    @property
    def useHTC(self):
        return self._useHTC

    @useHTC.setter
    def useHTC(self, useHTC) -> None:
        check_type("useHTC", useHTC, bool, none_ok=True)
        self._useHTC = useHTC

    @property
    def fluidKappa(self):
        return self._fluidKappa

    @fluidKappa.setter
    def fluidKappa(self, fluidKappa) -> None:
        check_type("fluidKappa", fluidKappa, str, none_ok=True)
        self._fluidKappa = fluidKappa

    @property
    def solidKappa(self):
        return self._solidKappa

    @solidKappa.setter
    def solidKappa(self, solidKappa) -> None:
        check_type("solidKappa", solidKappa, str, none_ok=True)
        self._solidKappa = solidKappa


class Coupling(CheckedList):
    """
    Coupling is a list of solvers/multi-physics loops.

    Parameters
    ----------
    solvers : list[Solver | MultiPhysicsLoop]
        List of solvers/multi-physics loops

    Attributes
    ----------
    mappings : list[FieldTransfer]
        List of field transfer
    """

    def __init__(
            self,
            solvers: list[Solver | MultiPhysicsLoop]=None
        ):
        super().__init__((Solver, MultiPhysicsLoop), "Coupling", solvers)

        self.mappings: list[FieldTransfer] = []


    def __repr__(self):
        text = underline("Coupling map:") + "\n"
        text += self.get_coupling_map_as_text()

        text += underline("Coupling loop:") + "\n"
        text += self.get_coupling_loop_as_text(self, depth=1)
        return(text)


    def add_field_transfer(
            self,
            sourceRegion: Solver,
            destinationRegion: Solver,
            sourceField: str,
            destinationField: str
        ):
        """
        Add field transfer link between 2 solvers/regions. The sourceField from
        the sourceRegion is mapped onto the destinationRegion, overwriting the
        destinationField.
        """
        check_type("sourceRegion", sourceRegion, Solver)
        check_type("destinationRegion", destinationRegion, Solver)
        check_type("sourceField", sourceField, str)
        check_type("destinationField", destinationField, str)
        self.mappings.append(FieldTransfer(
            sourceRegion,
            destinationRegion,
            sourceField,
            destinationField
        ))


    def get_coupling_map_as_text(self) -> str:
        text = f"{tab}None\n"
        if (len(self.mappings) > 0):
            text = ""
            regionNames = self.get_unique_regions(self)
            for srcRegionName in regionNames:
                for destRegionName in regionNames:
                    if (destRegionName == srcRegionName):
                        continue

                    reducedMapping = [
                        fieldTransfer for fieldTransfer in self.mappings
                        if
                            fieldTransfer.sourceRegion.region == srcRegionName
                            and
                            fieldTransfer.destinationRegion.region == destRegionName
                    ]
                    if (len(reducedMapping) == 0):
                        continue

                    nn = max([len(fieldTransfer.sourceField) for fieldTransfer in reducedMapping]) + 1

                    text += f"{tab}{srcRegionName} -> {destRegionName}\n"
                    for i, fieldTransfer in enumerate(reducedMapping):
                        lc = "└─" if i+1 == len(reducedMapping) else "├─"
                        text += f"{2*tab}{lc} {fieldTransfer.sourceField+' ':-<{nn}}-> {fieldTransfer.destinationField}\n"

        return(text)


    def get_coupling_loop_as_text(
            self,
            solvers: list[Solver | MultiPhysicsLoop],
            depth: int=0
        ) -> str:
        """
        Recursive method
        """
        text = ""
        for i, solver in enumerate(solvers):
            lc = "├─"
            if (i+1 == len(solvers) and isinstance(solver, Solver)):
                lc = "└─"
            elif (isinstance(solver, MultiPhysicsLoop)):
                lc = "├⥁" # Represents a loop

            prefix = f"{'│'.join(depth*[tab])}{lc} "

            text += f"{prefix}{solver.region}: {bold(solver.solver)}"
            if (isinstance(solver, Solver)):
                text += "\n"
            if (isinstance(solver, MultiPhysicsLoop)):
                text += f" (maxRes={solver.maxResidual}, maxIter={solver.maxIterations})\n"
                text += self.get_coupling_loop_as_text(solver, depth=depth+1)

        return(text)


    def get_unique_regions(self, solvers: list[Solver | MultiPhysicsLoop]) -> list[str]:
        uniqueRegions = []

        for solver in solvers:
            if (isinstance(solver, Solver) and solver.region not in uniqueRegions):
                uniqueRegions.append(solver.region)
            if (isinstance(solver, MultiPhysicsLoop)):
                uniqueRegions += self.get_unique_regions(solver)

        return(uniqueRegions)


    def get_unique_source_regions(self, destinationRegion: str) -> list[str]:
        uniqueRegions = []
        for link in self.mappings:
            if (
                link.destinationRegion.region == destinationRegion
                and link.sourceRegion.region not in uniqueRegions
            ):
                uniqueRegions.append(link.sourceRegion.region)
        return(uniqueRegions)


    def get_field_transfers(self, sourceRegion: str, destinationRegion: str) -> list[FieldTransfer]:
        return([
            link for link in self.mappings
                if link.destinationRegion.region == destinationRegion
                and link.sourceRegion.region == sourceRegion
        ])


    def write_multi_physics_solvers(self, solvers: list[Solver | MultiPhysicsLoop]) -> str:
        """
        Recursive writing
        """
        text = ""
        for solver in solvers:
            if (isinstance(solver, MultiPhysicsLoop)):
                text += f"{solver!r}"
                text += self.write_multi_physics_solvers(solver)
        return(text)


    def write_mesh_deformation(self, solvers: list[Solver | MultiPhysicsLoop]) -> str:
        """
        Recursive writing
        """
        text = ""
        for solver in solvers:
            if (isinstance(solver, Solver) and solver.isMeshDeformation):
                text += f"{tab}{solver.region}\n"
                text += tab + "{\n"
                text += f"{2*tab}{'displacementField':15} {solver.displacementFieldName};\n"
                text += tab + "}\n"
            elif (isinstance(solver, MultiPhysicsLoop)):
                text += self.write_mesh_deformation(solver)
        return(text)


    def write_mesh_deformation_to_openfoam(self) -> str:
        text = "meshDeformation\n{\n"
        text += self.write_mesh_deformation(self)
        text += "}\n\n"
        return(text)


    def write_mapping_to_openfoam(self) -> str:
        text = "mappings\n{\n"

        uniqueRegions = self.get_unique_regions(self)

        for destinationRegion in uniqueRegions:
            text += f"{tab}{destinationRegion}\n"
            text += tab+"{\n"

            uniqueSourceRegions = self.get_unique_source_regions(destinationRegion)
            for sourceRegion in uniqueSourceRegions:
                text += f"{2*tab}{sourceRegion} // -> {destinationRegion}\n"
                text += 2*tab+"{\n"

                sourceFields, targetFields = "", ""
                for link in self.get_field_transfers(sourceRegion, destinationRegion):
                    sourceFields += link.sourceField + " "
                    targetFields += link.destinationField + " "

                text += f"{3*tab}{'sourceFields':15} ( {sourceFields});\n"
                text += f"{3*tab}{'targetFields':15} ( {targetFields});\n"

                text += 2*tab+"}\n"

            text += tab+"}\n"

        text +=  "}\n"

        return(text)


    def write_level0_solvers(self) -> str:
        text = "regionSolvers\n"
        text += "{\n"
        text += tab + "Level_0\n"
        text += tab +"{\n"
        for solver in self:
            text += f"{2*tab}{solver.region:15} {solver.solver};\n"
        text += tab + "}\n"
        text += self.write_multi_physics_solvers(self)
        text += "}\n\n"
        return(text)


    def export_to_openfoam(self):
        with open(f"system/regionsDict", 'w') as f:
            f.write(openfoamHeader)
            f.write(openfoamFileHeader("regionsDict"))

            f.write(self.write_level0_solvers())
            f.write(self.write_mesh_deformation_to_openfoam())
            f.write(self.write_mapping_to_openfoam())

            f.write(openfoamFooterLine)


    def plot_coupling_graph(self):
        """
        Plot the coupling graph conneting solvers and field transfered.
        """
        from python_mermaid.diagram import MermaidDiagram, Node, Link

        nodes = {e: Node(e) for e in self.get_unique_regions(self)}

        links = [Link(
            nodes[link.sourceRegion.region],
            nodes[link.destinationRegion.region],
            message=link.sourceField
        ) for link in self.mappings]

        chart = MermaidDiagram(
            title="Coupling",
            nodes=list(nodes.values()),
            links=links
        )

        print(f"Generate coupling graph ... ", end="")
        status = generate_mermaid_graph_as_image(f"{chart}", "fig_graph_coupling.png")
        print(Keyword.DONE if status else Keyword.ERROR)

        return(chart)


    def plot_solving_graph(self):
        """
        Plot solving graph connecting solvers and multi-physics loops as a tree.
        """
        from python_mermaid.diagram import MermaidDiagram, Node, Link

        def extract(l, loopNode: str='root'):
            nodes, links = {loopNode: Node(loopNode, shape='circle')}, []
            for solver in l:
                if (isinstance(solver, Solver) and solver.region not in nodes):
                    nodes[solver.region] = Node(solver.region)
                if (isinstance(solver, MultiPhysicsLoop)):
                    newNodes, newLinks = extract(solver, loopNode=solver.region)
                    nodes = nodes | newNodes
                    links += newLinks

                links.append(Link(
                    nodes[loopNode], nodes[solver.region]
                ))
            return(nodes, links)

        nodes, links = extract(self)

        chart = MermaidDiagram(
            title="Solving",
            nodes=list(nodes.values()),
            links=links
        )

        print(f"Generate solving graph ... ", end="")
        status = generate_mermaid_graph_as_image(f"{chart}", "fig_graph_solving.png")
        print(Keyword.DONE if status else Keyword.ERROR)

        return(chart)


    def plot_solving_flowchart(self):
        """
        Plot solving flowchart connecting solvers and multi-physics loops.
        """
        from python_mermaid.diagram import MermaidDiagram, Node, Link

        def extract(l, loopName: str='Time Loop', info: str=None):
            if (info is None):
                info = loopName

            loopNode = Node(info, shape='trapezoid')
            previousNode = loopNode
            nodes = {loopName: previousNode}
            links = []
            for solver in l:
                if (isinstance(solver, Solver) and solver.region not in nodes):
                    nodes[solver.region] = Node(solver.region)

                    links.append(Link(
                        previousNode, nodes[solver.region]
                    ))
                    previousNode = nodes[solver.region]

                # Create an inner loop
                if (isinstance(solver, MultiPhysicsLoop)):
                    newNodes, newLinks = extract(
                        solver,
                        loopName=solver.region,
                        info=f"`{solver.region}\ntype: **{solver.solver}**\nmaxResidual: {solver.maxResidual}\nmaxIteration: {solver.maxIterations}`"
                    )
                    nodes = nodes | newNodes
                    links += newLinks

                    links.append(Link(
                        previousNode, nodes[solver.region]
                    ))
                    previousNode = list(newNodes.values())[-1]

            # Close the loop
            decisionNode = Node(f'{loopName} break', shape='trapezoid-alt')
            nodes[f'{loopName} break'] = decisionNode
            links.append(Link(previousNode, decisionNode))
            links.append(Link(decisionNode, loopNode))

            return(nodes, links)

        nodes, links = extract(self)

        nodes['start'] = Node("Start", shape='circle')
        nodes['end'] = Node("`End`", shape='double-circle')
        links.append(Link(nodes['start'], nodes['Time Loop']))
        links.append(Link(nodes['Time Loop break'], nodes['end']))

        chart = MermaidDiagram(
            title="Solving Flowchart",
            nodes=list(nodes.values()),
            links=links
        )

        print(f"Generate solving flowchart ... ", end="")
        status = generate_mermaid_graph_as_image(f"{chart}", "fig_graph_solving_flowchart.png")
        print(Keyword.DONE if status else Keyword.ERROR)

        return(chart)
