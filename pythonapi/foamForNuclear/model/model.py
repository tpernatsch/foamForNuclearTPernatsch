import os

import numpy as np
import matplotlib.pyplot as plt
import time as timePack
import pyvista as pv
import foamlib
import tqdm

from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import *
from foamForNuclear.coupling import Coupling
from foamForNuclear.externalCouplingDict import ExternalCouplingDict
from foamForNuclear.functionObjects import FMUSimulator, FunctionObject
from foamForNuclear.mesh.mesh import Mesh
from foamForNuclear.solver import Solver, Solvers
from foamForNuclear.controlDict import ControlDict
from foamForNuclear.timeFolder import TimeFolder


class Residuals:
    def __init__(self, parameters: list, title: str="", nLastIter: int=0) -> None:
        self.parameters: dict = {}
        for param in parameters:
            self.parameters[param] = []

        self.title: str = title
        self.nLastIter: int = nLastIter
        self.filename: str = ""

    def extractFromFile(
            self,
            filename: str,
            isResidual: bool=True,
            idx: int=-2,
            scale: float=1
        ):
        """
        Extraction of the residuals in the log file

        Parameters
        ----------
        filename : str
            Name of the log file
        isResidual : bool, default True
            Flag to indicate if parameters are residuals
        idx : int, default -2
            Index of the value when line are seprated by spaces
        scale : float, default 1
            Scale the value on the fly
        """

        self.filename = filename

        if (os.path.isfile(filename)):
            with open(filename, 'r') as file:
                print(f"Analyze {filename} ... ", end="")
                for line in file:
                    for param in self.parameters:
                        if (isResidual):
                            if (f"Solving for {param}," in line):
                                self.parameters[param].append(
                                    float(line.split()[7].split(",")[0])
                                )
                            elif (f"Residual {param}" in line):
                                self.parameters[param].append(
                                    float(line.split()[3].split(",")[0])
                                )
                        else:
                            if (param in line):
                                self.parameters[param].append(
                                    float(line.split()[idx]) * scale
                                )
                print(Keyword.DONE)
        else:
            print(f"{filename} does not exist")

    def addTitle(self, ax, suffix: str=""):
        title = self.title
        if (self.title == ""):
            title = self.filename.split(".")[-1]
        if (suffix != ""):
            title = f"{title} {suffix}"

        ax.set_title(title)

    def plot(
            self,
            ax,
            axLast=None,
            ls: str="-",
            prefix: str="",
            modulo: int=1,
            idx: int=0
        ) -> None:
        """ Plot all the residuals """

        self.addTitle(ax)
        if (axLast != None):
            self.addTitle(axLast, suffix=f"Last {self.nLastIter}")

        for param in self.parameters:
            if (len(self.parameters[param]) > 0):
                values = self.parameters[param][idx::modulo]
                ax.plot(
                    np.arange(1, len(values)+1, 1),
                    values,
                    ls=ls,
                    label=prefix+param.split("=")[0]
                )
                if (axLast != None):
                    axLast.plot(
                        np.arange(1, len(values[-self.nLastIter:])+1, 1),
                        values[-self.nLastIter:],
                        ls=ls,
                        label=prefix+param.split("=")[0]
                    )


class Model:
    """
    Base class to define an OpenFOAM-based case.

    Parameters
    ----------
    settings : ControlDict
        Main parameters of the simulation.
    solvers : Solvers
        Collection of solvers to be used.
    coupling : Coupling
        Coupling object used to couple multiple regions.
    timeFolders : list[TimeFolder]
        List of time folders.
    externalCouplingDict : ExternalCouplingDict
        External coupling dict object for FMU coupling.
    caseFolder : str, optional
        Case folder (default `"./"`).

    Attributes
    ----------
    settings : ControlDict
        Main parameters of the simulation (default `ControlDict()`)
    solvers : Solvers
        Collection of solvers to be used (default `Solvers()`)
    coupling : Coupling
        Coupling object used to couple multiple regions
    timeFolders : list[TimeFolder]
        List of time folders
    externalCouplingDict : ExternalCouplingDict
        External coupling dict object for FMU coupling
    caseFolder : str
        Case folder (default `"./"`).
    """

    def __init__(
            self,
            settings: ControlDict=None,
            solvers: Solvers=None,
            timeFolders: list[TimeFolder]=[],
            coupling: Coupling=None,
            externalCouplingDict: ExternalCouplingDict=None,
            caseFolder: str="./"
        ):
        self.settings: ControlDict = ControlDict() if settings is None else settings
        self.solvers: Solvers = Solvers() if solvers is None else solvers
        self.coupling: Coupling = coupling
        self.timeFolders: list[TimeFolder] = timeFolders
        self.externalCouplingDict: ExternalCouplingDict = externalCouplingDict
        self.caseFolder: str = caseFolder


    def __repr__(self):
        isParallel = self.is_parallel
        maxSubdomains = self.get_number_processors

        text = f"{self.settings}"
        text += underline('Solvers:') + "\n"
        if (len(self.solvers) > 0):
            for solver in self.solvers:
                text += f"{solver.__repr__(depth = 1)}"

                if (isParallel and solver.decomposeParDict.numberOfSubdomains != maxSubdomains):
                    text += f"{tab}  {Keyword.WARNING}: Number of subdomains is different than the max subdomains ({maxSubdomains})\n"

                # Check fields are present
                if (solver.solver == "twoPhase"):
                    continue

                requiredFields = solver.get_required_fields()
                createdFields = [field.name for field in self.timeFolders[0]]
                for field in requiredFields:
                    if (field not in createdFields):
                        text += f"{tab}  {Keyword.WARNING}: Field '{field}' has NOT been declared in the time folder.\n"

        else:
            text += f"{tab}None\n"

        # Add coupling info if present
        if (self.coupling is not None):
            text += f"{self.coupling}"

            if (len(self.solvers) != len(self.coupling.get_unique_regions(self.coupling))):
                text += f"{Keyword.WARNING}: Mismatch number of regions in 'Solvers' and 'Coupling' objects\n"

        # Add FMU simulator info is present
        FMUSimulator = self.get_FMU_simulator()
        if (FMUSimulator is not None):
            text += underline("FMU simulation:") + "\n"
            text += f"{tab}Found {FMUSimulator.name} functionObject\n"

            text += FMUSimulator.get_coupling_mapping_as_text(depth=1)

        return(text)


    @property
    def coupling(self):
        return self._coupling

    @coupling.setter
    def coupling(self, coupling) -> None:
        check_type("coupling", coupling, Coupling, none_ok=True)
        self._coupling = coupling


    @property
    def externalCouplingDict(self):
        return self._externalCouplingDict

    @externalCouplingDict.setter
    def externalCouplingDict(self, externalCouplingDict) -> None:
        check_type("externalCouplingDict", externalCouplingDict, ExternalCouplingDict, none_ok=True)
        self._externalCouplingDict = externalCouplingDict


    @property
    def caseFolder(self):
        return self._caseFolder

    @caseFolder.setter
    def caseFolder(self, caseFolder) -> None:
        check_type("caseFolder", caseFolder, str, none_ok=True)
        self._caseFolder = caseFolder


    def add_solver(self, solver: Solver):
        if solver is not None:
            check_type('solver', solver, Solver)
            self.solvers.append(solver)


    def add_time_folder(self, timeFolder: TimeFolder):
        if timeFolder is not None:
            check_type('timeFolder', timeFolder, TimeFolder)
            self.timeFolders.append(timeFolder)


    def add_function_object(self, functionObject: FunctionObject):
        check_type("functionObject", functionObject, FunctionObject)
        self.settings.add_function_object(functionObject)


    @property
    def is_parallel(self) -> bool:
        return(self.solvers.is_one_solver_parallel())


    @property
    def get_number_processors(self) -> int:
        return(self.solvers.get_max_subdomains())


    def get_FMU_simulator(self) -> FunctionObject:
        """
        Return the first FMUSimulator functionObject.
        """
        for functionObject in self.settings.functionObjects:
            if (isinstance(functionObject, FMUSimulator)):
                return(functionObject)
        return(None)


    def export_to_openfoam(self):
        # Create case folder
        for i, folder in enumerate(self.caseFolder.split("/")):
            folder = '/'.join(self.caseFolder.split('/')[:i+1])
            if (not os.path.exists(folder)):
                os.mkdir(folder)

        # Change directory to avoid manipulating all the file genration
        cwd = os.getcwd()
        os.chdir(self.caseFolder)

        if (not os.path.exists("constant")):
            os.mkdir("constant")
        if (not os.path.exists("system")):
            os.mkdir("system")

        if (len(self.solvers) == 1):
            self.coupling = Coupling(self.solvers)

        if (self.coupling is not None and self.settings.application == "GeN-Foam"):
            self.coupling.export_to_openfoam()
            self.settings.coupling = self.coupling
            self.settings.solvers = self.solvers

        # External coupling dict for FMI
        if (self.externalCouplingDict is not None):
            self.externalCouplingDict.export_to_openfoam()
        elif (self.settings.solveFMI is not None and self.settings.solveFMI):
            self.externalCouplingDict = ExternalCouplingDict()
            self.externalCouplingDict.export_to_openfoam()

        if (len(self.timeFolders) == 1):
            self.solvers.set_time_folder(self.timeFolders[0])

        self.settings.export_to_openfoam()
        self.solvers.export_to_openfoam()
        for timeFolder in self.timeFolders:
            timeFolder.export_to_openfoam()

        # Return to the current directory
        os.chdir(cwd)


    def keff(self, time: float=None) -> float:
        """
        Return the keff at the last time step or at the time step specified by
        the parameter `time`.

        Parameters
        ----------
        time : float
            Time step at which to extract the keff (default last time step,
            specified in the `settings`)
        """
        if (time is None):
            time = self.settings.endTime

        f = foamlib.FoamFile(f"{self.caseFolder}/{time:g}/uniform/reactorState")

        return(f['keff'])


    def get_parameters_from_log(
            self,
            parameters: list[tuple],
            writeInterval: float=None
        ) -> dict[str, list]:
        """
        Parameters
        ----------
        parameters : list[tuple]
            (paramName, paramName in log, position in line)
        writeInterval : float
            Write every `writeInterval` in the results list
        """
        filename = f"{self.caseFolder}/log.{self.settings.application}"

        timeList = []
        values = {paramName: [] for paramName, _, _ in parameters}
        valuesTemp = {paramName: None for paramName, _, _ in parameters}

        with open(filename, 'r') as f:
            print(f"Read {filename} ...")
            previousTime = 0
            for line in f:
                for paramName, paramNameInLog, valuePosition in parameters:
                    if (paramNameInLog in line):
                        valuesTemp[paramName] = float(line.split()[valuePosition])

                if ("Time =" == line[:6] and all([item is not None for item in valuesTemp.values()])):
                    newTime = float(line.split()[2])

                    # Skip writing
                    if (writeInterval is not None and previousTime+writeInterval > newTime):
                        continue

                    timeList.append(newTime)
                    for key, value in valuesTemp.items():
                        values[key].append(value)

                    previousTime = newTime

        return({'time': timeList} | values)

    def get_keff_from_log(self, writeInterval: float=None):
        return(self.get_parameters_from_log(
            parameters=[
                ("keff", "keff_ = ", 2)
            ],
            writeInterval=writeInterval
        ))

    def get_parameters_from_point_kinetics(
            self,
            writeInterval: float=None
        ) -> dict[str, list]:
        """
        Extract point-kinetics parameters evolution over time.
        """
        return(self.get_parameters_from_log(
            parameters=[
                ("totalPower", "totalPower = ", 2),
                ("fissionPower", "-> fission = ", 3),
                ("decayPower", "-> decay   = ", 3),
                ("totalReactivity", "totalReactivity", 2),
                ("extReactivity", "-> extReactivity", 3),
                ("dopplerReactivity", "-> Doppler (fast) =", 4),
                ("TFuelReactivity", "-> TFuel", 3),
                ("TCladReactivity", "-> TClad", 3),
                ("TCoolReactivity", "-> TCool", 3),
                ("rhoCoolReactivity", "-> rhoCool", 3),
                ("TStructReactivity", "-> TStruct", 3),
                ("TStructMechReactivity", "-> TStructMech", 3),
                ("drivelineReactivity", "-> driveline", 3),
                ("BoronReactivity", "-> Boron", 3),
                ("TFuel", "TFuel = ", 2),
                ("TClad", "TClad = ", 2),
                ("TCool", "TCool = ", 2),
                ("rhoCool", "rhoCool = ", 2),
                ("TStruct", "TStruct = ", 2),
                ("TStructMech", "TStructMech = ", 2),
                ("driveline", "driveline (T, exp) = ", 4),
            ],
            writeInterval=writeInterval
        ))


    def plot_residuals(
            self,
            parameters: list[str],
            title: str='',
            nLastIter: int=100
        ):
        t1 = timePack.time()

        # Residuals to plot
        ResidualsGFneutro = Residuals(
            parameters=parameters,
            title=title,
            nLastIter=nLastIter
        )

        # Extract from file
        ResidualsGFneutro.extractFromFile(filename=f"{self.caseFolder}/log.{self.settings.application}")


        # Create figure
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
        # fig.subplots_adjust(hspace=0.5)
        ax, axLast = axes.flatten()

        # Plot the resisuals
        ResidualsGFneutro.plot(ax=ax, axLast=axLast)

        # Options
        for ax in axes.flatten():
            ax.set_xlabel("Iterations")
            ax.set_ylabel("Residual")
            ax.set_yscale('log')
            ax.grid()
            ax.legend()

        # Save figures
        fig.tight_layout()
        fig.savefig(f"fig_results_residuals_{'-'.join(parameters)}.png")

        t2 = timePack.time()

        print(f"Processing time for plotting residuals = {t2-t1:.6g} s")


    def get_time_steps(self):
        """
        Return the list of available time step folders.
        """
        reader = pv.POpenFOAMReader(f"./{self.caseFolder}/foam.foam")
        return([t if '.' in f"{t:g}" else int(t) for t in reader.time_values])


    def get_data_range(
            self,
            region: str | Mesh,
            time: float,
            fieldName: str,
            removeZeros: bool=False
        ) -> tuple[float, float]:
        """
        Return the minimum and the maximum of a field.

        Parameters
        ----------
        region : str | Mesh
            Name of the region
        time : float
            Time folder
        fieldName : str
            Name of the field
        removeZeros : bool
            Avoid return 0 as mimimum of a field, useful for field such as
            `T.fuelAvForNeutronics` (default `False`).

        Return
        ------
        Return the minimum and the maximum of a field.
        """
        regionName = region
        if (isinstance(region, Mesh)):
            regionName = region.region

        reader = pv.POpenFOAMReader(f"./{self.caseFolder}/foam.foam")

        reader.set_active_time_value(time)
        reader.cell_to_point_creation = False

        if (regionName != ""):
            mesh = reader.read()[regionName]
        else:
            mesh = reader.read()

        internalMesh = mesh['internalMesh']

        field = internalMesh.cell_data[fieldName]
        if (removeZeros):
            field = [value for value in field if value > 0]

        return(np.min(field), np.max(field))


    def sample_over_line(
            self,
            region: str | Mesh,
            time: float,
            fieldName: str,
            point1: tuple,
            point2: tuple
        ):
        """
        Parameters
        ----------
        region : str | Mesh
            Name of the region
        time : float
            Time folder
        fieldName : str
            Name of the field
        point1 : tuple
            First point (e.g `(0, 0, 0)`)
        point2 : tuple
            Second point (e.g `(0, 0, 1)`)
        """
        regionName = region
        if (isinstance(region, Mesh)):
            regionName = region.region

        reader = pv.POpenFOAMReader(f"./{self.caseFolder}/foam.foam")

        reader.set_active_time_value(time)
        reader.cell_to_point_creation = False

        if (regionName != ""):
            mesh = reader.read()[regionName]
        else:
            mesh = reader.read()

        internalMesh = mesh['internalMesh']

        sample = internalMesh.sample_over_line(point1, point2, tolerance=1e-6)

        return(sample.points, sample[fieldName])


    def plot_slice(
            self,
            region: str | Mesh,
            time: float=0,
            fieldName: str=None,
            normal: str='y',
            offset: tuple=(0, 0, 0),
            scalingVector: list=[1, 1, 1],
            cmap: str='inferno',
            lighting: bool=True,
            show_edges: bool=False,
            unit: str='-',
            isVerticalLegend: bool=True,
            limits: list[float]=None
        ):
        """
        Parameters
        ----------
        region : str | Mesh | list[str] | list[Mesh]
            Name of the region or list of name of region
        time : float
            Time step (default `0`)
        fieldName : str
            Field name to color the mesh (default `None`)
        normal : str
            Change camera view along the specified axis (`x`, `y`, or `z`)
        offset : tuple[float]
            Offset vector (default `(0, 0, 0)`)
        scalingVector : list[float]
            Scaling vector (default `[1, 1, 1]`).
        cmap : str
            Colormap name based on the values of the `fieldName`
        lighting : bool
            Use light on the mesh (default `True`)
        show_edges : bool, default False
            Show mesh edges
        unit : str, default '-'
            Add unit to the color
        isVerticalLegend : bool, default True
            If true, the color bar is vertical, else it is horizontal
        limits : list[float]
            Limit the colormap range (default `None`).
        """
        self.plot_mesh(
            region=region,
            time=time,
            fieldName=fieldName,
            cmap=cmap,
            normal=normal,
            offset=offset,
            scalingVector=scalingVector,
            lighting=lighting,
            show_edges=show_edges,
            unit=unit,
            isVerticalLegend=isVerticalLegend,
            limits=limits,
            isSlice=True
        )


    def plot_mesh(
            self,
            region: str | Mesh,
            time: float=0,
            fieldName: str=None,
            cmap: str=None,
            normal: str=None,
            offset: tuple=(0, 0, 0),
            scalingVector: list=[1, 1, 1],
            lighting: bool=True,
            show_edges: bool=False,
            unit: str='-',
            isVerticalLegend: bool=True,
            limits: list[float]=None,
            isSlice: bool=False
        ):
        """
        Plot the mesh. Apply grey color if no field is provided.

        Parameters
        ----------
        region : str | Mesh | list[str] | list[Mesh]
            Name of the region or list of name of region
        time : float
            Time step
        fieldName : str
            Field name to color the mesh
        cmap : str
            Colormap name based on the values of the `fieldName`
        normal : str
            Change camera view along the specified axis (`x`, `y`, or `z`)
        offset : tuple[float]
            Offset vector (default `(0, 0, 0)`).
        scalingVector : list[float]
            Scaling vector (default `[1, 1, 1]`).
        lighting : bool
            Use light on the mesh (default True)
        show_edges : bool, default False
            Show mesh edges
        unit : str, default '-'
            Add unit to the color
        isVerticalLegend : bool, default True
            If true, the color bar is vertical, else it is horizontal.
        limits : list[float]
            Limit the colormap range (default `None`).
        """
        if (isinstance(region, str)):
            regionNames = [region]
        elif (isinstance(region, Mesh)):
            regionNames = [region.region]
        elif (isinstance(region, list)):
            regionNames = [obj.region if isinstance(obj, Mesh) else obj for obj in region]
        else:
            msg = f"'region' must be of type str, Mesh, list[str] or list[Mesh]. Found {type(region).__name__}"
            raise ValueError(msg)

        reader = pv.POpenFOAMReader(f"./{self.caseFolder}/foam.foam")

        reader.set_active_time_value(time)
        reader.cell_to_point_creation = False

        plotter = pv.Plotter(off_screen=True)

        for regionName in regionNames:

            if (regionName != ""):
                mesh = reader.read()[regionName]
            else:
                mesh = reader.read()

            internalMesh = mesh['internalMesh']
            internalMesh = internalMesh.scale(scalingVector)
            if (isSlice):
                slice_mesh = internalMesh.slice(normal=normal)
                internalMesh = slice_mesh.translate(offset, inplace=True)

            # Recompute limits
            if (fieldName is not None and isinstance(limits, (list, tuple))):
                if (limits[0] is None):
                    limits = [np.min(internalMesh.cell_data[fieldName]), limits[1]]
                if (limits[1] is None):
                    limits = [limits[0], np.max(internalMesh.cell_data[fieldName])]

            plotter.add_mesh(
                internalMesh,
                color="grey" if fieldName is None else None,
                scalars=fieldName,
                cmap=cmap,
                lighting=lighting,
                show_edges=show_edges,
                show_scalar_bar=fieldName is not None,
                scalar_bar_args={
                    'title': f'{fieldName} [{unit}]',
                    'vertical': isVerticalLegend,
                    'position_x': 0.05,
                    # 'position_y': 0.05
                },
                clim=limits
            )

        plotter.add_text(
            f"{time} s",
            position='upper_left',
        )
        if (isSlice):
            plotter.add_text(
                f"Offset = {offset} m",
                position='upper_right',
            )
        plotter.show_bounds(location='outer')

        if (isSlice):
            normal = normal.lower()

        if (normal == 'x'):
            plotter.view_yz()
            camera_position = 'yz'
        elif (normal == 'y'):
            plotter.view_xz()
            camera_position = 'xz'
        elif (normal == 'z'):
            plotter.view_xy()
            camera_position = 'xy'

        plotter.enable_parallel_projection()

        ext = ""
        figType = 'mesh'
        if (isSlice):
            figType = 'slice'
        if (normal is not None):
            ext += "_" + camera_position
        if (fieldName is not None):
            ext = "_" + fieldName + ext
            figType = 'results'

        plotter.screenshot(f"fig_{figType}_{'_'.join(regionNames)}_{time}{ext}.png")


    def plot_boundary(
            self,
            region: str | Mesh,
            boundaryName: str,
            time: float=0,
            normal: str=None,
            show_edges: bool=True,
        ):
        """
        Parameters
        ----------
        normal : str
            Change camera view along the specified axis (`x`, `y`, or `z`)
        """
        regionName = region
        if (isinstance(region, Mesh)):
            regionName = region.region

        if (normal is not None):
            normal = normal.lower()

        reader = pv.POpenFOAMReader(f"./{self.caseFolder}/foam.foam")

        reader.set_active_time_value(time)
        reader.cell_to_point_creation = False

        if (regionName != ""):
            mesh = reader.read()[regionName]
        else:
            mesh = reader.read()

        boundaries = mesh['boundary']

        plotter = pv.Plotter(off_screen=True)
        plotter.add_mesh(
            boundaries[boundaryName],
            color="grey",
            show_edges=show_edges
        )
        plotter.show_bounds(location='outer')

        if (normal == 'x'):
            plotter.view_yz()
            camera_position = 'yz'
        elif (normal == 'y'):
            plotter.view_xz()
            camera_position = 'xz'
        elif (normal == 'z'):
            plotter.view_xy()
            camera_position = 'xy'

        plotter.enable_parallel_projection()

        ext = ""
        if (normal is not None):
            ext += "_" + camera_position
        plotter.screenshot(f"fig_boundary_{regionName}_{boundaryName}{ext}.png")


    def plot_animation(
            self,
            region: str | Mesh,
            fieldName: str=None,
            cmap: str=None,
            normal: str=None,
            lighting: bool=True,
            show_edges: bool=False,
            unit: str='-',
            isVerticalLegend: bool=True,
            threshold: list=None,
            thresholdFieldName: str=None,
            fps: int=None,
            offset: tuple=(0, 0, 0),
            limits: list[float]=None
        ):
        """
        Plot animation over all time steps. A uniform time stepping is
        recommanded for a nice render.

        Parameters
        ----------
        region : str | Mesh
            Name of the region
        fieldName : str
            Field name to color the mesh
        cmap : str
            Colormap name based on the values of the `fieldName`
        normal : str
            Change camera view along the specified axis (`x`, `y`, or `z`)
        lighting : bool
            Use light on the mesh (default `True`)
        show_edges : bool, default False
            Show mesh edges
        unit : str, default '-'
            Add unit to the color
        isVerticalLegend : bool, default True
            If true, the color bar is vertical, else it is horizontal
        threshold : list, default None
            If provided, apply a threshold onto the mesh (e.g `[0.9, 1.1]`)
        thresholdFieldName : str, default None
            Name of the threshold field, can be different from `fieldName`
        fps : int, default None
            Number of frame per second in the gif. If default to `None`, the fps
            will be computed to as `1/(timeStep[1] - timeStep[0])`
        limits : list[float]
            Limit the colormap range (default `None`).
        """
        regionName = region
        if (isinstance(region, Mesh)):
            regionName = region.region

        if (normal is not None):
            normal = normal.lower()

        isUseThreshold = threshold is not None

        reader = pv.POpenFOAMReader(f"./{self.caseFolder}/foam.foam")
        reader.cell_to_point_creation = False

        # Create a plotter object and set the scalars to the Z height
        plotter = pv.Plotter(notebook=False, off_screen=True)


        def render(time, limits):
            plotter.clear_actors()

            reader.set_active_time_value(time)

            if (regionName != ""):
                mesh = reader.read()[regionName]
            else:
                mesh = reader.read()


            internalMesh = mesh['internalMesh']

            if (normal is None):
                meshToRender = internalMesh
            else:
                slice_mesh = internalMesh.slice(normal=normal)
                slice_internal_mesh = slice_mesh.translate(offset, inplace=True)
                meshToRender = slice_internal_mesh

            # Recompute limits
            if (fieldName is not None and isinstance(limits, (list, tuple))):
                if (limits[0] is None):
                    limits = [np.min(meshToRender.cell_data[fieldName]), limits[1]]
                if (limits[1] is None):
                    limits = [limits[0], np.max(meshToRender.cell_data[fieldName])]

            plotter.add_mesh(
                meshToRender,
                color='white' if isUseThreshold else None,
                scalars=None if isUseThreshold else fieldName,
                cmap=None if isUseThreshold else cmap,
                lighting=lighting,
                opacity=0.2 if isUseThreshold else 1,
                show_edges=show_edges,
                show_scalar_bar=not isUseThreshold,
                scalar_bar_args={
                    'title': f'{fieldName} [{unit}]',
                    'vertical': isVerticalLegend,
                    'position_x': 0.05,
                    # 'position_y': 0.05
                },
                clim=limits
            )
            plotter.add_text(
                f"{time} s",
                position='upper_left',
            )

            if (threshold is not None):
                plotter.add_mesh(
                    meshToRender.threshold(threshold, scalars=thresholdFieldName),
                    scalars=fieldName,
                    cmap=cmap,
                    lighting=lighting,
                    scalar_bar_args={
                        'title': f'{fieldName} [{unit}]',
                        'vertical': isVerticalLegend,
                        'position_x': 0.05,
                        # 'position_y': 0.05
                    }
                )

            if (normal == 'x'):
                plotter.view_yz()
            elif (normal == 'y'):
                plotter.view_xz()
            elif (normal == 'z'):
                plotter.view_xy()

            plotter.enable_parallel_projection()
            plotter.show_bounds(location='outer')


        if (fps is None):
            fps = 1/(reader.time_values[1] - reader.time_values[0])

        ext = ""
        if (normal is not None):
            if (normal == 'x'):
                camera_position = 'yz'
            elif (normal == 'y'):
                camera_position = 'xz'
            elif (normal == 'z'):
                camera_position = 'xy'
            ext += "_" + camera_position

        # Open a gif
        plotter.open_gif(f'fig_results_{regionName}_{fieldName}{ext}.gif', fps=fps)

        for time in tqdm.tqdm(reader.time_values, desc=f"Render {fieldName}{ext}", unit='frame'):
            try:
                render(time, limits=limits)
                # Write a frame. This triggers a render.
                plotter.write_frame()
            except:
                pass

        # Closes and finalizes movie
        plotter.close()
