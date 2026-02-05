from __future__ import annotations

from foamForNuclear.common import *
from foamForNuclear.numerics import fvSolutionFMI
from foamForNuclear.openfoamFile import OpenFOAMFile

from foamForNuclear.checkvalue import check_type

import attrs as attr
from attrs import define, field, validators as v
from foamForNuclear._attrs_tools import auto_type_validator
from foamForNuclear.solvers._solvers import Solvers


_START_FROM_TYPES = {'startTime', 'firstTime', 'latestTime'}
_STOP_AT_TYPES = {'endTime', 'noWriteNow', 'writeNow', 'nextWrite'}
_WRITE_CONTROL_TYPES = {
    'none', 'timeStep', 'runTime', 'adjustable', 'adjustableRunTime',
    'writeTime', 'clockTime', 'cpuTime'
}
_WRITE_FORMAT_TYPES = {'ascii', 'binary'}
_TIME_FORMAT_TYPES = {'general', 'fixed', 'scientific'}


class SimulationParametersFMI(OpenFOAMFile):
    def __init__(
            self,
            host: str="\"127.0.0.1\"",
            port: int=8000,
            STDGRAD: str="pointCellsLeastSquares",
            STDLAP: str="Gauss linear corrected"
        ):
        super().__init__("simulationParameters", "system")

        self.host = host
        self.port = port
        self.STDGRAD = STDGRAD
        self.STDLAP = STDLAP

    @property
    def host(self):
        return self._host

    @host.setter
    def host(self, host) -> None:
        check_type("host", host, str)
        self._host = host

    @property
    def port(self):
        return self._port

    @port.setter
    def port(self, port) -> None:
        check_type("port", port, int)
        self._port = port


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""

        text += addParameter("host", self.host, isAddExtraLine=True)
        text += addParameter("port", self.port, isAddExtraLine=True)
        text += addParameter("STDGRAD", self.STDGRAD, isAddExtraLine=True)
        text += addParameter("STDLAP", self.STDLAP, isAddExtraLine=True)

        return(text)


@define(
    slots=True,
    on_setattr=[attr.setters.convert, attr.setters.validate],
    field_transformer=auto_type_validator,
    repr=False,
    kw_only=True
)
class ControlDict(OpenFOAMFile):
    """
    The controlDict file contains the main parameter settings.

    Parameters
    ----------
    application: str
        Name of the solver application (e.g. `GeN-Foam`, `OFFBEAT`, ...)
    endTime : float
        End time
    deltaT : float
        Time step
    startFrom : float {"startTime", "firstTime", "latestTime"}
        Start from keyword.
    stopAt : str {"endTime", "noWriteNow", "writeNow", "nextWrite"}
        Stop at keyword.
    startTime : float
        Start time.
    writeControl : str {"none", "timeStep", "runTime", "adjustable", "adjustableRunTime", "clockTime", "cpuTime"}
        Write control keyword.
    writeInterval : float
        Write interval.
    purgeWrite : float
        Purge write.
    writeFormat : str {"ascii", "binary"}
        Write format keyword.
    writePrecision : int
        Write precision.
    writeCompression : bool
        Write compression flag.
    timeFormat : str {"general", "fixed", "scientific"}
        Time format keyword.
    timePrecision : int
        Time precision.
    runTimeModifiable : bool
        Run time modifiable flag (default `False`).
    libs: OpenFOAMList[str]
        List of libraries to be added at runtime.
    fvSolutionFMI: fvSolutionFMI
    simulationParametersFMI: SimulationParametersFMI
    """

    application: str = ''
    # time step
    startFrom: str = field(
        default="latestTime",
        validator=v.and_(v.instance_of(str), v.in_(_START_FROM_TYPES))
    )
    startTime: float | int = 0.0
    stopAt: str = field(
        default="endTime",
        validator=v.and_(v.instance_of(str), v.in_(_STOP_AT_TYPES))
    )
    endTime: float | int = 1.0
    deltaT: float | int = field(
        default=1.0,
        validator=v.and_(v.instance_of((float, int)), v.ge(0.0))
    )

    # write options
    writeControl: str = field(
        default="timeStep",
        validator=v.and_(v.instance_of(str), v.in_(_WRITE_CONTROL_TYPES))
    )
    writeInterval: float | int = field(
        default=1.0,
        validator=v.and_(v.instance_of((float, int)), v.ge(0.0))
    )
    purgeWrite: int = field(
        default=0,
        validator=v.and_(v.instance_of(int), v.ge(0))
    )
    writeFormat: str = field(
        default="ascii",
        validator=v.and_(v.instance_of(str), v.in_(_WRITE_FORMAT_TYPES))
    )
    writePrecision: int = field(
        default=8,
        validator=v.and_(v.instance_of(int), v.ge(0.0))
    )
    writeCompression: bool = False
    timeFormat: str = field(
        default="general",
        validator=v.and_(v.instance_of(str), v.in_(_TIME_FORMAT_TYPES))
    )
    timePrecision: int = field(
        default=8,
        validator=v.and_(v.instance_of(int), v.ge(0.0))
    )
    runTimeModifiable: bool = False

    libs: OpenFOAMList = field(factory=lambda: OpenFOAMList(str, "libs"))
    solvers: Solvers | None = None
    fvSolutionFMI: fvSolutionFMI = field(factory=fvSolutionFMI)
    simulationParametersFMI: SimulationParametersFMI = field(factory=SimulationParametersFMI)

    def __attrs_post_init__(self):
        super().__init__(name="controlDict", folder="system")

    def __repr__(self):
        text  = underline("Application settings:") + "\n"
        text += f"{tab}Application: {self.application}\n"
        if (self.startFrom == "startTime"):
            text += f"{tab}Start time: {self.startTime} s ({convertTimeUnit(self.startTime)})\n"
        text += f"{tab}End time: {self.endTime} s ({convertTimeUnit(self.endTime)})\n"

        return(text)

    def _export_body(self) -> str:
        text = ""

        if (not self.libs.is_empty):
            text += f"{self.libs!r}\n"

        text += addParameter('application', self.application, isAddExtraLine=True)
        text += addParameter('startFrom', self.startFrom, isAddExtraLine=True)
        text += addParameter('stopAt', self.stopAt, isAddExtraLine=True)
        text += addParameter('startTime', self.startTime, isAddExtraLine=True)
        text += addParameter('endTime', self.endTime, isAddExtraLine=True)
        text += addParameter('deltaT', self.deltaT, isAddExtraLine=True)
        text += addParameter('writeControl', self.writeControl, isAddExtraLine=True)
        text += addParameter('writeInterval', self.writeInterval, isAddExtraLine=True)
        text += addParameter('purgeWrite', self.purgeWrite, isAddExtraLine=True)
        text += addParameter('writeFormat', self.writeFormat, isAddExtraLine=True)
        text += addParameter('writePrecision', self.writePrecision, isAddExtraLine=True)
        text += addParameter('writeCompression', self.writeCompression, isAddExtraLine=True)
        text += addParameter('timeFormat', self.timeFormat, isAddExtraLine=True)
        text += addParameter('timePrecision', self.timePrecision, isAddExtraLine=True)
        text += addParameter('runTimeModifiable', self.runTimeModifiable, isAddExtraLine=True)

        if (self.solvers is not None):
            text += self.solvers.write_to_controlDict()

        return text

    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self) -> str:
        # IMPORTANT: call the virtual body, so subclasses extend it
        return self._export_body()

    def import_from_openfoam(self):
        foamFile = super().import_from_openfoam()
        attributes = self.get_attributes_as_list()
        for attr in attributes:
            if (attr in foamFile.keys()):
                setattr(self, attr, foamFile[attr])


@define(
    slots=True,
    on_setattr=[attr.setters.convert, attr.setters.validate],
    field_transformer=auto_type_validator,
    repr=False,
    kw_only=True
)
class ControlDictFfn(ControlDict):
    """
    The controlDict file for foamForNuclear (i.e. GeN-Foam) simulations that
    contains the main parameter settings.
    """
    adjustTimeStep: bool = False
    maxDeltaT: float | int = field(
        default=1.0,
        validator=v.and_(v.instance_of((float, int)), v.ge(0.0))
    )
    minDeltaT: float | int = field(
        default=None,
        validator=v.optional(v.and_(v.instance_of((float, int)), v.ge(0.0)))
    )
    maxCo: float | int = field(
        default=1.0,
        validator=v.and_(v.instance_of((float, int)), v.ge(0.0))
    )
    maxAlphaCo: float | int = field(
        default=None,
        validator=v.optional(v.and_(v.instance_of((float, int)), v.ge(0.0)))
    )
    maxPowerVariation: float | int = field(
        default=0.01,
        validator=v.and_(v.instance_of((float, int)), v.ge(0.0))
    )
    maxCoTwoPhase: float | int | None = None

    liquidFuel: bool = False
    solveFMI: bool = field(default=None, validator=v.optional(v.instance_of(bool)))
    marginToPhaseChange: bool = field(default=None, validator=v.optional(v.instance_of(bool)))
    writeRestartFields: bool = field(default=None, validator=v.optional(v.instance_of(bool)))
    writeContinuityErrors: bool = field(default=None, validator=v.optional(v.instance_of(bool)))
    includeKineticEnergy: bool = field(default=None, validator=v.optional(v.instance_of(bool)))

    # def __attrs_post_init__(self):
    #     super().__init__()

    def __repr__(self):
        text = super().__repr__()

        if (self.adjustTimeStep):
            text += f"{tab}Adatative time step (max ΔT: {convertTimeUnit(self.maxDeltaT)}, max Courant: {self.maxCo})\n"
        else:
            text += f"{tab}Fix time step: {self.deltaT} s ({convertTimeUnit(self.deltaT)})\n"

        return(text)


    def _export_body(self) -> str:
        text = super()._export_body()

        text += addParameter('adjustTimeStep', self.adjustTimeStep, isAddExtraLine=True)

        if (self.includeKineticEnergy is not None):
            text += addParameter('includeKineticEnergy', self.includeKineticEnergy, isAddExtraLine=True)
        if (self.writeRestartFields is not None):
            text += addParameter('writeRestartFields', self.writeRestartFields, isAddExtraLine=True)
        if (self.writeContinuityErrors is not None):
            text += addParameter('writeContinuityErrors', self.writeContinuityErrors, isAddExtraLine=True)
        text += addParameter('maxDeltaT', self.maxDeltaT, isAddExtraLine=True)
        if (self.minDeltaT is not None):
            text += addParameter('minDeltaT', self.minDeltaT, isAddExtraLine=True)
        text += addParameter('maxCo', self.maxCo, isAddExtraLine=True)
        if (self.maxAlphaCo is not None):
            text += addParameter('maxAlphaCo', self.maxAlphaCo, isAddExtraLine=True)
        if (self.maxCoTwoPhase is not None):
            text += addParameter('maxCoTwoPhase', self.maxCoTwoPhase, isAddExtraLine=True)

        if (self.marginToPhaseChange is not None):
            text += addParameter('marginToPhaseChange', self.marginToPhaseChange, isAddExtraLine=True)
        text += addParameter('maxPowerVariation', self.maxPowerVariation, isAddExtraLine=True)

        if (self.liquidFuel):
            text += addParameter('liquidFuel', self.liquidFuel, isAddExtraLine=True)

        if (self.solveFMI is not None):
            text += addParameter('solveFMI', self.solveFMI, isAddExtraLine=True)

            if (self.solveFMI):
                self.fvSolutionFMI.export_to_openfoam()
                self.simulationParametersFMI.export_to_openfoam()

                with open('Allrun', 'w') as f:
                    allrunText = "#!/bin/sh\n"
                    allrunText += r'cd "${0%/*}" || exit' + "\n"
                    allrunText += ". ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions\n"
                    allrunText += "runApplication dictModifier dictMod.json\n"
                    allrunText += f"rm log.{self.application}\n"
                    allrunText += f"runApplication {self.application}\n"
                    f.write(allrunText)

        return(text)


@define(
    slots=True,
    on_setattr=[attr.setters.convert, attr.setters.validate],
    field_transformer=auto_type_validator,
    repr=False,
    kw_only=True
)
class ControlDictOffbeat(ControlDict):
    """
    The controlDict file for OFFBEAT simulations that contains the main
    parameter settings.
    """
    userTime: str = "seconds"
    adjustableTimeStep: bool = False
    maxDeltaT: float | int = field(
        default=None,
        validator=v.optional(v.and_(v.instance_of((float, int)), v.ge(0.0)))
    )
    minDeltaT: float | int = field(
        default=None,
        validator=v.optional(v.and_(v.instance_of((float, int)), v.ge(0.0)))
    )
    maxRelativeDeltaTIncrease: float | int = field(
        default=None,
        validator=v.optional(v.and_(v.instance_of((float, int)), v.ge(0.0)))
    )
    minRelativeDeltaTDecrease: float | int = field(
        default=None,
        validator=v.optional(v.and_(v.instance_of((float, int)), v.ge(0.0)))
    )
    maxRelativePowerIncrease: float | int = field(
        default=None,
        validator=v.optional(v.and_(v.instance_of((float, int)), v.ge(0.0)))
    )
    maxRelativePowerDecrease: float | int = field(
        default=None,
        validator=v.optional(v.and_(v.instance_of((float, int)), v.ge(0.0)))
    )
    maxBurnupIncrease: float | int = field(
        default=None,
        validator=v.optional(v.and_(v.instance_of((float, int)), v.ge(0.0)))
    )
    maxAverageCreep: float | int = field(
        default=None,
        validator=v.optional(v.and_(v.instance_of((float, int)), v.ge(0.0)))
    )
    maxMaximumCreep: float | int = field(
        default=None,
        validator=v.optional(v.and_(v.instance_of((float, int)), v.ge(0.0)))
    )
    maxFGR: float | int = field(
        default=None,
        validator=v.optional(v.and_(v.instance_of((float, int)), v.ge(0.0)))
    )

    # def __attrs_post_init__(self):
    #     super().__init__()

    def __repr__(self):
        text = super().__repr__()

        if (self.adjustableTimeStep):
            text += f"{tab}Adatative time step (max ΔT: {convertTimeUnit(self.maxDeltaT)}, min ΔT: {convertTimeUnit(self.minDeltaT)})\n"
        else:
            text += f"{tab}Fix time step: {self.deltaT} s ({convertTimeUnit(self.deltaT)})\n"

        return(text)

    def _export_body(self) -> str:
        text = super()._export_body()

        if (self.userTime is not None):
            text += addParameter('userTime', OpenFOAMDict({'type': self.userTime}), isAddExtraLine=False)
        text += addParameter('adjustableTimeStep', self.adjustableTimeStep, isAddExtraLine=True)

        if (self.maxDeltaT is not None):
            text += addParameter('maxDeltaT', self.maxDeltaT, isAddExtraLine=True)
        if (self.minDeltaT is not None):
            text += addParameter('minDeltaT', self.minDeltaT, isAddExtraLine=True)
        if (self.maxRelativeDeltaTIncrease is not None):
            text += addParameter('maxRelativeDeltaTIncrease', self.maxRelativeDeltaTIncrease, isAddExtraLine=True)
        if (self.minRelativeDeltaTDecrease is not None):
            text += addParameter('minRelativeDeltaTDecrease', self.minRelativeDeltaTDecrease, isAddExtraLine=True)
        if (self.maxRelativePowerIncrease is not None):
            text += addParameter('maxRelativePowerIncrease', self.maxRelativePowerIncrease, isAddExtraLine=True)
        if (self.maxRelativePowerDecrease is not None):
            text += addParameter('maxRelativePowerDecrease', self.maxRelativePowerDecrease, isAddExtraLine=True)
        if (self.maxBurnupIncrease is not None):
            text += addParameter('maxBurnupIncrease', self.maxBurnupIncrease, isAddExtraLine=True)
        if (self.maxAverageCreep is not None):
            text += addParameter('maxAverageCreep', self.maxAverageCreep, isAddExtraLine=True)
        if (self.maxMaximumCreep is not None):
            text += addParameter('maxMaximumCreep', self.maxMaximumCreep, isAddExtraLine=True)
        if (self.maxFGR is not None):
            text += addParameter('maxFGR', self.maxFGR, isAddExtraLine=True)

        return text
