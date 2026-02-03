from foamForNuclear.common import *
from foamForNuclear.coupling import Coupling
from foamForNuclear.functionObjects import FMUSimulator, FunctionObject, FunctionObjects
from foamForNuclear.fvSolution import fvSolutionFMI
from foamForNuclear.openfoamFile import OpenFOAMFile
from foamForNuclear.solver import Solvers

from .checkvalue import check_positive, check_type, check_value


_START_FROM_TYPES = {'startTime', 'firstTime', 'latestTime'}
_STOP_AT_TYPES = {'endTime', 'noWriteNow', 'writeNow', 'nextWrite'}
_WRITE_CONTROL_TYPES = {
    'none', 'timeStep', 'runTime', 'adjustable', 'adjustableRunTime',
    'writeTime', 'clockTime', 'cpuTime'
}
_WRITE_FORMAT_TYPES = {'ascii', 'binary'}
_TIME_FORMAT_TYPES = {'general', 'fixed', 'scientific'}


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
    adjustTimeStep : bool
        Adjust time step flag (default `False`).
    solveFMI : bool
        Flag to indicate to solve considering the solver as an FMU. Must be set
        to `False` if the solver is managing an FMU using `OMSimulatorContainer`
        or `PyFMIContainer`
    """

    def __init__(
            self,
            application: str='',
            endTime: float=1,
            deltaT: float=1,
            startFrom: str='latestTime',
            stopAt: str='endTime',
            startTime: float=0.0,
            writeControl: str='timeStep',
            writeInterval: float=1.0,
            purgeWrite: float=0.0,
            writeFormat: str='ascii',
            writePrecision: int=8,
            writeCompression: bool=False,
            timeFormat: str='general',
            timePrecision: int=8,
            runTimeModifiable: bool=False,
            adjustTimeStep: bool=False,
            maxDeltaT: float=1,
            maxCo: float=1,
            maxAlphaCo: float=None,
            maxPowerVariation: float=0.01,
            maxRelativeDeltaTIncrease: float=None,
            minRelativeDeltaTDecrease: float=None,
            maxRelativePowerIncrease: float=None,
            maxRelativePowerDecrease: float=None,
            maxBurnupIncrease: float=None,
            maxAverageCreep: float=None,
            maxMaximumCreep: float=None,
            maxFGR: float=None,
            solveFMI: bool=None,
            minDeltaT: float=None,
            maxCoTwoPhase: float=None,
            marginToPhaseChange: float=None,
            writeRestartFields: bool=None,
            writeContinuityErrors: bool=None,
            includeKineticEnergy: bool=None,
            libs: OpenFOAMList=OpenFOAMList(str, "libs"),
            functionObjects: FunctionObjects=FunctionObjects(),
            **kwargs
        ):
        super().__init__(name="controlDict", folder="system", **kwargs)

        self.application = application
        self.startFrom = startFrom
        self.stopAt = stopAt
        self.startTime = startTime
        self.endTime = endTime
        self.deltaT = deltaT
        self.writeControl = writeControl
        self.writeInterval = writeInterval
        self.purgeWrite = purgeWrite
        self.writeFormat = writeFormat
        self.writePrecision = writePrecision
        self.writeCompression = writeCompression
        self.timeFormat = timeFormat
        self.timePrecision = timePrecision
        self.runTimeModifiable = runTimeModifiable
        self.adjustTimeStep = adjustTimeStep
        self.maxDeltaT = maxDeltaT
        self.maxCo = maxCo
        self.maxAlphaCo = maxAlphaCo
        self.maxPowerVariation = maxPowerVariation
        self.maxRelativeDeltaTIncrease = maxRelativeDeltaTIncrease
        self.minRelativeDeltaTDecrease = minRelativeDeltaTDecrease
        self.maxRelativePowerIncrease = maxRelativePowerIncrease
        self.maxRelativePowerDecrease = maxRelativePowerDecrease
        self.maxBurnupIncrease = maxBurnupIncrease
        self.maxAverageCreep = maxAverageCreep
        self.maxMaximumCreep = maxMaximumCreep
        self.maxFGR = maxFGR
        self.solveFMI = solveFMI
        self.minDeltaT = minDeltaT
        self.maxCoTwoPhase = maxCoTwoPhase
        self.marginToPhaseChange = marginToPhaseChange
        self.writeRestartFields = writeRestartFields
        self.writeContinuityErrors = writeContinuityErrors
        self.includeKineticEnergy = includeKineticEnergy

        self.solvers = None
        self.coupling = None
        self.libs = libs
        self.functionObjects = functionObjects

        self.fvSolutionFMI = fvSolutionFMI()
        self.simulationParametersFMI = SimulationParametersFMI()


    def __repr__(self):
        text  = underline("Application settings:") + "\n"
        text += f"{tab}Application: {self.application}\n"
        if (self.startFrom == "startTime"):
            text += f"{tab}Start time: {self.startTime} s ({convertTimeUnit(self.startTime)})\n"
        text += f"{tab}End time: {self.endTime} s ({convertTimeUnit(self.endTime)})\n"

        if (self.adjustTimeStep):
            text += f"{tab}Adatative time step (max ΔT: {convertTimeUnit(self.maxDeltaT)}, max Courant: {self.maxCo})\n"
        else:
            text += f"{tab}Fix time step: {self.deltaT} s ({convertTimeUnit(self.deltaT)})\n"

        return(text)


    @property
    def application(self):
        return self._application

    @application.setter
    def application(self, application):
        check_type("application", application, str)
        self._application = application

    @property
    def startFrom(self):
        return self._startFrom

    @startFrom.setter
    def startFrom(self, startFrom):
        check_type('startFrom', startFrom, str)
        check_value('startFrom', startFrom, _START_FROM_TYPES)
        self._startFrom = startFrom

    @property
    def stopAt(self):
        return self._stopAt

    @stopAt.setter
    def stopAt(self, stopAt):
        check_type('stopAt', stopAt, str)
        check_value('stopAt', stopAt, _STOP_AT_TYPES)
        self._stopAt = stopAt

    @property
    def startTime(self):
        return self._startTime

    @startTime.setter
    def startTime(self, startTime) -> None:
        check_type("startTime", startTime, (float, int))
        self._startTime = startTime

    @property
    def endTime(self):
        return self._endTime

    @endTime.setter
    def endTime(self, endTime) -> None:
        check_type("endTime", endTime, (float, int))
        self._endTime = endTime

    @property
    def deltaT(self):
        return self._deltaT

    @deltaT.setter
    def deltaT(self, deltaT) -> None:
        check_type("deltaT", deltaT, (float, int))
        check_positive("deltaT", deltaT, is_strict=True)
        self._deltaT = deltaT

    @property
    def writeControl(self):
        return self._writeControl

    @writeControl.setter
    def writeControl(self, writeControl):
        check_type('writeControl', writeControl, str)
        check_value('writeControl', writeControl, _WRITE_CONTROL_TYPES)
        self._writeControl = writeControl

    @property
    def writeInterval(self):
        return self._writeInterval

    @writeInterval.setter
    def writeInterval(self, writeInterval) -> None:
        check_type("writeInterval", writeInterval, (float, int))
        check_positive("writeInterval", writeInterval, is_strict=True)
        self._writeInterval = writeInterval

    @property
    def purgeWrite(self):
        return self._purgeWrite

    @purgeWrite.setter
    def purgeWrite(self, purgeWrite) -> None:
        check_type("purgeWrite", purgeWrite, (float, int))
        check_positive("purgeWrite", purgeWrite)
        self._purgeWrite = purgeWrite

    @property
    def writeFormat(self):
        return self._writeFormat

    @writeFormat.setter
    def writeFormat(self, writeFormat):
        check_type('writeFormat', writeFormat, str)
        check_value('writeFormat', writeFormat, _WRITE_FORMAT_TYPES)
        self._writeFormat = writeFormat

    @property
    def writePrecision(self):
        return self._writePrecision

    @writePrecision.setter
    def writePrecision(self, writePrecision) -> None:
        check_type("writePrecision", writePrecision, int)
        check_positive("writePrecision", writePrecision, is_strict=True)
        self._writePrecision = writePrecision

    @property
    def writeCompression(self):
        return self._writeCompression

    @writeCompression.setter
    def writeCompression(self, writeCompression) -> None:
        check_type("writeCompression", writeCompression, bool)
        self._writeCompression = writeCompression

    @property
    def timeFormat(self):
        return self._timeFormat

    @timeFormat.setter
    def timeFormat(self, timeFormat):
        check_type('timeFormat', timeFormat, str)
        check_value('timeFormat', timeFormat, _TIME_FORMAT_TYPES)
        self._timeFormat = timeFormat

    @property
    def timePrecision(self):
        return self._timePrecision

    @timePrecision.setter
    def timePrecision(self, timePrecision) -> None:
        check_type("timePrecision", timePrecision, int)
        check_positive("timePrecision", timePrecision, is_strict=True)
        self._timePrecision = timePrecision

    @property
    def runTimeModifiable(self):
        return self._runTimeModifiable

    @runTimeModifiable.setter
    def runTimeModifiable(self, runTimeModifiable) -> None:
        check_type("runTimeModifiable", runTimeModifiable, bool)
        self._runTimeModifiable = runTimeModifiable

    @property
    def adjustTimeStep(self):
        return self._adjustTimeStep

    @adjustTimeStep.setter
    def adjustTimeStep(self, adjustTimeStep) -> None:
        check_type("adjustTimeStep", adjustTimeStep, bool)
        self._adjustTimeStep = adjustTimeStep

    @property
    def maxDeltaT(self):
        return self._maxDeltaT

    @maxDeltaT.setter
    def maxDeltaT(self, maxDeltaT) -> None:
        check_type("maxDeltaT", maxDeltaT, (float, int))
        check_positive("maxDeltaT", maxDeltaT, is_strict=True)
        self._maxDeltaT = maxDeltaT

    @property
    def maxCo(self):
        return self._maxCo

    @maxCo.setter
    def maxCo(self, maxCo) -> None:
        check_type("maxCo", maxCo, (float, int))
        check_positive("maxCo", maxCo, is_strict=True)
        self._maxCo = maxCo

    @property
    def maxAlphaCo(self):
        return self._maxAlphaCo

    @maxAlphaCo.setter
    def maxAlphaCo(self, maxAlphaCo) -> None:
        check_type("maxAlphaCo", maxAlphaCo, (float, int), none_ok=True)
        if (maxAlphaCo is not None):
            check_positive("maxAlphaCo", maxAlphaCo, is_strict=True)
        self._maxAlphaCo = maxAlphaCo

    @property
    def maxPowerVariation(self):
        return self._maxPowerVariation

    @maxPowerVariation.setter
    def maxPowerVariation(self, maxPowerVariation) -> None:
        check_type("maxPowerVariation", maxPowerVariation, (float, int))
        check_positive("maxPowerVariation", maxPowerVariation, is_strict=True)
        self._maxPowerVariation = maxPowerVariation

    @property
    def maxRelativeDeltaTIncrease(self):
        return self._maxRelativeDeltaTIncrease

    @maxRelativeDeltaTIncrease.setter
    def maxRelativeDeltaTIncrease(self, maxRelativeDeltaTIncrease) -> None:
        check_type("maxRelativeDeltaTIncrease", maxRelativeDeltaTIncrease, (float, int), none_ok=True)
        if (maxRelativeDeltaTIncrease is not None):
            check_positive("maxRelativeDeltaTIncrease", maxRelativeDeltaTIncrease, is_strict=True)
        self._maxRelativeDeltaTIncrease = maxRelativeDeltaTIncrease

    @property
    def minRelativeDeltaTDecrease(self):
        return self._minRelativeDeltaTDecrease

    @minRelativeDeltaTDecrease.setter
    def minRelativeDeltaTDecrease(self, minRelativeDeltaTDecrease) -> None:
        check_type("minRelativeDeltaTDecrease", minRelativeDeltaTDecrease, (float, int), none_ok=True)
        if (minRelativeDeltaTDecrease is not None):
            check_positive("minRelativeDeltaTDecrease", minRelativeDeltaTDecrease, is_strict=True)
        self._minRelativeDeltaTDecrease = minRelativeDeltaTDecrease

    @property
    def maxRelativePowerIncrease(self):
        return self._maxRelativePowerIncrease

    @maxRelativePowerIncrease.setter
    def maxRelativePowerIncrease(self, maxRelativePowerIncrease) -> None:
        check_type("maxRelativePowerIncrease", maxRelativePowerIncrease, (float, int), none_ok=True)
        if (maxRelativePowerIncrease is not None):
            check_positive("maxRelativePowerIncrease", maxRelativePowerIncrease, is_strict=True)
        self._maxRelativePowerIncrease = maxRelativePowerIncrease

    @property
    def maxRelativePowerDecrease(self):
        return self._maxRelativePowerDecrease

    @maxRelativePowerDecrease.setter
    def maxRelativePowerDecrease(self, maxRelativePowerDecrease) -> None:
        check_type("maxRelativePowerDecrease", maxRelativePowerDecrease, (float, int), none_ok=True)
        if (maxRelativePowerDecrease is not None):
            check_positive("maxRelativePowerDecrease", maxRelativePowerDecrease, is_strict=True)
        self._maxRelativePowerDecrease = maxRelativePowerDecrease

    @property
    def maxBurnupIncrease(self):
        return self._maxBurnupIncrease

    @maxBurnupIncrease.setter
    def maxBurnupIncrease(self, maxBurnupIncrease) -> None:
        check_type("maxBurnupIncrease", maxBurnupIncrease, (float, int), none_ok=True)
        if (maxBurnupIncrease is not None):
            check_positive("maxBurnupIncrease", maxBurnupIncrease, is_strict=True)
        self._maxBurnupIncrease = maxBurnupIncrease

    @property
    def maxAverageCreep(self):
        return self._maxAverageCreep

    @maxAverageCreep.setter
    def maxAverageCreep(self, maxAverageCreep) -> None:
        check_type("maxAverageCreep", maxAverageCreep, (float, int), none_ok=True)
        if (maxAverageCreep is not None):
            check_positive("maxAverageCreep", maxAverageCreep, is_strict=True)
        self._maxAverageCreep = maxAverageCreep

    @property
    def maxMaximumCreep(self):
        return self._maxMaximumCreep

    @maxMaximumCreep.setter
    def maxMaximumCreep(self, maxMaximumCreep) -> None:
        check_type("maxMaximumCreep", maxMaximumCreep, (float, int), none_ok=True)
        if (maxMaximumCreep is not None):
            check_positive("maxMaximumCreep", maxMaximumCreep, is_strict=True)
        self._maxMaximumCreep = maxMaximumCreep

    @property
    def maxFGR(self):
        return self._maxFGR

    @maxFGR.setter
    def maxFGR(self, maxFGR) -> None:
        check_type("maxFGR", maxFGR, (float, int), none_ok=True)
        if (maxFGR is not None):
            check_positive("maxFGR", maxFGR, is_strict=True)
        self._maxFGR = maxFGR

    @property
    def solveFMI(self):
        return self._solveFMI

    @solveFMI.setter
    def solveFMI(self, solveFMI) -> None:
        check_type("solveFMI", solveFMI, bool, none_ok=True)
        self._solveFMI = solveFMI

    @property
    def includeKineticEnergy(self):
        return self._includeKineticEnergy

    @includeKineticEnergy.setter
    def includeKineticEnergy(self, includeKineticEnergy) -> None:
        check_type("includeKineticEnergy", includeKineticEnergy, bool, none_ok=True)
        self._includeKineticEnergy = includeKineticEnergy

    @property
    def writeRestartFields(self):
        return self._writeRestartFields

    @writeRestartFields.setter
    def writeRestartFields(self, writeRestartFields) -> None:
        check_type("writeRestartFields", writeRestartFields, bool, none_ok=True)
        self._writeRestartFields = writeRestartFields

    @property
    def writeContinuityErrors(self):
        return self._writeContinuityErrors

    @writeContinuityErrors.setter
    def writeContinuityErrors(self, writeContinuityErrors) -> None:
        check_type("writeContinuityErrors", writeContinuityErrors, bool, none_ok=True)
        self._writeContinuityErrors = writeContinuityErrors

    @property
    def minDeltaT(self):
        return self._minDeltaT

    @minDeltaT.setter
    def minDeltaT(self, minDeltaT) -> None:
        check_type("minDeltaT", minDeltaT, (float, int), none_ok=True)
        if (minDeltaT is not None):
            check_positive("minDeltaT", minDeltaT, is_strict=True)
        self._minDeltaT = minDeltaT

    @property
    def maxCoTwoPhase(self):
        return self._maxCoTwoPhase

    @maxCoTwoPhase.setter
    def maxCoTwoPhase(self, maxCoTwoPhase) -> None:
        check_type("maxCoTwoPhase", maxCoTwoPhase, (float, int), none_ok=True)
        self._maxCoTwoPhase = maxCoTwoPhase

    @property
    def marginToPhaseChange(self):
        return self._marginToPhaseChange

    @marginToPhaseChange.setter
    def marginToPhaseChange(self, marginToPhaseChange) -> None:
        check_type("marginToPhaseChange", marginToPhaseChange, (float, int), none_ok=True)
        self._marginToPhaseChange = marginToPhaseChange

    @property
    def host(self):
        return self.simulationParametersFMI.host

    @host.setter
    def host(self, host) -> None:
        check_type("host", host, str)
        self.simulationParametersFMI.host = host

    @property
    def port(self):
        return self.simulationParametersFMI.port

    @port.setter
    def port(self, port) -> None:
        check_type("port", port, int)
        self.simulationParametersFMI.port = port

    @property
    def coupling(self):
        return self._coupling

    @coupling.setter
    def coupling(self, coupling):
        check_type("coupling", coupling, Coupling, none_ok=True)
        self._coupling = coupling

    @property
    def solvers(self):
        return self._solvers

    @solvers.setter
    def solvers(self, solvers):
        check_type("solvers", solvers, Solvers, none_ok=True)
        self._solvers = solvers

    @property
    def functionObjects(self):
        return self._functionObjects

    @functionObjects.setter
    def functionObjects(self, functionObjects):
        check_type("functionObjects", functionObjects, FunctionObjects, none_ok=True)
        self._functionObjects = functionObjects


    def add_function_object(self, functionObject: FunctionObject):
        check_type("functionObject", functionObject, FunctionObject)
        self.functionObjects.append(functionObject)


    def import_from_openfoam(self):
        foamFile = super().import_from_openfoam()
        attributes = self.get_attributes_as_list()
        for attr in attributes:
            if (attr in foamFile.keys()):
                setattr(self, attr, foamFile[attr])


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""

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

        if (self.application == "offbeat"):
            text += addParameter('adjustableTimeStep', self.adjustTimeStep, isAddExtraLine=True)
        else:
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

        if (self.marginToPhaseChange is not None):
            text += addParameter('marginToPhaseChange', self.marginToPhaseChange, isAddExtraLine=True)
        text += addParameter('maxPowerVariation', self.maxPowerVariation, isAddExtraLine=True)

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

        # if (self.coupling is not None):
        #     text += self.coupling.write_to_controlDict()

        if (self.solvers is not None):
            text += self.solvers.write_to_controlDict()

        if (self.functionObjects is not None):
            for functionObject in self.functionObjects:
                if (isinstance(functionObject, FMUSimulator)):
                    self.libs.append('externalComm')
            if (not self.libs.is_empty):
                text += f"{self.libs!r}\n"
            text += f"{self.functionObjects!r}\n"

        return(text)


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
