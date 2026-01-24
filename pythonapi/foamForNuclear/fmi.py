import multiprocessing
import os
import time
import pandas as pd
import numpy as np
import psutil
import uncertainties as un

from foamForNuclear.checkvalue import check_type

try:
    import FMU4FOAM
    from pyfmi import load_fmu # type: ignore


    class PyFMIContainer(FMU4FOAM.FMUBase):
        """
        PyFMI Container used to simulate an FMU using the pyfmi package.

        Python API to manage the FMU in GeN-Foam. Use PyFMI instead of OMSimulator to
        be able to reset the FMU state.
        Inspired by https://github.com/DLR-RY/FMU4FOAM/tree/master/examples.

        Usage::

            from foamForNuclear.fmi import PyFMIContainer

            class ControlerExample(PyFMIContainer):
                def __init__(self, endTime, jsonFile):
                    super().__init__(
                        endTime,
                        jsonFile,
                        # Define by the user
                        fmuName='ControlerExample.fmu',
                        outputFilename='ControlerExample.csv'
                    )

        Parameters
        ----------
        endTime : float
            Filled by the simulator
        jsonFile: str
            Filled by the simulator
        fmuName : str
            Name of the FMU file
        outputFilename : str
            Name of the output file to store FMU results, extension: `.csv`
        writeInterval : float
            Time distance between 2 writing in CSV file (default to 0, which is
            equivalent to writing every time step)
        substep : float
            Number of sub step performed by the FMU (default `1`)


        Author: Thomas Guilbaud, EPFL/Transmutex SA, 04/2023
        """
        def __init__(
                self,
                endTime,
                jsonFile,
                fmuName: str,
                outputFilename: str,
                log_level: int=2,
                writeInterval: float=0,
                substep: float=1,
                parametersToRecord: list[str]=None
            ) -> None:
            super().__init__(endTime, jsonFile)

            # Some useful variables
            self.outputFilename = outputFilename
            # self.step_size = 0.1
            self.current_time = 0
            self.previous_time = self.current_time
            self.indexWrite = 0
            self.writeInterval = writeInterval
            self.substep = substep

            # Create the FMU
            self.fmu1 = load_fmu(fmuName)
            # opts = self.fmu1.simulate_options()
            # opts['result_handling'] = "none"
            # print("keys", opts.keys(), flush=True)
            # opts['result_max_size'] = 1 * 1e9
            # self.fmu1.options = opts
            # help(self.fmu1)
            self.fmu1.set_log_level(log_level)

            # Start initialization mode
            self.fmu1.initialize()

            # Extract all the parameters and variables composing the FMU
            self.parameterList = parametersToRecord
            if (self.parameterList is None):
                self.parameterList = [key for key in self.fmu1.get_model_variables() if key[:3] != "_D_"]

            # Create the output buffer (write into a csv)
            self.results = pd.DataFrame(columns=['time'] + self.parameterList)

            # If a simulation has already been performed -> restart
            if (os.path.exists(self.outputFilename)):
                # Read the last line of the output file
                data = pd.read_csv(self.outputFilename)
                lastIdx = data.index.max()

                # Set all the variables for the restart based on the previous simulation
                for d in data.columns:
                    val = np.float64(data[d][lastIdx])
                    if (d != "time"):
                        self.setVar(d, val)
                    else:
                        self.current_time = val

                # Reset buffer
                self.results.loc[0] = data.loc[lastIdx]


        def setVar(self, key: str, val: float) -> None:
            valueRef = self.getValueRef(key)
            self.fmu1.set_real([valueRef], [val])


        def getVar(self, key: str) -> float:
            valueRef = self.getValueRef(key)
            return(self.fmu1.get_real([valueRef])[0])


        def getValueRef(self, key: str) -> int:
            return(self.fmu1.get_variable_valueref(key))


        def stepUntil(self, t: float) -> None:
            # Compute the time step for the FMU
            step_size = (t - self.current_time) / 1.0
            while (self.current_time < t):
                # Make a step with the FMU
                self.fmu1.do_step(current_t=self.current_time, step_size=step_size, new_step=True)
                self.current_time += step_size

            # Update the current time
            self.current_time = t

            # Write the results
            if (self.previous_time + self.writeInterval < self.current_time):
                isFirstWrite = pd.isnull(self.results.index.max())
                self.results.loc[0] = [self.current_time] + [self.getVar(key) for key in self.parameterList]
                if (isFirstWrite):
                    self.results.to_csv(self.outputFilename, index=False)
                else:
                    self.results.to_csv(self.outputFilename, index=False, mode='a', header=False)

                self.previous_time = self.current_time


        def __del__(self) -> None:
            pass


except ModuleNotFoundError:
    # Error handling
    # print("No module FMU4FOAM or pyfmi found.")
    # sys.exit(0)
    pass



try:
    import FMU4FOAM
    from OMSimulator import OMSimulator


    class OMSimulatorContainer(FMU4FOAM.FMUBase):
        """
        OMSimulator Container used to simulate an FMU using the OMSimulator
        package.

        Usage::

            from foamForNuclear.fmi import OMSimulator

            class ControlerExample(OMSimulatorContainer):
                def __init__(self, endTime, jsonFile):
                    super().__init__(
                        endTime,
                        jsonFile,
                        # Define by the user
                        fmuName='ControlerExample.fmu',
                        outputFilename='ControlerExample.csv'
                    )

        Parameters
        ----------
        endTime : float
            Filled by the simulator
        jsonFile: str
            Filled by the simulator
        fmuName : str
            Name of the FMU file
        outputFilename : str
            Name of the output file to store FMU results, extension: `.csv`
        """
        def __init__(
                self,
                endTime,
                jsonFile,
                fmuName: str,
                outputFilename: str
            ):
            super().__init__(endTime, jsonFile)

            self.oms = OMSimulator()
            self.oms.setTempDirectory("./temp/")
            self.oms.newModel("model")
            self.oms.addSystem("model.root", self.oms.system_wc)
            self.oms.setSolver("model", self.oms.solver_sc_cvode)

            # instantiate FMUs
            self.oms.addSubModel("model.root.system1", fmuName)

            # simulation settings
            self.oms.setResultFile("model", outputFilename)
            self.oms.setStopTime("model", endTime)
            # self.oms.setFixedStepSize("model.root", 0.001)
            self.oms.removeSignalsFromResults("model", ".*_D_.*")

            self.oms.instantiate("model")
            #self.oms.setReal("model.root.system1.momentumSource", 1)

            self.oms.initialize("model")
            #self.oms.simulate("model")

        def setVar(self, key: str, val: float) -> None:
            return self.oms.setReal(key,val)

        def getVar(self, key: str) -> float:
            return self.oms.getReal(key)[0]

        def stepUntil(self, t):
            self.oms.stepUntil("model", t)

        def __del__(self):
            self.oms.terminate("model")
            self.oms.delete("model")


except ModuleNotFoundError:
    # Error handling
    pass

except ImportError:
    # Error handling
    pass


try:
    import FMU4FOAM
    import fmpy


    class FMPyContainer(FMU4FOAM.FMUBase):
        """
        FMPy Container used to simulate an FMU using the pyfmi package.

        Python API to manage the FMU in GeN-Foam. Use FMPy instead of OMSimulator to
        be able to reset the FMU state.
        Inspired by https://github.com/DLR-RY/FMU4FOAM/tree/master/examples.

        Usage::

            from foamForNuclear.fmi import FMPyContainer

            class ControlerExample(FMPyContainer):
                def __init__(self, endTime, jsonFile):
                    super().__init__(
                        endTime,
                        jsonFile,
                        # Define by the user
                        fmuName='ControlerExample.fmu',
                        outputFilename='ControlerExample.csv'
                    )

        Parameters
        ----------
        endTime : float
            Filled by the simulator
        jsonFile: str
            Filled by the simulator
        fmuName : str
            Name of the FMU file
        outputFilename : str
            Name of the output file to store FMU results, extension: `.csv`
        writeInterval : float
            Time distance between 2 writing in CSV file (default `0`, which is
            equivalent to writing every time step)
        substep : float
            Number of sub step performed by the FMU (default `1`)
        maxGb : float
            Maximum memory allocation of the FMU (default `2` Gb). If the memory
            allocation of the FMU exceed this limit, the FMU state is saved, the
            FMU is killed and restarted, and the FMU state restored.
        parametersToRecord : list[str]
            List of parameter to store on disk (default `None` which is
            equivalent of writing all variables, can be time consumming).
        initialParametersBeforeFMUInit : dict
            Initialization of some specific parameters in the FMU before the FMU
            initialization. E.g set the port number: `{"port": 8001}`.
        """
        def __init__(
                self,
                endTime,
                jsonFile,
                fmuName: str,
                outputFilename: str,
                startTime: float=0,
                log_level: int=2,
                writeInterval: float=0,
                substep: float=1,
                maxGb: float=2,
                parametersToRecord: list[str]=None,
                initialParametersBeforeFMUInit: dict=None
            ) -> None:
            if (jsonFile is not None):
                super().__init__(endTime, jsonFile)

            # Some useful variables
            self.endTime = endTime
            self.outputFilename = outputFilename
            # self.step_size = 0.1
            self.current_time = startTime
            self.previous_time = self.current_time
            self.indexWrite = 0
            self.writeInterval = writeInterval
            self.substep = substep
            self.maxGb = maxGb

            # Store restart time step if memory leak
            self.restartPoints = []

            # Extract the FMU to a temporary directory
            self.unzipdir = fmpy.extract(fmuName)

            # Read the model description
            self.model_description = fmpy.read_model_description(self.unzipdir)

            # Instantiate the FMU beforehand, so we can keep it alive
            self.fmu1 = fmpy.instantiate_fmu(unzipdir=self.unzipdir, model_description=self.model_description)

            # Extract value reference and types
            self.valueRef = {
                variable.name: [variable.valueReference]
                for variable in self.model_description.modelVariables
            }
            self.valueType = {
                variable.name: variable.type
                for variable in self.model_description.modelVariables
            }

            # Init values before initialization of the FMU
            if (initialParametersBeforeFMUInit is not None):
                for parameterName, value in initialParametersBeforeFMUInit.items():
                    self.setVar(parameterName, value)

            # Start initialization mode
            self.fmu1.setupExperiment(stopTime=endTime)
            self.fmu1.enterInitializationMode()
            self.fmu1.exitInitializationMode()


            # Extract all the parameters and variables composing the FMU
            self.parameterList = parametersToRecord
            if (self.parameterList is None):
                self.parameterList = [key for key in self.valueRef.keys() if key[:3] != "_D_"]

            # Create the output buffer (write into a csv)
            self.results = pd.DataFrame(columns=['time'] + self.parameterList)

            # If a simulation has already been performed -> restart
            if (os.path.exists(self.outputFilename)):
                print(f"Read {self.outputFilename} ...", flush=True)
                # Read the last line of the output file
                data = pd.read_csv(self.outputFilename)
                lastIdx = data.index.max()

                # Set all the variables for the restart based on the previous simulation
                for d in data.columns:
                    val = data[d][lastIdx]
                    if (d != "time"):
                        try:
                            self.setVar(d, val)
                        except:
                            print(f"Cannot set {d} to {val}", flush=True)
                    else:
                        self.current_time = val

                # Reset buffer
                self.results.loc[0] = data.loc[lastIdx]


        def getVar(self, key: str):
            valueRef = self.valueRef[key]
            valueType = self.valueType[key]

            if (valueType == 'Real'):
                return(self.fmu1.getReal(valueRef)[0])
            elif (valueType == 'Integer'):
                return(self.fmu1.getInteger(valueRef)[0])
            elif (valueType == 'Boolean'):
                return(self.fmu1.getBoolean(valueRef)[0])
            elif (valueType == 'String'):
                return(self.fmu1.getString(valueRef)[0].decode("utf-8"))
            elif (valueType == 'Enumeration'):
                return(None)

            msg = f"Type {valueType} not supported. Available types: 'Real', 'Integer', 'Boolean', 'String'"
            raise ValueError(msg)


        def setVar(self, key: str, val: float | bool | int | str) -> None:
            valueRef = self.valueRef[key]
            valueType = self.valueType[key]

            if (valueType == 'Real'):
                self.fmu1.setReal(valueRef, [val])
            elif (valueType == 'Integer'):
                self.fmu1.setInteger(valueRef, [int(val)])
            elif (valueType == 'Boolean'):
                self.fmu1.setBoolean(valueRef, [bool(val)])
            elif (valueType == 'String'):
                self.fmu1.setString(valueRef, [val])
            elif (valueType == 'Enumeration'):
                pass
            else:
                msg = f"Type {valueType} not supported. Available types: 'Real', 'Integer', 'Boolean', 'String', 'Enumeration'"
                raise ValueError(msg)


        def stepUntil(self, t: float, newStep: bool=True) -> None:
            # Compute the time step for the FMU
            step_size = (t - self.current_time) / self.substep
            current_time_temp = self.current_time
            while (current_time_temp < t):
                # Make a step with the FMU
                self.fmu1.doStep(
                    currentCommunicationPoint=current_time_temp,
                    communicationStepSize=step_size,
                    noSetFMUStatePriorToCurrentPoint=newStep
                )
                current_time_temp += step_size

            # If memory leak, max 1 Gb, save state, kill the FMU and restart it
            if psutil.Process().memory_info().rss > self.maxGb * 1024**3:
                # Save
                # state = self.fmu1.getFMUstate()
                allValues = {
                    variable.name: self.getVar(variable.name)
                    for variable in self.model_description.modelVariables
                }

                # Kill
                self.fmu1.terminate()
                self.fmu1.freeInstance()

                # Restart
                self.fmu1 = fmpy.instantiate_fmu(self.unzipdir, self.model_description)
                # self.fmu1.setFMUstate(state)
                self.fmu1.setupExperiment(
                    startTime=self.current_time,
                    stopTime=self.endTime
                )
                self.fmu1.enterInitializationMode()
                self.fmu1.exitInitializationMode()

                # Reset all the values
                for variableName, value in allValues.items():
                    self.setVar(variableName, value)

                print(f"Restart FMU at t = {self.current_time:g}", flush=True)
                self.restartPoints.append(self.current_time)


            if (newStep):
                # Update the current time
                self.current_time = t

                # Write the results
                if (self.previous_time + self.writeInterval < self.current_time):
                    isFirstWrite = pd.isnull(self.results.index.max())
                    self.results.loc[0] = [self.current_time] + [self.getVar(key) for key in self.parameterList]
                    if (isFirstWrite):
                        self.results.to_csv(self.outputFilename, index=False)
                    else:
                        self.results.to_csv(self.outputFilename, index=False, mode='a', header=False)

                    self.previous_time = self.current_time


        def __del__(self) -> None:
            self.fmu1.terminate()
            self.fmu1.freeInstance()


    class FMPyMasterRunner:
        """
        Master runner for coupled FMU simulations. Uses an adaptative time step
        and check model convergence.

        Example ::

            FMPyMasterRunner(
                fmuNames={
                    "core": "Core.fmu",
                    "controller": "Controller.fmu"
                },
                connections=[
                    # Source                            | Destination
                    # FMU name     FMI port             | FMU name      FMI port
                    ("core",       "powerOutput",         "controller", "powerInput"        ),
                    ("controller", "extReactivityOutput", "core",       "extReactivityInput")
                ]
            )

        Parameters
        ----------
        fmuNames : dict[str]
            Dict of FMU names. Order matter as it is also the initialization
            order and execution order
        connections : list[tuple]
            The tuple must have the following form:
            `(fmuNameSrc, fmiPortNameSrc, fmuNameDest, fmiPortNameDest)`
        """
        def __init__(
                self,
                fmuNames: dict[str],
                connections: list[tuple],
                startTime: float,
                endTime: float,
                stepSize: float,
                minStepSize: float=0.05,
                maxStepSize: float=5,
                minNImplicitSteps: int=5,
                maxNImplicitSteps: int=100,  # if see oscillations when reducing time step, increase n implicit
                adjustStepSize: bool=True,
                writeInterval: float=0,
                relaxationFactorTimeStep: float=1,
                parametersToRecord: dict[list[str]]=None,
                initialParametersBeforeFMUInit: dict[dict]=None,
                outputFilenames: dict[str, str]=None
            ):
            self.startTime = startTime
            self.endTime = endTime
            self.stepSize = stepSize
            self.minStepSize = minStepSize
            self.maxStepSize = maxStepSize
            self.minNImplicitSteps = minNImplicitSteps
            self.maxNImplicitSteps = maxNImplicitSteps
            self.adjustStepSize = adjustStepSize
            self.relaxationFactorTimeStep = relaxationFactorTimeStep

            self.fmus = {}

            # Initialization of FMUs
            nFMU = 0
            for fmuName, fmuFileName in fmuNames.items():
                print(f"Instantiate FMU {fmuName}")
                parametersToRecord_i = None
                if (parametersToRecord is not None and fmuName in parametersToRecord.keys()):
                    parametersToRecord_i = parametersToRecord[fmuName]

                initialParametersBeforeFMUInit_i = None
                if (initialParametersBeforeFMUInit is not None and fmuName in initialParametersBeforeFMUInit.keys()):
                    initialParametersBeforeFMUInit_i = initialParametersBeforeFMUInit[fmuName]

                outputFilename_i = f"{fmuFileName[:-4]}_{nFMU}.csv"
                if (outputFilenames is not None and fmuName in outputFilenames.keys()):
                    outputFilename_i = outputFilenames[fmuName]

                self.fmus[fmuName] = FMPyContainer(
                    endTime=endTime,
                    jsonFile=None,
                    fmuName=fmuFileName,
                    outputFilename=outputFilename_i,
                    startTime=startTime,
                    writeInterval=writeInterval,
                    parametersToRecord=parametersToRecord_i,
                    initialParametersBeforeFMUInit=initialParametersBeforeFMUInit_i
                )

                # # Overwrite the oscmd that is by default "bash -i", when using
                # # multiple OpenFOAM FMUs, the code crash
                # if ("oscmd" in list(self.fmus[fmuName].valueRef.keys())):
                #     self.fmus[fmuName].setVar("oscmd", "bash")

                nFMU += 1


            self.connections = connections


        def add_connection(
                self,
                fmuNameSrc: str,
                fmiPortNameSrc: str,
                fmuNameDest: str,
                fmiPortNameDest: str
            ):
            check_type("fmuNameSrc", fmuNameSrc, str)
            check_type("fmiPortNameSrc", fmiPortNameSrc, str)
            check_type("fmuNameDest", fmuNameDest, str)
            check_type("fmiPortNameDest", fmiPortNameDest, str)

            self.connections.append(
                (fmuNameSrc, fmiPortNameSrc, fmuNameDest, fmiPortNameDest)
            )

        def step_all(self, currentTime, stepSize, isLastStep):
            for fmu in self.fmus.values():
                status = fmu.stepUntil(
                    t=currentTime + stepSize,
                    newStep=isLastStep
                )

        def exchange_variables(self):
            for connection in self.connections:
                srcFMU = self.fmus[connection[0]]
                destFMU = self.fmus[connection[2]]
                value = srcFMU.getVar(connection[1])
                destFMU.setVar(connection[3], value)

        def is_converge_flag_present(self):
            for name, conn in self.pipes.items():
                conn.send({"cmd": "valueRef"})

            return any(["isConverged" in conn.recv()["valueRef"].keys() for name, conn in self.pipes.items()])

        def evaluate_min_step_size(self, evaluated_minStepSize):
            for fmu in self.fmus.values():
                if ("dt" in fmu.valueRef.keys()):
                    evaluated_minStepSize = min(evaluated_minStepSize, fmu.getVar("dt"))

            return evaluated_minStepSize

        def is_all_model_converged(self):
            isAllModelConverged = True
            for fmu in self.fmus.values():
                if ("isConverged" in fmu.valueRef.keys()):
                    isAllModelConverged = isAllModelConverged and fmu.getVar("isConverged") == 1

            return(isAllModelConverged)

        def simulate(self):
            exec_start_time = time.time()

            # Check at least one FMU has isConvergedFlag
            isConvergedFlag = False
            for fmu in self.fmus.values():
                if ("isConverged" in fmu.valueRef.keys()):
                    isConvergedFlag = True
                    break
            print("Converged flag present in at least one FMUs", isConvergedFlag)

            # Loop inspired by https://github.com/modelon-community/PyFMI/blob/master/src/pyfmi/master.pyx
            current_time = self.startTime # Why it works? OpenFOAM current time is the time + the time step?
            new_step_size = self.stepSize
            convergedInNiter = 0
            nLastPoints = 10
            exec_time_list = np.zeros(nLastPoints)
            iter_time = 0
            # Time loop
            while (current_time + self.stepSize < self.endTime):
                # Time recording
                exec_time = time.time() - exec_start_time
                exec_time_list[iter_time] = exec_time
                exec_time_diff = [tf-ti for ti, tf in zip(exec_time_list[:-1], exec_time_list[1:]) if tf-ti > 0] # Approximation, missing (end-first) term
                end_time = (self.endTime-self.startTime)/(current_time-self.startTime) * exec_time if current_time != self.startTime else exec_time
                avgTimeStep = un.ufloat(np.mean(exec_time_diff), np.std(exec_time_diff))
                print("Time= {:.3f}s (dt={:<8.3e}s)  Exec= {:.1f}s   End= {:.1f}s   Remain= {:.1f}s   Eff= {} sim/real (N iter= {})".format(
                    current_time,
                    self.stepSize,
                    exec_time,
                    end_time,
                    end_time - exec_time,
                    self.stepSize/avgTimeStep,
                    convergedInNiter
                ))
                iter_time = iter_time + 1 if nLastPoints > iter_time+1 else 0

                # Implicit loop
                evaluated_minStepSize = 1e32 if self.adjustStepSize else self.stepSize
                if (self.maxNImplicitSteps > 0):
                    iter_i = 0
                    while (iter_i <= self.maxNImplicitSteps - 1):
                        # Each FMU do one step with their data
                        # This do_step function receive data from FMUs and sends data from
                        # buffer. Calling do_step release the FMUs because they received a
                        # response from do_step
                        isLastStep = iter_i == self.maxNImplicitSteps - 1 # True only at last step
                        self.step_all(currentTime=current_time, stepSize=self.stepSize, isLastStep=isLastStep)

                        # Exchange info between FMUs for the next step
                        self.exchange_variables()

                        # Look for smallest time step among solvers
                        if (self.adjustStepSize):
                            evaluated_minStepSize = self.evaluate_min_step_size(evaluated_minStepSize)

                        # Optional convergence flags
                        if (
                            isConvergedFlag
                            and not isLastStep
                            and self.is_all_model_converged()
                            and self.minNImplicitSteps <= iter_i
                        ):
                            convergedInNiter = iter_i+1
                            # print(f"Converged in {iter_i+1} iterations")
                            iter_i = self.maxNImplicitSteps - 2
                            # No break, needs to do a last time step with
                            # 'new_step' to True
                        elif (not isLastStep):
                            convergedInNiter = iter_i+1

                        iter_i += 1

                # Adapt new step to time step boundaries provided by the user
                if (self.adjustStepSize):
                    if (evaluated_minStepSize > 0):
                        new_step_size = evaluated_minStepSize
                        new_step_size = min(max(new_step_size, self.minStepSize), self.maxStepSize) # Correct if new step size is too small

                    # Update time step
                    if (current_time + new_step_size < self.endTime):
                        current_time += self.stepSize # new_step_size #
                        alpha = self.relaxationFactorTimeStep
                        self.stepSize = alpha * new_step_size + (1-alpha) * self.stepSize
                        # print(f'Min step size {evaluated_minStepSize:.5g}, step size {step_size:.5g}, new_step size {new_step_size:.5g}; {current_time}={connection[0].get("t")[0]}\n')
                    else:
                        # In case step-size exceeds the end time
                        self.stepSize = self.endTime - current_time
                        current_time = self.endTime
                        print('Last step', self.stepSize, current_time)

                else:
                    current_time += self.stepSize

            print(f"Simulation done, total time = {time.time() - exec_start_time:.3f} s")


    def fmu_worker(fmu_args, conn):
        """
        FMU worker used only for multi-processing. This function is alive using
        the `while True` loop and a pipe connection like TCP.

        Parameters
        ----------
        conn : multiprocessing.Pipe
        """
        fmu = FMPyContainer(**fmu_args)

        while True:
            msg = conn.recv()
            if msg["cmd"] == "step":
                result = fmu.stepUntil(t=msg["time"], newStep=msg["newStep"])
                conn.send({"status": result})
            elif msg["cmd"] == "getVar":
                value = fmu.getVar(msg["varName"])
                conn.send({"value": value})
            elif msg["cmd"] == "setVar":
                fmu.setVar(msg["varName"], msg["value"])
                conn.send({"status": "ok"})
            elif msg["cmd"] == "valueRef":
                valueRef = fmu.valueRef
                conn.send({"valueRef": valueRef})
            elif msg["cmd"] == "terminate":
                # fmu.terminate()
                conn.close()
                break


    class FMPyMasterRunnerParallel:
        """
        Master runner for coupled FMU simulations in parallel. Uses an
        adaptative time step and check model convergence.

        Example ::

            FMPyMasterRunnerParallel(
                fmuNames={
                    "core": "Core.fmu",
                    "controller": "Controller.fmu"
                },
                connections=[
                    # Source                            | Destination
                    # FMU name     FMI port             | FMU name      FMI port
                    ("core",       "powerOutput",         "controller", "powerInput"        ),
                    ("controller", "extReactivityOutput", "core",       "extReactivityInput")
                ]
            )

        Parameters
        ----------
        fmuNames : dict[str]
            Dict of FMU names. Order matter as it is also the initialization
            order and execution order
        connections : list[tuple]
            The tuple must have the following form:
            `(fmuNameSrc, fmiPortNameSrc, fmuNameDest, fmiPortNameDest)`
        """
        def __init__(
                self,
                fmuNames: dict[str],
                connections: list[tuple],
                startTime: float,
                endTime: float,
                stepSize: float,
                minStepSize: float=0.05,
                maxStepSize: float=5,
                minNImplicitSteps: int=5,
                maxNImplicitSteps: int=100,  # if see oscillations when reducing time step, increase n implicit
                adjustStepSize: bool=True,
                writeInterval: float=0,
                relaxationFactorTimeStep: float=1,
                parametersToRecord: dict[list[str]]=None,
                initialParametersBeforeFMUInit: dict={},
                outputFilenames: dict[str, str]=None,
            ):
            self.startTime = startTime
            self.endTime = endTime
            self.stepSize = stepSize
            self.minStepSize = minStepSize
            self.maxStepSize = maxStepSize
            self.minNImplicitSteps = minNImplicitSteps
            self.maxNImplicitSteps = maxNImplicitSteps
            self.adjustStepSize = adjustStepSize
            self.relaxationFactorTimeStep = relaxationFactorTimeStep

            # For parallel computing management
            self.processes = {}
            self.pipes = {}

            # Initialization of FMUs
            nFMU = 0
            for fmuName, fmuFileName in fmuNames.items():
                print(f"Instantiate FMU {fmuName}")
                parametersToRecord_i = None
                if (parametersToRecord is not None and fmuName in parametersToRecord.keys()):
                    parametersToRecord_i = parametersToRecord[fmuName]

                initialParametersBeforeFMUInit_i = None
                if (initialParametersBeforeFMUInit is not None and fmuName in initialParametersBeforeFMUInit.keys()):
                    initialParametersBeforeFMUInit_i = initialParametersBeforeFMUInit[fmuName]

                outputFilename_i = f"{fmuFileName[:-4]}_{nFMU}.csv"
                if (outputFilenames is not None and fmuName in outputFilenames.keys()):
                    outputFilename_i = outputFilenames[fmuName]

                parent_conn, child_conn = multiprocessing.Pipe()
                fmu_args = {
                    "endTime": endTime,
                    "jsonFile": None,
                    "fmuName": fmuFileName,
                    "outputFilename": outputFilename_i,
                    "startTime": startTime,
                    "writeInterval": writeInterval,
                    "parametersToRecord": parametersToRecord_i,
                    "initialParametersBeforeFMUInit": initialParametersBeforeFMUInit_i
                }
                p = multiprocessing.Process(target=fmu_worker, args=(fmu_args, child_conn))
                p.start()
                self.processes[fmuName] = p
                self.pipes[fmuName] = parent_conn


                # # Overwrite the oscmd that is by default "bash -i", when using
                # # multiple OpenFOAM FMUs, the code crash
                # if ("oscmd" in list(self.fmus[fmuName].valueRef.keys())):
                #     self.fmus[fmuName].setVar("oscmd", "bash")

                nFMU += 1


            self.connections = connections


        def add_connection(
                self,
                fmuNameSrc: str,
                fmiPortNameSrc: str,
                fmuNameDest: str,
                fmiPortNameDest: str
            ):
            check_type("fmuNameSrc", fmuNameSrc, str)
            check_type("fmiPortNameSrc", fmiPortNameSrc, str)
            check_type("fmuNameDest", fmuNameDest, str)
            check_type("fmiPortNameDest", fmiPortNameDest, str)

            self.connections.append(
                (fmuNameSrc, fmiPortNameSrc, fmuNameDest, fmiPortNameDest)
            )

        def step_all(self, currentTime, stepSize, isLastStep):
            for name, conn in self.pipes.items():
                conn.send({"cmd": "step", "time": currentTime + stepSize, "newStep": isLastStep})

            return {name: conn.recv() for name, conn in self.pipes.items()}

        def exchange_variables(self, connections):
            for src_name, src_var, dest_name, dest_var in connections:
                self.pipes[src_name].send({"cmd": "getVar", "varName": src_var})
                value = self.pipes[src_name].recv()["value"]
                self.pipes[dest_name].send({"cmd": "setVar", "varName": dest_var, "value": value})
                self.pipes[dest_name].recv()

        def terminate_all(self):
            for conn in self.pipes.values():
                conn.send({"cmd": "terminate"})

        def is_converge_flag_present(self):
            for name, conn in self.pipes.items():
                conn.send({"cmd": "valueRef"})

            return any(["isConverged" in conn.recv()["valueRef"].keys() for name, conn in self.pipes.items()])

        def evaluate_min_step_size(self, evaluated_minStepSize):
            for name, conn in self.pipes.items():
                conn.send({"cmd": "valueRef"})
                keys = conn.recv()["valueRef"].keys()

                if ("dt" in keys):
                    conn.send({"cmd": "getVar", "varName": "dt"})
                    dt = conn.recv()["value"]
                    evaluated_minStepSize = min(evaluated_minStepSize, dt)

            return evaluated_minStepSize

        def is_all_model_converged(self):
            isAllModelConverged = True
            # nConverged = 0
            # nTotal = 0

            for name, conn in self.pipes.items():
                conn.send({"cmd": "valueRef"})
                keys = conn.recv()["valueRef"].keys()

                if ("isConverged" in keys):
                    conn.send({"cmd": "getVar", "varName": "isConverged"})
                    isConverged = conn.recv()["value"]

                    isAllModelConverged = isAllModelConverged and isConverged == 1

                    # nConverged += 1 if isConverged == 1 else 0
                    # nTotal += 1

            return(isAllModelConverged)

        def simulate(self):
            exec_start_time = time.time()

            isConvergedFlag = self.is_converge_flag_present()
            print("Converged flag present in at least one FMUs", isConvergedFlag)

            # Loop inspired by https://github.com/modelon-community/PyFMI/blob/master/src/pyfmi/master.pyx
            current_time = self.startTime # + self.stepSize # Why it works? OpenFOAM current time is the time + the time step?
            new_step_size = self.stepSize
            convergedInNiter = 0
            nLastPoints = 10
            exec_time_list = np.zeros(nLastPoints)
            iter_time = 0
            # Time loop
            while (current_time + self.stepSize < self.endTime):
                # Time information
                exec_time = time.time() - exec_start_time
                exec_time_list[iter_time] = exec_time
                exec_time_diff = [tf-ti for ti, tf in zip(exec_time_list[:-1], exec_time_list[1:]) if tf-ti > 0] # Approximation, missing (end-first) term
                end_time = (self.endTime-self.startTime)/(current_time-self.startTime) * exec_time if current_time != self.startTime else exec_time
                avgTimeStep = un.ufloat(np.mean(exec_time_diff), np.std(exec_time_diff))
                print("Time= {:.3f}s (dt={:<8.3e}s)  Exec= {:.1f}s   End= {:.1f}s   Remain= {:.1f}s   Eff= {} sim/real (N iter= {})".format(
                    current_time,
                    self.stepSize,
                    exec_time,
                    end_time,
                    end_time - exec_time,
                    self.stepSize/avgTimeStep,
                    convergedInNiter
                ))
                iter_time = iter_time + 1 if nLastPoints > iter_time+1 else 0

                # Implicit loop
                evaluated_minStepSize = 1e32 if self.adjustStepSize else self.stepSize
                if (self.maxNImplicitSteps > 0):
                    iter_i = 0
                    while (iter_i <= self.maxNImplicitSteps - 1):
                        # Each FMU do one step with their data
                        # This do_step function receive data from FMUs and sends data from
                        # buffer. Calling do_step release the FMUs because they received a
                        # response from do_step
                        isLastStep = iter_i == self.maxNImplicitSteps - 1 # True only at last step

                        self.step_all(currentTime=current_time, stepSize=self.stepSize, isLastStep=isLastStep)

                        # Exchange info between FMUs for the next step
                        self.exchange_variables(self.connections)

                        # Look for smallest time step among solvers
                        if (self.adjustStepSize):
                            evaluated_minStepSize = self.evaluate_min_step_size(evaluated_minStepSize)

                        # Optional convergence flags
                        if (
                            isConvergedFlag
                            and not isLastStep
                            and self.is_all_model_converged()
                            and self.minNImplicitSteps <= iter_i
                        ):
                            convergedInNiter = iter_i+1
                            # print(f"Converged in {iter_i+1} iterations")
                            iter_i = self.maxNImplicitSteps - 2
                            # No break, needs to do a last time step with
                            # 'new_step' to True
                        elif (not isLastStep):
                            convergedInNiter = iter_i+1

                        iter_i += 1

                # Adapt new step to time step boundaries provided by the user
                if (self.adjustStepSize):
                    if (evaluated_minStepSize > 0):
                        new_step_size = evaluated_minStepSize
                        new_step_size = min(max(new_step_size, self.minStepSize), self.maxStepSize) # Correct if new step size is too small

                    # Update time step
                    if (current_time + new_step_size < self.endTime):
                        current_time += self.stepSize # new_step_size #
                        alpha = self.relaxationFactorTimeStep
                        self.stepSize = alpha * new_step_size + (1-alpha) * self.stepSize
                        # print(f'Min step size {evaluated_minStepSize:.5g}, step size {self.stepSize:.5g}, new_step size {new_step_size:.5g}; {current_time}\n')

                    else:
                        # In case step-size exceeds the end time
                        self.stepSize = self.endTime - current_time
                        current_time = self.endTime
                        print('Last step', self.stepSize, current_time)

                else:
                    current_time += self.stepSize

            print(f"Simulation done, total time = {time.time() - exec_start_time:.3f} s")

            self.terminate_all()


except ModuleNotFoundError:
    # Error handling
    # print("No module FMU4FOAM or fmpy found.")
    # sys.exit(0)
    pass
