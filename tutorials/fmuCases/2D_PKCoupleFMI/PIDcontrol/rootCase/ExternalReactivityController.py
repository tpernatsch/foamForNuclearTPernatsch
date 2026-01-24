"""
Python API to manage the FMU in GeN-Foam. Use PyFMI instead of OMSimulator to
be able to reset the FMU state.
Inspired by https://github.com/DLR-RY/FMU4FOAM/tree/master/examples.

Author: Thomas Guilbaud, EPFL/Transmutex SA, 04/2023
"""

#=============================================================================*
# Imports
#=============================================================================*

import FMU4FOAM

import fmpy
import os
import pandas as pd
import psutil

#=============================================================================*


class ExternalReactivityController(FMU4FOAM.FMUBase):
    """
    Duplicated from foamForNuclear.FMPyContainer
    """

    def __init__(self, endTime, jsonFile) -> None:
        super().__init__(endTime, jsonFile)

        # Some useful variables
        self.endTime = endTime
        self.outputFilename = "ExternalReactivityController.csv"
        # self.step_size = 0.1
        self.current_time = 0
        self.previous_time = self.current_time
        self.indexWrite = 0
        self.writeInterval = 0
        self.substep = 1
        self.maxGb = 2
        parametersToRecord: list[str] = None
        initialParametersBeforeFMUInit: dict = None

        # Store restart time step if memory leak
        self.restartPoints = []

        # Extract the FMU to a temporary directory
        self.unzipdir = fmpy.extract("ExternalReactivityController.fmu")

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

#=============================================================================*
