"""
Python API to manage the FMU in GeN-Foam. Use PyFMI (v2.10.0) instead of
OMSimulator to be able to reset the FMU state.
Inspired by https://github.com/DLR-RY/FMU4FOAM/tree/master/examples.

Author: Thomas Guilbaud, EPFL/Transmutex SA, 11/2022
"""

#=============================================================================*
# Imports
#=============================================================================*

import FMU4FOAM

from pyfmi import load_fmu
from pyfmi.master import Master
import os
import pandas as pd
import numpy as np
import math

#=============================================================================*


class secondaryCircuit(FMU4FOAM.FMUBase):

    def __init__(self, endTime, jsonFile) -> None:
        jsonFile = "./FMU.json"
        super().__init__(endTime, jsonFile)

        # Some useful variables
        self.outputFilenameSecCir = 'secondaryCircuit.csv'
        self.outputFilenamePIDs = 'pidControllers.csv'
        self.step_size = 0.1
        self.current_time = 0
        self.isRestart = os.path.exists(self.outputFilenameSecCir) and os.path.exists(self.outputFilenamePIDs)

        print("Restart simulation" if self.isRestart else "First simulation", flush=True)

        # Create the FMUs
        logLevel = 2
        self.secondaryCircuitFMU = load_fmu('modelica/SecondaryCircuitWithoutControl.fmu', 'cs', log_level=logLevel)
        self.pidControllersFMU = load_fmu('simulink/pidControllers.fmu', 'cs', log_level=logLevel)
        # self.pidControllersFMU = load_fmu('modelica/PIDfreq.fmu', log_level=logLevel) # Modelica model for test, should be equivalent to Simulink
        
        models = [
            self.secondaryCircuitFMU, self.pidControllersFMU
        ]

        # Connections from output to input between FMUs
        connections = [
            (self.secondaryCircuitFMU, "frequencyOutput", self.pidControllersFMU, "measFrequency"),
            (self.secondaryCircuitFMU, "electricalPowerOutput", self.pidControllersFMU, "measElectricalPower"),
            (self.pidControllersFMU, "pidFrequencyResponse", self.secondaryCircuitFMU, "valveControlInput")
        ]

        # Read the last line of the output file to get the latest time if there
        # was a restart
        if (self.isRestart):
            data = pd.read_csv(self.outputFilenameSecCir)
            self.current_time = np.float64(data["time"][data.index.max()])

        # Start initialization models and set start time 
        for model in models:
            model.initialize()
            if ("time" in [key for key in model.get_model_variables()]):
                model.set("time", self.current_time)

        # Create the Master simulator to handle all the FMUs
        self.master_simulator = Master(models, connections)
        opts = self.master_simulator.simulate_options()
        opts["step_size"] = self.step_size
        opts["result_handling"] = "file"
        # Don't intialize if restart the simulation, initialization will be
        # done below
        opts['initialize'] = not self.isRestart
        # Add the options to the master to enable the FMUs communication
        self.master_simulator.opts = opts

            
        # Extract all the parameters and variables composing the FMUs
        self.paramsSecondaryCircuit = [key for key in self.secondaryCircuitFMU.get_model_variables()]
        self.paramsPidControllers = [key for key in self.pidControllersFMU.get_model_variables()]

        # Create the output buffer (write into a csv) if first simulation
        # if (not self.isRestart):
        self.resultsSecondaryCircuit = pd.DataFrame(columns=['time'] + self.paramsSecondaryCircuit)
        self.resultsPidControllers = pd.DataFrame(columns=['time'] + self.paramsPidControllers)

        # If a simulation has already been performed -> restart
        if (self.isRestart):
            # Read the last line of the output file
            data = pd.read_csv(self.outputFilenameSecCir)
            lastIdx = data.index.max()

            # Set all the variables for the restart based on the previous 
            # simulation
            for d in data.columns:
                val = np.float64(data[d][lastIdx])
                if (d != "time"):
                    self.setVar(d, val)

            # Reset buffer
            self.resultsSecondaryCircuit.loc[0] = data.loc[lastIdx]

        # If a simulation has already been performed -> restart
        if (self.isRestart):
            # Read the last line of the output file
            data = pd.read_csv(self.outputFilenamePIDs)
            lastIdx = data.index.max()

            # Set all the variables for the restart based on the previous 
            # simulation
            self.pidControllersFMU.set("measFrequency", data["measFrequency"][lastIdx])
            self.pidControllersFMU.set("measElectricalPower", data["measElectricalPower"][lastIdx])
            self.pidControllersFMU.set("innerCorePower", data["innerCorePower"][lastIdx])
            self.pidControllersFMU.set("outerCorePower", data["outerCorePower"][lastIdx])
            self.pidControllersFMU.set("pidFrequencyResponse", data["pidFrequencyResponse"][lastIdx])

            # Force reset of PID frequency (https://ch.mathworks.com/help/simulink/slref/pidcontroller.html)
            self.pidControllersFMU.set("integralInitValue", data["integralOutput"][lastIdx]) # P=0.4; I=4
            self.pidControllersFMU.set("resetPidFrequency", 1)

            # Reset buffer
            self.resultsPidControllers.loc[0] = data.loc[lastIdx]

        else:
            # Force initialization of PID to avoid error
            self.pidControllersFMU.set("measFrequency", 50)


        # Set parameters relevant for load follow. Be carefull when changing
        # the name isLoadFollow as it is used to set the load follow transient
        isLoadFollowSource = False
        if (isLoadFollowSource):
            self.setVar("stepActivateExtSource.startTime", 1000)

        isLoadFollowReactivity = False
        if (isLoadFollowReactivity or isLoadFollowSource):
            self.setVar("rampLoad.startTime", 1100)
            self.setVar("rampLoad.height", -20e6)
            self.setVar("rampLoad2.startTime", 1600)
            self.setVar("rampLoad2.height", 20e6)


    def setVar(self, key: str, val: float) -> None:
        valueRef = self.getValueRef(key)
        if (key in self.paramsSecondaryCircuit):
            self.secondaryCircuitFMU.set_real([valueRef], [val])
            # self.secondaryCircuitFMU.set(key, val)
        elif (key in self.paramsPidControllers):
            self.pidControllersFMU.set_real([valueRef], [val])


    def getVar(self, key: str) -> float:
        valueRef = self.getValueRef(key)
        if (key in self.paramsSecondaryCircuit):
            return(self.secondaryCircuitFMU.get_real([valueRef])[0])
            # return(self.secondaryCircuitFMU.get(key)[0])
        elif (key in self.paramsPidControllers):
            return(self.pidControllersFMU.get_real([valueRef])[0])
    

    def getValueRef(self, key: str) -> int:
        if (key in self.paramsSecondaryCircuit):
            return(self.secondaryCircuitFMU.get_variable_valueref(key))
        elif (key in self.paramsPidControllers):
            return(self.pidControllersFMU.get_variable_valueref(key))
            

    def stepUntil(self, t: float) -> None:
        # Compute the time step for the FMU
        step_size = (t - self.current_time) / 2.0 # 2 improves the statbility

        # Loop inspired by https://github.com/modelon-community/PyFMI/blob/master/src/pyfmi/master.pyx
        while (self.current_time < t):
            # Make a step with each FMU
            for model in self.master_simulator.models:
                model.do_step(self.current_time, step_size, new_step=True)

            # Exchange the information between FMUs
            self.master_simulator.exchange_connection_data()

            self.current_time += step_size

        # Update the current time
        self.current_time = round(t, 8)
        
        # Write the results
        isFirstWrite = pd.isnull(self.resultsSecondaryCircuit.index.max())
        self.resultsSecondaryCircuit.loc[0] = [self.current_time] + [self.getVar(key) for key in self.paramsSecondaryCircuit]
        self.resultsPidControllers.loc[0] = [self.current_time] + [self.pidControllersFMU.get(key)[0] for key in self.paramsPidControllers]
        if (isFirstWrite):
            self.resultsSecondaryCircuit.to_csv(self.outputFilenameSecCir, index=False)
            self.resultsPidControllers.to_csv(self.outputFilenamePIDs, index=False)
        else:
            # Temporary fix to limit memory leakage because the output csv 
            # file is held in the RAM, it writes every seconds of simulation
            if (math.floor(self.current_time) == self.current_time):
                self.resultsSecondaryCircuit.to_csv(
                    self.outputFilenameSecCir, 
                    index=False, mode='a', header=False
                )
                self.resultsPidControllers.to_csv(
                    self.outputFilenamePIDs, 
                    index=False, mode='a', header=False
                )

    def __del__(self) -> None:
        pass

#=============================================================================*
