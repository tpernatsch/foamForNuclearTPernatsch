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
import os
import pandas as pd
import numpy as np
import math

#=============================================================================*


class secondaryCircuit(FMU4FOAM.FMUBase):

    def __init__(self, endTime: float, jsonFile: str) -> None:
        print("Start Secondary Circuit Modelica construction", flush=True, end='\n')
        jsonFile = "FMU.json"
        super().__init__(endTime, jsonFile)

        # Some useful variables
        self.outputFilename: str = 'secondaryCircuit.csv'
        self.step_size: float = 1
        self.current_time: float = 0
        self.writerCounter: int = 0

        # Create the FMU
        self.fmu = load_fmu('modelica/SecondaryCircuitWithControl.fmu')
        # self.fmu.set_log_level(7) # Print useful information in case of bugs
        print("Finish loading FMU", flush=True, end='\n')

        # Start initialization mode
        self.fmu.initialize()
        print("Finish initialize FMU", flush=True, end='\n')

        # Extract all the parameters and variables composing the FMU
        self.parameterList = [key for key in self.fmu.get_model_variables()]

        # Create the output buffer (write into a csv)
        self.results = pd.DataFrame(columns=['time'] + self.parameterList)

        # If a simulation has already been performed -> restart
        if (os.path.exists(self.outputFilename)):
            # Read the last line of the output file
            data = pd.read_csv(self.outputFilename)
            lastIdx = data.index.max()

            # Set all the variables for the restart based on the previous 
            # simulation
            for d in data.columns:
                val = np.float64(data[d][lastIdx])
                if (d != "time"):
                    self.setVar(d, val)
                else:
                    self.current_time = val

                print(f"  {d}", flush=True, end='\n')

            # Reset buffer
            self.results.loc[0] = data.loc[lastIdx]

        # Set parameters relevant for load follow. Be carefull when changing
        # the name isLoadFollow as it is used to set the load follow transient
        isLoadFollowReactivity = False
        if (isLoadFollowReactivity):
            self.setVar("stepActivateExtReactivity.startTime", 1000)

        isLoadFollowSource = False
        if (isLoadFollowSource):
            self.setVar("stepActivateExtSource.startTime", 1000)

        if (isLoadFollowReactivity or isLoadFollowSource):
            self.setVar("rampLoad.startTime", 1100)
            self.setVar("rampLoad.height", -20e6)
            self.setVar("rampLoad2.startTime", 1600)
            self.setVar("rampLoad2.height", 20e6)

        

    def setVar(self, key: str, val: float) -> None:
        valueRef = self.getValueRef(key)
        self.fmu.set_real([valueRef], [val])


    def getVar(self, key: str) -> float:
        valueRef = self.getValueRef(key)
        return(self.fmu.get_real([valueRef])[0])
    

    def getValueRef(self, key: str) -> int:
        return(self.fmu.get_variable_valueref(key))


    def stepUntil(self, t: float) -> None:
        # Compute the time step for the FMU
        step_size = (t - self.current_time) / 1.0
        while (self.current_time < t):
            # Make a step with the FMU
            self.fmu.do_step(
                current_t=self.current_time,
                step_size=step_size,
                new_step=True
            )
            self.current_time += step_size

        # Update the current time
        self.current_time = round(t, 8)
        
        # Write the results
        isFirstWrite = pd.isnull(self.results.index.max())
        self.results.loc[0] = [self.current_time] + [self.getVar(key) for key in self.parameterList]
        if (isFirstWrite):
            self.results.to_csv(self.outputFilename, index=False)
        else:
            # Temporary fix to limit memory leakage because the output csv 
            # file is held in the RAM, it writes every seconds of simulation
            if (math.floor(self.current_time) == self.current_time):
            #     self.writerCounter += 1
            
            # if (self.writerCounter >= 10):
                self.results.to_csv(
                    self.outputFilename, index=False, mode='a', header=False
                )
                # self.writerCounter = 0


    def __del__(self) -> None:
        pass

#=============================================================================*
