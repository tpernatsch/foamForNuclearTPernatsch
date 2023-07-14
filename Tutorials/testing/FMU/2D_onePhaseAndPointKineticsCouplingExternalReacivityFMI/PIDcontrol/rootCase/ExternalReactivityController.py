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

from pyfmi import load_fmu
import os
import pandas as pd
import numpy as np

#=============================================================================*


class ExternalReactivityController(FMU4FOAM.FMUBase):

    def __init__(self, endTime, jsonFile) -> None:
        super().__init__(endTime, jsonFile)

        # Some useful variables
        self.outputFilename = 'ExternalReactivityController.csv'
        self.step_size = 0.1
        self.current_time = 0

        # Create the FMU
        self.fmu1 = load_fmu('ExternalReactivityController.fmu')
        # self.fmu1.set_log_level(7)

        # Start initialization mode
        self.fmu1.initialize()

        # Extract all the parameters and variables composing the FMU
        self.parameterList = [key for key in self.fmu1.get_model_variables()]

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
        isFirstWrite = pd.isnull(self.results.index.max())
        self.results.loc[0] = [self.current_time] + [self.getVar(key) for key in self.parameterList]
        if (isFirstWrite):
            self.results.to_csv(self.outputFilename, index=False)
        else:
            self.results.to_csv(self.outputFilename, index=False, mode='a', header=False)


    def __del__(self) -> None:
        pass

#=============================================================================*
