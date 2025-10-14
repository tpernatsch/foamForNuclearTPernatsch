"""

"""
#=============================================================================*
# Imports

from foamForNuclear.fmi import OMSimulatorContainer, FMPyContainer

#=============================================================================*

# class ExternalReactivityController(PyFMIContainer):
# class ExternalReactivityController(OMSimulatorContainer):
class ExternalReactivityController(FMPyContainer):
    def __init__(
            self,
            endTime,
            jsonFile
        ):
        super().__init__(
            endTime,
            jsonFile,
            fmuName='ExternalReactivityController.fmu',
            outputFilename='ExternalReactivityController.csv'
        )


#=============================================================================*
