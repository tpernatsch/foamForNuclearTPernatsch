"""

"""
#=============================================================================*
# Imports

from foamForNuclear.fmi import OMSimulatorContainer, FMPyContainer

#=============================================================================*

# class ExternalReactivityController(PyFMIContainer):
# class ExternalReactivityController(OMSimulatorContainer):
class Controller(FMPyContainer):
    def __init__(
            self,
            endTime,
            jsonFile
        ):
        super().__init__(
            endTime,
            jsonFile,
            fmuName='Controller.fmu',
            outputFilename='Controller.csv',
            parametersToRecord=[
                'massFlow',
                # 'inletP',
                # 'outletP',
                'massFlowCmd.y',
                'pid.u_s',
                'pid.u_m',
                'pumpMomentum'
            ],
            writeInterval=0.1
        )


#=============================================================================*
