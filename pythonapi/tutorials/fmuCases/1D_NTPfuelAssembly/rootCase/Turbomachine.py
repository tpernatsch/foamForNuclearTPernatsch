"""

"""
#=============================================================================*
# Imports

from foamForNuclear.fmi import OMSimulatorContainer, FMPyContainer

#=============================================================================*

# class Turbomachine(PyFMIContainer):
#     def __init__(
#             self,
#             endTime,
#             jsonFile
#         ):
#         super().__init__(
#             endTime,
#             jsonFile,
#             fmuName='modelica/NTPTurbomachine2.fmu',
#             outputFilename='Turbomachine.csv',
#             log_level=0,
#             writeInterval=0.05,
#             parametersToRecord=[
#                 'pTurbine_in',
#                 'SinkP1.in_p0',
#                 'pNozzleChamber_out',
#                 'sensTnozzleWall_in.T',
#                 'sensTnozzleWall_out.T',
#                 'Tturbine_out',
#                 'Tnozzle_in',
#                 'sourceMassFlow.in_T',
#                 'NozzleChamber.T[1]',
#                 'NozzleChamber.T[11]',
#                 'sensW.w',
#                 'mFlowTurbine_out',
#                 'mFlowCoreOutlet_in',
#                 'controlSystem.massFlowCmd.y',
#                 'controlSystem.extReactivity',
#                 'corePower_in',
#                 'controlSystem.powerCmd.y',
#                 'controlSystem.corePower.y',
#                 'Turbine1.omega',
#                 'valveLin.cmd'
#             ]
#         )


# class Turbomachine(OMSimulatorContainer):
#     def __init__(
#             self,
#             endTime,
#             jsonFile
#         ):
#         super().__init__(
#             endTime,
#             jsonFile,
#             fmuName='modelica/NTPTurbomachine2.fmu',
#             outputFilename='Turbomachine.csv',
#         )


class Turbomachine(FMPyContainer):
    def __init__(
            self,
            endTime,
            jsonFile
        ):
        super().__init__(
            endTime,
            jsonFile,
            fmuName='../modelica/NTPTurbomachine.fmu',
            outputFilename='Turbomachine.csv',
            writeInterval=0.1,
            substep=1,
            startTime=10,
            parametersToRecord=[
                'pTurbine_in',
                'SinkP1.in_p0',
                'pNozzleChamber_out',
                'sensTnozzleWall_in.T',
                'sensTnozzleWall_out.T',
                'Tturbine_out',
                'Tnozzle_in',
                'sourceMassFlow.in_T',
                'NozzleChamber.T[1]',
                'NozzleChamber.T[11]',
                'sensW.w',
                'corePower_in',
                'mFlowTurbine_out',
                'mFlowCoreOutlet_in',
                'Turbine1.omega',
                'valveLin.cmd',
                'controlSystem.PI_massFlow.u_s',
                'controlSystem.PI_massFlow.u_m',
                'controlSystem.extReactivity',
                'controlSystem.controlDrumReactivity.angle',
                'controlSystem.PI_Power.u_s',
                'controlSystem.PI_Power.u_m',
                'controlSystem.corePower.y',
                'controlSystem.PI_Power.y',
                'controlSystem.PI_DoublingTime.y',
                'controlSystem.PI_DoublingTime.u_m',
                # 'controlSystem.minCDangle.y',
                'controlSystem.PI_Thrust.y',
                'controlSystem.PI_Thrust.u_s',
                'controlSystem.PI_Thrust.u_m',
                'controlSystem.PI_DeltaV.u_s',
                'controlSystem.PI_DeltaV.u_m',
                'controlSystem.thrustCalculator.thrust',
                'controlSystem.thrustCalculator.exhaustVelocity',
                'controlSystem.doublingTimeCalculator.doublingTime',
                'controlSystem.thrustCmd.y',
                'controlSystem.deltaVCalculator.deltaV',
                'controlSystem.deltaVCalculator.endPropellant.activePort',
                'controlSystem.deltaVCalculator.totalMassPropellant.y',
                'controlSystem.deltaVCalculator.totalMassPropellant.u1',
                'controlSystem.deltaVCalculator.totalMassPropellant.u2',
                # 'controlSystem.feedbackThrust.u1',
            ]
        )

#=============================================================================*
