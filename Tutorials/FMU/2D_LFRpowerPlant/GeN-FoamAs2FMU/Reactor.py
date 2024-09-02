from pythonfmu.variables import Boolean
from FMU4FOAM import Fmi2Causality, Fmi2Variability, Real, Boolean , Integer, String
from FMU4FOAM import OF2Fmu


class Reactor(OF2Fmu):

    author = "Thomas Guilbaud, EPFL/Transmutex, 2023/05/09"
    description = "FMU class to handle the reactor modeled in GeN-Foam."

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.gfExternalSourceMod = 0
        self.gfExternalReactivity = 0

        self.gfSteamGenTemperature1 = 631
        self.gfSteamGenTemperature2 = 634
        self.gfSteamGenTemperature3 = 634
        self.gfSteamGenTemperature4 = 634
        self.gfSteamGenTemperature5 = 634
        self.gfSteamGenTemperature6 = 649
        self.gfSteamGenTemperature7 = 699
        self.gfSteamGenTemperature8 = 722

        self.gfInnerCorePowerOut = 6487
        self.gfOuterCorePowerOut = 13755

        self.gfSteamGenPowerOut1 = 9080
        self.gfSteamGenPowerOut2 = 7510
        self.gfSteamGenPowerOut3 = 4230
        self.gfSteamGenPowerOut4 = 6110
        self.gfSteamGenPowerOut5 = 12130
        self.gfSteamGenPowerOut6 = 15640
        self.gfSteamGenPowerOut7 = 14100
        self.gfSteamGenPowerOut8 = 4660

        self.outputPath = "Reactor" # will extract the OfCase files in Testing
        # However in this test it will be overwrite by pyfmi

        self.register_variable(Real("gfExternalSourceMod", causality=Fmi2Causality.input))
        self.register_variable(Real("gfExternalReactivity", causality=Fmi2Causality.input))
        self.register_variable(Real("gfSteamGenTemperature1", causality=Fmi2Causality.input))
        self.register_variable(Real("gfSteamGenTemperature2", causality=Fmi2Causality.input))
        self.register_variable(Real("gfSteamGenTemperature3", causality=Fmi2Causality.input))
        self.register_variable(Real("gfSteamGenTemperature4", causality=Fmi2Causality.input))
        self.register_variable(Real("gfSteamGenTemperature5", causality=Fmi2Causality.input))
        self.register_variable(Real("gfSteamGenTemperature6", causality=Fmi2Causality.input))
        self.register_variable(Real("gfSteamGenTemperature7", causality=Fmi2Causality.input))
        self.register_variable(Real("gfSteamGenTemperature8", causality=Fmi2Causality.input))
        self.register_variable(Real("gfInnerCorePowerOut", causality=Fmi2Causality.output))
        self.register_variable(Real("gfOuterCorePowerOut", causality=Fmi2Causality.output))
        self.register_variable(Real("gfSteamGenPowerOut1", causality=Fmi2Causality.output))
        self.register_variable(Real("gfSteamGenPowerOut2", causality=Fmi2Causality.output))
        self.register_variable(Real("gfSteamGenPowerOut3", causality=Fmi2Causality.output))
        self.register_variable(Real("gfSteamGenPowerOut4", causality=Fmi2Causality.output))
        self.register_variable(Real("gfSteamGenPowerOut5", causality=Fmi2Causality.output))
        self.register_variable(Real("gfSteamGenPowerOut6", causality=Fmi2Causality.output))
        self.register_variable(Real("gfSteamGenPowerOut7", causality=Fmi2Causality.output))
        self.register_variable(Real("gfSteamGenPowerOut8", causality=Fmi2Causality.output))

        self.writeInterval = 1.0
        self.deltaT = 1.0
        self.STDLAP = "Gauss linear corrected"
        self.STDGRAD = "pointCellsLeastSquares"
        

        self.register_variable(String("STDGRAD", 
            causality=Fmi2Causality.parameter, variability=Fmi2Variability.tunable))
        self.register_variable(String("STDLAP", 
            causality=Fmi2Causality.parameter, variability=Fmi2Variability.tunable))
        self.register_variable(Real("writeInterval", 
            causality=Fmi2Causality.parameter, variability=Fmi2Variability.tunable))
        self.register_variable(Real("deltaT", 
            causality=Fmi2Causality.parameter, variability=Fmi2Variability.tunable))
        