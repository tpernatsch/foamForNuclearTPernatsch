from pythonfmu.variables import Boolean
from FMU4FOAM import Fmi2Causality, Fmi2Variability, Real, Boolean , Integer, String
from FMU4FOAM import OF2Fmu
import numpy as np


class FuelPinFMU(OF2Fmu):

    author = "Thomas Guilbaud, EPFL/Transmutex, 2023/05/03"
    description = "Fuel pin(s) model for multi-scale nuclear reactor simulations."

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        print("Init OFFBEAT fuel pin FMU", flush=True)

        nPoints = 15 # 63 # 43

        self.offbeatAxialLocations: str = ""
        self.offbeatQavg: float = 0 # 230e2
        self.offbeatHeatSrcAxial: str = ""
        # self.offbeatTstructAxial: str = ""
        self.offbeatTfluidAxial: str = ""
        self.offbeatHtcAxial: str = ""
        self.offbeatSurfacePowerAxial: str = ""
        self.offbeatRhoCpdTdtAxial: str = ""
        self.offbeatTFuelAxial: str = ""
        self.offbeatTCladAxial: str = ""

        for z in np.linspace(-1.06, 0.55+1.02+0.12 - 1.060, nPoints): # np.linspace(0, 1.847, nPoints):
            self.offbeatAxialLocations += f"{z:.6f} "
            self.offbeatHeatSrcAxial += "1 "
            # self.offbeatTstructAxial += "673.15 "
            self.offbeatTfluidAxial += "673.15 "
            self.offbeatHtcAxial += "0.0 "
            self.offbeatSurfacePowerAxial += "0 "
            self.offbeatRhoCpdTdtAxial += "0 "
            self.offbeatTFuelAxial += "673.15 "
            self.offbeatTCladAxial += "673.15 "

        self.register_variable(String("offbeatAxialLocations",
            causality=Fmi2Causality.input, variability=Fmi2Variability.discrete))
        self.register_variable(Real("offbeatQavg",
            causality=Fmi2Causality.input, variability=Fmi2Variability.continuous))
        self.register_variable(String("offbeatHeatSrcAxial",
            causality=Fmi2Causality.input, variability=Fmi2Variability.discrete))
        # self.register_variable(String("offbeatTstructAxial",
        #     causality=Fmi2Causality.input, variability=Fmi2Variability.discrete))
        self.register_variable(String("offbeatTfluidAxial",
            causality=Fmi2Causality.input, variability=Fmi2Variability.discrete))
        self.register_variable(String("offbeatHtcAxial",
            causality=Fmi2Causality.input, variability=Fmi2Variability.discrete))
        self.register_variable(String("offbeatSurfacePowerAxial",
            causality=Fmi2Causality.output, variability=Fmi2Variability.discrete))
        self.register_variable(String("offbeatTFuelAxial",
            causality=Fmi2Causality.output, variability=Fmi2Variability.discrete))
        self.register_variable(String("offbeatTCladAxial",
            causality=Fmi2Causality.output, variability=Fmi2Variability.discrete))
        self.register_variable(String("offbeatRhoCpdTdtAxial",
            causality=Fmi2Causality.output, variability=Fmi2Variability.discrete))

        # General parameters
        self.outputPath = "FuelPin" # will extract the OfCase files in self.outputPath

        # However in this test, self.port will be overwrite by pyfmi
        self.host = "\"127.0.0.1\""
        self.port = 8003
        self.oscmd = 'bash'
        self.startTime = 0
        self.writeInterval = 1
        self.timeStep = 0.1
        self.STDLAP = "Gauss linear corrected"
        self.STDGRAD = "pointCellsLeastSquares"

        # self.register_variable(String("host",
        #     causality=Fmi2Causality.parameter, variability=Fmi2Variability.tunable))
        self.register_variable(String("STDGRAD",
            causality=Fmi2Causality.parameter, variability=Fmi2Variability.tunable))
        self.register_variable(String("STDLAP",
            causality=Fmi2Causality.parameter, variability=Fmi2Variability.tunable))
        self.register_variable(Real("startTime",
            causality=Fmi2Causality.parameter, variability=Fmi2Variability.tunable))
        self.register_variable(Real("writeInterval",
            causality=Fmi2Causality.parameter, variability=Fmi2Variability.tunable))
        self.register_variable(Real("timeStep",
            causality=Fmi2Causality.parameter, variability=Fmi2Variability.tunable))

        # Recover delta time from solver
        self.t = self.startTime
        self.register_variable(Real("t", causality=Fmi2Causality.output))
        self.dt = self.timeStep
        self.register_variable(Real("dt", causality=Fmi2Causality.output))
        self.isConverged = 1
        self.register_variable(Real("isConverged", causality=Fmi2Causality.output))
