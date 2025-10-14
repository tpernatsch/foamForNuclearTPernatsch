from pythonfmu.variables import Boolean
from FMU4FOAM import Fmi2Causality, Fmi2Variability, Real, String
from FMU4FOAM import OF2Fmu


class CoreFMU(OF2Fmu):

    author = "Thomas Guilbaud"
    description = "Core model."

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        print("Init Core FMU", flush=True)

        self.gfExtReact: float = 0
        self.gfPower: float = 0

        # Inputs
        self.register_variable(Real("gfExtReact",
            causality=Fmi2Causality.input, variability=Fmi2Variability.continuous))

        # Outputs
        self.register_variable(Real("gfPower",
            causality=Fmi2Causality.output, variability=Fmi2Variability.continuous))

        # General parameters
        self.outputPath = "CoreFMU" # will extract the OfCase files in self.outputPath

        # However in this test, self.port will be overwrite by pyfmi
        self.host = "\"127.0.0.1\""
        self.port = 8003
        self.oscmd = "bash"
        self.startTime = 0
        self.writeInterval = 1
        self.timeStep = 0.02
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
        # self.register_variable(Real("oscmd",
        #     causality=Fmi2Causality.parameter, variability=Fmi2Variability.tunable))

        # Recover delta time from solver
        self.t = self.startTime
        self.register_variable(Real("t", causality=Fmi2Causality.output))
        self.dt = self.timeStep
        self.register_variable(Real("dt", causality=Fmi2Causality.output))
        self.isConverged = 1
        self.register_variable(Real("isConverged", causality=Fmi2Causality.output))
