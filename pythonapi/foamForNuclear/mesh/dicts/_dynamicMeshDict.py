

from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import List, OpenFOAMDict, Table, Vector, addParameter, tab
from foamForNuclear.openfoamFile import OpenFOAMFile

_DYNAMIC_FV_MESH_TYPES = {
    "staticFvMesh", "solidBodyMotionFvMesh", "dynamicRefineFvMesh",
    "dynamicMotionSolverFvMesh"
}
_MOTION_SOLVER_TYPES = {
    "solidBody",
    "coded", "displacementInterpolation", "displacementLayeredMotion",
    "displacementPointSmoothing", "displacementSmartPointSmoothing",
    "hexMeshSmoother", "multiDisplacement", "multiSolidBodyMotionSolver",
    "velocityDisplacement",

    "displacementSBRStress", "displacementLaplacian",
    "velocityLaplacian",
}
_SOLID_BODY_MOTION_FUNCTION_TYPES = {
    "axisRotationMotion", "linearMotion", "multiMotion",
    "oscillatingLinearMotion", "oscillatingRotatingMotion", "rotatingMotion",
    "SDA", "solidBodyMotionFunction", "tabulated6DoFMotion"
}
_DIFFUSIVITY_MOTION_TYPES = {
    "uniform", "directional", "motionDirectional", "quadratic", "file"
}
_QUADRATIC_MOTION_TYPES = {"inverseDistance"}
_INTERPOLATION_SCHEME_TYPES = {"spline", "linear"}


class MotionDict(OpenFOAMDict):
    def __init__(
            self,
            name: str,
            solidBodyMotionFunction: str
        ):
        super().__init__(name=name)

        self.solidBodyMotionFunction = solidBodyMotionFunction


    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name) -> None:
        check_type("name", name, str)
        self._name = name


    @property
    def solidBodyMotionFunction(self):
        return self._solidBodyMotionFunction

    @solidBodyMotionFunction.setter
    def solidBodyMotionFunction(self, solidBodyMotionFunction) -> None:
        check_type("solidBodyMotionFunction", solidBodyMotionFunction, str, none_ok=True)
        if (solidBodyMotionFunction is not None):
            check_value("solidBodyMotionFunction", solidBodyMotionFunction, _SOLID_BODY_MOTION_FUNCTION_TYPES)
        self._solidBodyMotionFunction = solidBodyMotionFunction


class LinearMotion(MotionDict):
    """
    Linear Motion used in dynamicMeshDict. The motion is defined by
    the following formula:

        displacement = velocity * t

    Parameters
    ----------
    name : str
        Name of the motion.
    velocity : Vector
        Motion velocity in m/s.
    """
    def __init__(
            self,
            name: str,
            velocity: Vector,
        ):
        super().__init__(name, solidBodyMotionFunction="linearMotion")

        self.velocity = velocity

    def __repr__(self, depth = 0):
        linearMotionCoeffs = OpenFOAMDict(
            items={
                'velocity': self.velocity,
            }
        )
        self.__setitem__("solidBodyMotionFunction", self.solidBodyMotionFunction)
        self.__setitem__("linearMotionCoeffs", linearMotionCoeffs)

        return super().__repr__(depth)

    @property
    def description(self):
        return(f"{self.name}: {self.solidBodyMotionFunction} (velocity={self.velocity} m/s)")

    @property
    def velocity(self):
        return self._velocity

    @velocity.setter
    def velocity(self, velocity) -> None:
        check_type("velocity", velocity, Vector)
        self._velocity = velocity


class OscillatingLinearMotion(MotionDict):
    """
    Oscillating Linear Motion used in dynamicMeshDict. The motion is defined by
    the following formula:

        displacement(t) = amplitude * sin(omega * (t + timeShift)) + amplitudeShift

    Parameters
    ----------
    name : str
        Name of the motion.
    amplitude : Vector | Table
        Amplitude vector of the oscilation in m.
    omega : float | Table
        Angular frequency in rad/s.
    timeShift : float | Table
        Phase shift in s (default `None`).
    amplitudeShift : Vector | Table
        Vertical shift in m (default `None`).
    """
    def __init__(
            self,
            name: str,
            amplitude: Vector | Table,
            omega: float | Table,
            timeShift: float | Table=None,
            amplitudeShift: Vector | Table=None
        ):
        super().__init__(name, solidBodyMotionFunction="oscillatingLinearMotion")

        self.amplitude = amplitude
        self.omega = omega
        self.timeShift = timeShift
        self.amplitudeShift = amplitudeShift

    def __repr__(self, depth = 0):
        oscillatingLinearMotionCoeffs = OpenFOAMDict(
            items={
                'amplitude': self.amplitude,
                'omega': self.omega
            }
        )
        if (self.timeShift is not None):
            oscillatingLinearMotionCoeffs['timeShift'] = self.timeShift
        if (self.amplitudeShift is not None):
            oscillatingLinearMotionCoeffs['amplitudeShift'] = self.amplitudeShift

        self.__setitem__("solidBodyMotionFunction", self.solidBodyMotionFunction)
        self.__setitem__("oscillatingLinearMotionCoeffs", oscillatingLinearMotionCoeffs)

        return super().__repr__(depth)

    @property
    def description(self):
        if (isinstance(self.amplitude, Table) and isinstance(self.omega, Table)):
            return(f"{self.name}: {self.solidBodyMotionFunction} (amplitude=Table m; omega=Table rad/s)")
        elif (isinstance(self.amplitude, Table)):
            return(f"{self.name}: {self.solidBodyMotionFunction} (amplitude=Table m; omega={self.omega} rad/s)")
        elif (isinstance(self.omega, Table)):
            return(f"{self.name}: {self.solidBodyMotionFunction} (amplitude={self.amplitude} m; omega=Table rad/s)")

        return(f"{self.name}: {self.solidBodyMotionFunction} (amplitude={self.amplitude} m; omega={self.omega} rad/s)")

    @property
    def amplitude(self):
        return self._amplitude

    @amplitude.setter
    def amplitude(self, amplitude) -> None:
        check_type("amplitude", amplitude, (Vector, Table))
        self._amplitude = amplitude

    @property
    def omega(self):
        return self._omega

    @omega.setter
    def omega(self, omega) -> None:
        check_type("omega", omega, (float, int, Table))
        self._omega = omega

    @property
    def timeShift(self):
        return self._timeShift

    @timeShift.setter
    def timeShift(self, timeShift) -> None:
        check_type("timeShift", timeShift, (float, int, Table), none_ok=True)
        self._timeShift = timeShift

    @property
    def amplitudeShift(self):
        return self._amplitudeShift

    @amplitudeShift.setter
    def amplitudeShift(self, amplitudeShift) -> None:
        check_type("amplitudeShift", amplitudeShift, (Vector, Table), none_ok=True)
        self._amplitudeShift = amplitudeShift


class RotatingMotion(MotionDict):
    """
    Rotating Motion used in dynamicMeshDict.

    Parameters
    ----------
    name : str
        Name of the motion.
    origin : Vector
        Origin position of the rotation in m (default `Vector(0, 0, 0)`).
    axis : Vector
        Axis rotation (default `Vector(0, 0, 1)`, along the Z-axis).
    omega : float | Table
        Angular frequency in rad/s.
    """
    def __init__(
            self,
            name: str,
            omega: float | Table,
            origin: Vector=Vector(0, 0, 0),
            axis: Vector=Vector(0, 0, 1),
        ):
        super().__init__(name, solidBodyMotionFunction="rotatingMotion")

        self.omega = omega
        self.origin = origin
        self.axis = axis

    def __repr__(self, depth = 0):
        rotatingMotionCoeffs = OpenFOAMDict(
            items={
                'origin': self.origin,
                'axis': self.axis,
                'omega': self.omega
            }
        )
        self.__setitem__("solidBodyMotionFunction", self.solidBodyMotionFunction)
        self.__setitem__("rotatingMotionCoeffs", rotatingMotionCoeffs)

        return super().__repr__(depth)

    @property
    def description(self):
        if (isinstance(self.omega, Table)):
            return(f"{self.name}: {self.solidBodyMotionFunction} (origin={self.origin} m; axis={self.axis}; omega=Table rad/s)")

        return(f"{self.name}: {self.solidBodyMotionFunction} (origin={self.origin} m; axis={self.axis}; omega={self.omega} rad/s)")

    @property
    def omega(self):
        return self._omega

    @omega.setter
    def omega(self, omega) -> None:
        check_type("omega", omega, (float, int, Table))
        self._omega = omega

    @property
    def origin(self):
        return self._origin

    @origin.setter
    def origin(self, origin) -> None:
        check_type("origin", origin, Vector)
        self._origin = origin

    @property
    def axis(self):
        return self._axis

    @axis.setter
    def axis(self, axis) -> None:
        check_type("axis", axis, Vector)
        self._axis = axis


class OscillatingRotatingMotion(MotionDict):
    """
    Oscillating Rotating Motion used in dynamicMeshDict. The motion is defined
    by the following formula:

        theta(t) = amplitude * sin(omega*t)

    Parameters
    ----------
    name : str
        Name of the motion.
    omega : float
        Angular frequency in rad/s.
    amplitude : Vector
        Amplitude in degrees.
    origin : Vector
        Origin position of the rotation in m (default `Vector(0, 0, 0)`).
    """
    def __init__(
            self,
            name: str,
            omega: float,
            amplitude: Vector,
            origin: Vector=Vector(0, 0, 0),
        ):
        super().__init__(name, solidBodyMotionFunction="oscillatingRotatingMotion")

        self.omega = omega
        self.origin = origin
        self.amplitude = amplitude

    def __repr__(self, depth = 0):
        oscillatingRotatingMotionCoeffs = OpenFOAMDict(
            items={
                'origin': self.origin,
                'amplitude': self.amplitude,
                'omega': self.omega
            }
        )
        self.__setitem__("solidBodyMotionFunction", self.solidBodyMotionFunction)
        self.__setitem__("oscillatingRotatingMotionCoeffs", oscillatingRotatingMotionCoeffs)

        return super().__repr__(depth)

    @property
    def description(self):
        return(f"{self.name}: {self.solidBodyMotionFunction} (origin={self.origin} m; amplitude={self.amplitude} deg; omega={self.omega} rad/s)")

    @property
    def omega(self):
        return self._omega

    @omega.setter
    def omega(self, omega) -> None:
        check_type("omega", omega, (float, int))
        self._omega = omega

    @property
    def origin(self):
        return self._origin

    @origin.setter
    def origin(self, origin) -> None:
        check_type("origin", origin, Vector)
        self._origin = origin

    @property
    def amplitude(self):
        return self._amplitude

    @amplitude.setter
    def amplitude(self, amplitude) -> None:
        check_type("amplitude", amplitude, Vector)
        self._amplitude = amplitude


class Tabulated6DoFMotion(MotionDict):
    """
    Tabulated 6 Degrees of Freedom (DoF) Motion used in dynamicMeshDict. The
    motion is defined by the table of time, position vector and rotation vector.

    Parameters
    ----------
    name : str
        Name of the motion.
    tableDoF : list
        Table containing the 6 DoF displacement.
    CofG : Vector
        Center of Gravity (default: `Vector(0, 0, 0)`).
    interpolationScheme : str, {'spline', 'linear'}
        Interpolation scheme between data points (default `linear`).
    """
    def __init__(
            self,
            name: str,
            tableDoF: list | None=None,
            CofG: Vector=Vector(0, 0, 0),
            interpolationScheme: str='linear'
        ):
        super().__init__(name, solidBodyMotionFunction="tabulated6DoFMotion")

        self.tableDoF = tableDoF
        self.CofG = CofG
        self.interpolationScheme = interpolationScheme

    def __repr__(self, depth: int=0, region: str=""):
        with open(f"constant/{region}/{self.name}.dat", 'w') as f:
            tabulated6DoFMotionTable = f"{len(self.tableDoF)}\n(\n"

            for t, position, rotation in self.tableDoF:
                tabulated6DoFMotionTable += f"{tab}( {t} ( {position} {rotation} ) )\n"

            tabulated6DoFMotionTable += ")\n"

            f.write(tabulated6DoFMotionTable)

        self.__setitem__("solidBodyMotionFunction", self.solidBodyMotionFunction)
        self.__setitem__("timeDataFileName", f'"<constant>/{region}/{self.name}.dat"')
        self.__setitem__("CofG", self.CofG)
        self.__setitem__("interpolationScheme", self.interpolationScheme)

        return super().__repr__(depth)

    @property
    def description(self):
        return(f"{self.name}: {self.solidBodyMotionFunction} ({len(self.tableDoF)} data points)")

    @property
    def tableDoF(self):
        return self._tableDoF

    @tableDoF.setter
    def tableDoF(self, tableDoF) -> None:
        check_type("tableDoF", tableDoF, list, none_ok=True)
        self._tableDoF = [] if tableDoF is None else tableDoF

    @property
    def CofG(self):
        return self._CofG

    @CofG.setter
    def CofG(self, CofG) -> None:
        check_type("CofG", CofG, Vector)
        self._CofG = CofG

    @property
    def interpolationScheme(self):
        return self._interpolationScheme

    @interpolationScheme.setter
    def interpolationScheme(self, interpolationScheme) -> None:
        check_type("interpolationScheme", interpolationScheme, str)
        check_value("interpolationScheme", interpolationScheme, _INTERPOLATION_SCHEME_TYPES)
        self._interpolationScheme = interpolationScheme

    def add_point(
            self,
            t: float,
            position: Vector=Vector(0, 0, 0),
            rotation: Vector=Vector(0, 0, 0)
        ):
        check_type("t", t, (float, int))
        check_type("position", position, Vector)
        check_type("rotation", rotation, Vector)
        self.tableDoF.append((t, position, rotation))


class DiffusivityMotion:
    """
    Diffusivity Motion used in `dynamicMeshDict`.

    Parameters
    ----------
    type : str
        Type of diffusivity (options: `uniform`, `directional`,
        `motionDirectional`, `quadratic`, `file`)
    """
    def __init__(
            self,
            type: str
        ):
        self.type = type

    def __repr__(self):
        return(f"{self.type}")

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _DIFFUSIVITY_MOTION_TYPES)
        self._type = type


class UniformDiffusivityMotion(DiffusivityMotion):
    def __init__(self):
        super().__init__(type="uniform")


class DirectionalDiffusivityMotion(DiffusivityMotion):
    def __init__(self, direction: Vector):
        super().__init__(type="directional")

        self.direction = direction

    def __repr__(self):
        text = super().__repr__()
        text += f" {self.direction}"
        return text

    @property
    def direction(self):
        return self._direction

    @direction.setter
    def direction(self, direction) -> None:
        check_type("direction", direction, Vector)
        self._direction = direction


class MotionDirectionalDiffusivityMotion(DiffusivityMotion):
    def __init__(self, direction: Vector):
        super().__init__(type="motionDirectional")

        self.direction = direction

    def __repr__(self):
        text = super().__repr__()
        text += f" {self.direction}"
        return text

    @property
    def direction(self):
        return self._direction

    @direction.setter
    def direction(self, direction) -> None:
        check_type("direction", direction, Vector)
        self._direction = direction


class QuadraticDiffusivityMotion(DiffusivityMotion):
    """
    Quadratic diffusivity motion

    Parameters
    ----------
    quadraticType : str
        Options: `inverseDistance`.
    patches: List | list[str]
        List of patch names (default `List()`).
    """
    def __init__(
            self,
            quadraticType: str,
            patches: List | list=[]
        ):
        super().__init__(type="quadratic")

        self.quadraticType = quadraticType
        self.patches: List = patches

    def __repr__(self):
        text = super().__repr__()
        text += f" {self.quadraticType} {self.patches}"
        return text

    @property
    def quadraticType(self):
        return self._quadraticType

    @quadraticType.setter
    def quadraticType(self, quadraticType) -> None:
        check_type("quadraticType", quadraticType, str)
        check_value("quadraticType", quadraticType, _QUADRATIC_MOTION_TYPES)
        self._quadraticType = quadraticType

    @property
    def patches(self):
        return self._patches

    @patches.setter
    def patches(self, patches) -> None:
        check_type("patches", patches, (list, List))
        if (isinstance(patches, list)):
            self._patches = List(patches)
        else:
            self._patches = patches

    def add_patch(self, patchName: str):
        check_type("patchName", patchName, str)
        self.patches.append(patchName)


class FileDiffusivityMotion(DiffusivityMotion):
    def __init__(self, filename: str):
        super().__init__(type="file")

        self.filename = filename

    def __repr__(self):
        text = super().__repr__()
        text += f" {self.filename}"
        return text

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, filename) -> None:
        check_type("filename", filename, str)
        self._filename = filename



class DynamicMeshDict(OpenFOAMFile):
    """
    Dictionary used for mesh motion.

    Parameters
    ----------
    dynamicFvMesh : str {"dynamicRefineFvMesh", "dynamicMotionSolverFvMesh"}
        Dynamic mesh type.
    motionSolver : str {"solidBody"}
        Motion solver type.
    solidBodyMotionFunction : str
        Options: `"axisRotationMotion"`, `"linearMotion"`, `"multiMotion"`,
        `"oscillatingLinearMotion"`, `"oscillatingRotatingMotion"`,
        `"rotatingMotion"`, `"SDA"`, `"solidBodyMotionFunction"`,
        `"tabulated6DoFMotion"`
    region : str
        Name of the region
    """
    def __init__(
            self,
            dynamicFvMesh: str,
            motionSolver: str=None,
            solidBodyMotionFunction: str=None,
            region: str="",
            motionSolverLibs: List=List(["fvMotionSolvers"]),
            extraParameters: dict | None=None
        ):
        super().__init__("dynamicMeshDict", folder="constant", region=region)

        self.dynamicFvMesh = dynamicFvMesh
        self.motionSolver = motionSolver
        self.solidBodyMotionFunction = solidBodyMotionFunction
        self.motionSolverLibs = motionSolverLibs
        self.extraParameters = extraParameters if extraParameters is not None else {}

        self.motions = []
        self.diffusivity = None

    def __repr__(self, depth: int=0):
        text = f"{depth*tab}dynamicFvMesh {self.dynamicFvMesh}\n"

        if (self.motionSolver is not None):
            text += f"{depth*tab}motionSolver {self.motionSolver}\n"

            if (self.solidBodyMotionFunction is not None):
                text += f"{depth*tab}solidBodyMotionFunction {self.solidBodyMotionFunction}\n"

                for motion in self.motions:
                    text += f"{(depth+1)*tab}{motion.description}\n"

        return text

    @property
    def dynamicFvMesh(self):
        return self._dynamicFvMesh

    @dynamicFvMesh.setter
    def dynamicFvMesh(self, dynamicFvMesh) -> None:
        check_type("dynamicFvMesh", dynamicFvMesh, str)
        check_value("dynamicFvMesh", dynamicFvMesh, _DYNAMIC_FV_MESH_TYPES)
        self._dynamicFvMesh = dynamicFvMesh


    @property
    def motionSolver(self):
        return self._motionSolver

    @motionSolver.setter
    def motionSolver(self, motionSolver) -> None:
        check_type("motionSolver", motionSolver, str, none_ok=True)
        if (motionSolver is not None):
            check_value("motionSolver", motionSolver, _MOTION_SOLVER_TYPES)
        self._motionSolver = motionSolver


    @property
    def solidBodyMotionFunction(self):
        return self._solidBodyMotionFunction

    @solidBodyMotionFunction.setter
    def solidBodyMotionFunction(self, solidBodyMotionFunction) -> None:
        check_type("solidBodyMotionFunction", solidBodyMotionFunction, str, none_ok=True)
        if (solidBodyMotionFunction is not None):
            check_value("solidBodyMotionFunction", solidBodyMotionFunction, _SOLID_BODY_MOTION_FUNCTION_TYPES)
        self._solidBodyMotionFunction = solidBodyMotionFunction


    @property
    def motionSolverLibs(self):
        return self._motionSolverLibs

    @motionSolverLibs.setter
    def motionSolverLibs(self, motionSolverLibs) -> None:
        check_type("motionSolverLibs", motionSolverLibs, (list, List))
        if (isinstance(motionSolverLibs, list)):
            self._motionSolverLibs = List(motionSolverLibs)
        else:
            self._motionSolverLibs = motionSolverLibs


    @property
    def diffusivity(self):
        return self._diffusivity

    @diffusivity.setter
    def diffusivity(self, diffusivity) -> None:
        check_type("diffusivity", diffusivity, DiffusivityMotion, none_ok=True)
        self._diffusivity = diffusivity

    @property
    def is_empty(self) -> bool:
        return(self.motionSolver == None)


    def add_motion(self, motion: MotionDict):
        check_type("motion", motion, MotionDict)
        self.motions.append(motion)

    def add_rotating_motion(
            self,
            name: str,
            omega: float,
            origin: Vector=Vector(0, 0, 0),
            axis: Vector=Vector(0, 0, 1),
        ) -> None:
        """
        Add a Rotating Motion.

        Parameters
        ----------
        name : str
            Name of the motion.
        origin : Vector
            Origin position of the rotation in m (default `Vector(0, 0, 0)`).
        axis : Vector
            Axis rotation (default `Vector(0, 0, 1)`, along the Z-axis).
        omega : float
            Angular frequency in rad/s.
        """
        self.add_motion(RotatingMotion(
            name=name,
            omega=omega,
            origin=origin,
            axis=axis
        ))

    def add_oscillating_rotating_motion(
            self,
            name: str,
            omega: float,
            amplitude: Vector,
            origin: Vector=Vector(0, 0, 0),
        ):
        """
        Oscillating Rotating Motion used in dynamicMeshDict.

        Parameters
        ----------
        name : str
            Name of the motion.
        omega : float
            Angular frequency in rad/s.
        amplitude : Vector
            Amplitude in degrees.
        origin : Vector
            Origin position of the rotation in m (default `Vector(0, 0, 0)`).
        """
        self.add_motion(OscillatingRotatingMotion(
            name=name,
            omega=omega,
            amplitude=amplitude,
            origin=origin
        ))

    def add_linear_motion(
            self,
            name: str,
            velocity: Vector
        ):
        """
        Add a linear Motion used in dynamicMeshDict.

        Parameters
        ----------
        name : str
            Name of the motion.
        velocity : Vector
            Motion velocity in m/s.
        """
        self.add_motion(LinearMotion(
            name=name,
            velocity=velocity
        ))

    def add_oscillating_linear_motion(
            self,
            name: str,
            amplitude: Vector,
            omega: float,
            timeShift: float = None,
            amplitudeShift: Vector = None
        ) -> None:
        """
        Add an Oscillating Linear Motion used in dynamicMeshDict. The motion is
        defined by the following formula:

            displacement(t) = amplitude * sin(omega * (t + timeShift)) + amplitudeShift

        Parameters
        ----------
        name : str
            Name of the motion.
        amplitude : Vector
            Amplitude vector of the oscilation in m.
        omega : float
            Angular frequency in rad/s.
        timeShift : float
            Phase shift in s (default `None`).
        amplitudeShift : Vector
            Vertical shift in m (default `None`).
        """
        self.add_motion(OscillatingLinearMotion(
            name=name,
            amplitude=amplitude,
            omega=omega,
            timeShift=timeShift,
            amplitudeShift=amplitudeShift
        ))


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""
        text += addParameter("dynamicFvMesh", self.dynamicFvMesh, isAddExtraLine=True)
        text += addParameter("motionSolverLibs", self.motionSolverLibs, isAddExtraLine=True)

        for key, value in self.extraParameters.items():
            text += addParameter(key, value, isAddExtraLine=True)

        if (self.dynamicFvMesh == "dynamicMotionSolverFvMesh"):
            check_type("motionSolver", self.motionSolver, str)
            check_value("motionSolver", self.motionSolver, _MOTION_SOLVER_TYPES)

            text += addParameter("motionSolver", self.motionSolver, isAddExtraLine=True)

            if (self.motionSolver == "solidBody"):
                check_type("solidBodyMotionFunction", self.solidBodyMotionFunction, str)
                check_value("solidBodyMotionFunction", self.solidBodyMotionFunction, _SOLID_BODY_MOTION_FUNCTION_TYPES)

                text += addParameter("solidBodyMotionFunction", self.solidBodyMotionFunction, isAddExtraLine=True)

                for motion in self.motions:
                    text += f"{motion!r}\n"

            elif (self.motionSolver == "multiSolidBodyMotionSolver"):
                text += "multiSolidBodyMotionSolverCoeffs\n{\n"

                for motion in self.motions:
                    if (isinstance(motion, Tabulated6DoFMotion)):
                        text += f"{tab}{motion.__repr__(depth=1, region=self.region)}\n"
                    else:
                        text += f"{tab}{motion.__repr__(depth=1)}\n"

                text += "}\n"

            elif (
                self.motionSolver == "displacementSBRStress" or
                self.motionSolver == "displacementLaplacian" or
                self.motionSolver == "velocityLaplacian"
            ):
                text += addParameter("diffusivity", self.diffusivity, isAddExtraLine=True)


        else:
            msg = f"dynamicFvMesh named '{self.dynamicFvMesh}' is not yet implemented"
            raise NotImplementedError(msg)

        return(text)
