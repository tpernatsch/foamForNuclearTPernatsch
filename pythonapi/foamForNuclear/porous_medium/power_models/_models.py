import numpy as np
import matplotlib.pyplot as plt
from foamForNuclear.checkvalue import check_positive, check_type, check_value
from foamForNuclear.common import List, OpenFOAMDict, Polynome, Table, OffbeatDict
from foamForNuclear.timeProfile import TimeProfile

from foamForNuclear._attrs_tools import offbeat_define, _to_List_str
from typing import ClassVar
from attrs import field, validators as v


@offbeat_define
class PowerModel(OffbeatDict):
    TYPE: ClassVar[str] = "none"


@offbeat_define
class FixedPower(PowerModel):
    """
    Describes a generic structure with a constant or time depdendent internal
    power production.
    Solves an energy equation for the structure, updating the average
    structure temperature (assumed equal to the structure surface temperature
    for fluid heat transfer purposes).
    """
    TYPE: ClassVar[str] = "fixedPower"
    volumetricArea: int | float
    T:  int | float
    Cp:  int | float
    rho:  int | float
    powerDensity:  int | float
    powerTimeProfile: TimeProfile | None = None


@offbeat_define
class FixedTemperature(PowerModel):
    """
    Describes a structure with a constant surface temperature
    """
    TYPE: ClassVar[str] = "fixedTemperature"
    volumetricArea: int | float
    T:  int | float
    Cp:  int | float
    rho:  int | float
    powerDensity:  int | float
    temperatureTimeProfile: TimeProfile | None = None


@offbeat_define
class FixedTemperatureFMU(PowerModel):
    """
    Describes a structure with a constant surface temperature controlled through
    an FMI port.

    Parameters
    ----------
    temperatureNameFromFMU : str
        Name of the FMI port
    volumetricArea : float
        Volumetric area in m2/m3
    T : float
        Initial temperature
    """
    TYPE: ClassVar[str] = "fixedTemperatureFMU"
    temperatureNameFromFMU: str
    volumetricArea: int | float
    T:  int | float


@offbeat_define
class HeatedPin(PowerModel):
    """
    Model for representing a heated pin with constant material properties that
    is coupled to the fluid(s) via a convective boundary condition. The
    equation is solved via the finite volume method.

    Parameters
    ----------
    innerRadius : float
        Pin inner radius.
    outerRadius : float
        Pin outer radius.
    meshSize : int
        Pin radial mesh.
    k : float
        Pin thermal conductivity.
    T : float
        Pin initial temperature.
    rho : float
        Pin density in kg/m3.
    Cp : float
        Pin specific heat capacity in J/kg/K
    rhoCp : float
        If `rho` and `Cp` are not provided `rhoCp` can be used.
    """
    TYPE: ClassVar[str] = "heatedPin"
    innerRadius: int | float
    outerRadius: int | float
    meshSize: int | float
    k:  int | float
    T:  int | float
    rho:  int | float | None = None
    Cp:  int | float | None = None
    rhoCp:  int | float | None = None


@offbeat_define
class NuclearFuelPin(PowerModel):
    """
    Model for representing a nuclear fuel pin with constant material
    properties. Models fuel, gap, cladding and solves for the heat equation
    with a finite differences scheme. Only radial heat diffusion is considered.
    Derivation of the equations is provided here:
    https://gitlab.com/foam-for-nuclear/GeN-Foam/-/blob/develop/src/classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/nuclearFuelPin/nuclearFuelPin.H?ref_type=heads

    Parameters
    ----------
    powerDensity : float
        Initial uniform power density if not computed by the neutronics
        sub-solver.
    fuelInnerRadius : float
        Fuel pellet inner radius.
    fuelOuterRadius : float
        Fuel pellet outer radius.
    cladInnerRadius : float
        Cladding inner radius.
    cladOuterRadius : float
        Cladding outer radius.
    fuelMeshSize : int
        Fuel pellet radial mesh size.
    cladMeshSize : int
        Cladding radial mesh size.
    fuelRho : float
        Fuel constant density in kg/m3.
    fuelCp : float
        Fuel constant specific heat capacity in J/kg/K.
    fuelK : float
        Fuel constant thermal conductivity in W/m/K.
    cladRho : float
        Cladding constant density in kg/m3.
    cladCp : float
        Cladding constant specific heat capacity in J/kg/K.
    cladK : float
        Cladding constant thermal conductivity in W/m/K.
    fuelT : float
        Fuel initial temperature.
    cladT : float
        Cladding initial temperature.
    gapH : float
        Gap gas constant heat transfer coefficient.
    gapHPowerDensity : Table
        Gap gas power-dependent heat transfer coefficient. Can be used in
        replacement of `gapH`.
    """
    TYPE: ClassVar[str] = "nuclearFuelPin"
    fuelInnerRadius: int | float
    fuelOuterRadius: int | float
    cladInnerRadius: int | float
    cladOuterRadius: int | float
    fuelMeshSize: int | float
    cladMeshSize: int | float
    fuelK: int | float
    fuelRho: int | float
    fuelCp: int | float
    cladK: int | float
    cladRho: int | float
    cladCp: int | float
    fuelT: int | float
    cladT: int | float            
    powerDensity: int | float = 0
    gapH: int | float | None = None
    gapHPowerDensity: int | float | None = None


@offbeat_define
class NuclearSteadyStatePebble(PowerModel):
    """
    Describes triso and pebble of a pebble-bed reactor. Only steady-state.
    """
    TYPE: ClassVar[str] = "nuclearSteadyStatePebble"
    powerDensityNeutronics:  int | float
    pebbleSurfaceTemperatureGuess:  int | float
    pebbleCoreRadius:  int | float
    pebbleMatrixRadius:  int | float
    pebbleShellRadius:  int | float
    trisoFuelRadius:  int | float
    trisoBufferRadius:  int | float
    trisoInnerPyCRadius:  int | float
    trisoSiCRadius:  int | float
    trisoOuterPyCRadius:  int | float
    nTRISO: int
    trisoFuelKCoeffs: list[int | float] = field(factory=list)
    trisoBufferKCoeffs: list[int | float] = field(factory=list)
    trisoPyCKCoeffs: list[int | float] = field(factory=list)
    trisoSiCKCoeffs: list[int | float] = field(factory=list)
    pebbleGraphiteKCoeffs: list[int | float] = field(factory=list)
    trisoFuelDensityCoeffs: list[int | float] = field(factory=list)
    trisoBufferDensityCoeffs: list[int | float] = field(factory=list)
    trisoPyCDensityCoeffs: list[int | float] = field(factory=list)
    trisoSiCDensityCoeffs: list[int | float] = field(factory=list)
    pebbleGraphiteDensityCoeffs: list[int | float] = field(factory=list)
    coolantDensityCoeffs: list[int | float] = field(factory=list)


# @offbeat_define
# class LumpedNuclearStructure(PowerModel):
#     r"""
#     Model for representing a lumped-parameter structure characterized
#     by a user selectable number of nodes (indicated below as T[idx]) ::

#         Tmax            T[0]           T[1]          T[n-1]         Tsurface
#           | --- H[0] --- | --- H[1] --- | --- H[2] --- |  --- H[n] --- |

#     Note that only T[idx] are state variables while Tmax and Tsurface can
#     be found after matrix solution.
#     Note also that one needs n+1 conductances when solving for n nodes
#     The equation in the class can be obtained starting from a simple
#     energy balance for the ith node:

#     .. math::
#         V_i \rho_i c_{p,i} \frac{dT}{dt} =
#            V Q q_f + H_i (T_{i-1} - T_i) - H_{i+1}(T_i - T_{i+1})

#     and dividing all terms by the volume of the structure.
#     This means that the heat conductances  should be calculated as
#     the power flowing between two nodes divided by the diffence in
#     temperatures at steadystate, abd divided by the volume of the
#     structure.
#     Note that Q is the power density from neutronincs, which is
#     assumed to spead over all structure. The power fraction qf
#     can be used to concentrate this power over certain nodes.
#     Zero-gradient BCs are used at the inner boundary and convective
#     BCs are used at the external boundary (in contact with the coolant)
#     """
#     TYPE: ClassVar[str] = "lumpedNuclearStructure"

#     def __init__(
#             self,
#             volumetricArea,
#             powerDensity,
#             nodeFuel,
#             nodeClad,
#             nodeMatrix,
#             kappaMatrix,
#             heatConductances,
#             rhoCp,
#             volumeFractions,
#             powerFractions,
#             T0,
#         ):
#         super().__init__("lumpedNuclearStructure")

#         self.volumetricArea = volumetricArea
#         self.powerDensity = powerDensity
#         self.nodesNumber = len(volumeFractions)
#         self.nodeFuel = nodeFuel
#         self.nodeClad = nodeClad
#         self.nodeMatrix = nodeMatrix
#         self.kappaMatrix = kappaMatrix
#         self.heatConductances = heatConductances
#         self.rhoCp = rhoCp
#         self.volumeFractions = volumeFractions
#         self.powerFractions = powerFractions
#         self.T0 = T0


#     def __repr__(self, depth=0):
#         self.nodesNumber = len(self.volumeFractions)

#         if (self.nodesNumber != len(self.rhoCp)):
#             msg = f"volumeFractions must have the same length as rhoCp, {self.nodesNumber} != {len(self.rhoCp)}"
#             raise ValueError(msg)
#         if (self.nodesNumber != len(self.powerFractions)):
#             msg = f"volumeFractions must have the same length as powerFractions, {self.nodesNumber} != {len(self.powerFractions)}"
#             raise ValueError(msg)
#         if (self.nodesNumber+1 != len(self.heatConductances)):
#             msg = f"volumeFractions must be one less than the length of heatConductances, {self.nodesNumber}+1 != {len(self.heatConductances)}"
#             raise ValueError(msg)

#         self.__setitem__("volumetricArea", self.volumetricArea)
#         self.__setitem__("powerDensity", self.powerDensity)
#         self.__setitem__("nodesNumber", self.nodesNumber)
#         self.__setitem__("nodeFuel", self.nodeFuel)
#         self.__setitem__("nodeClad", self.nodeClad)
#         self.__setitem__("nodeMatrix", self.nodeMatrix)
#         self.__setitem__("kappaMatrix", self.kappaMatrix)
#         self.__setitem__("heatConductances", self.heatConductances)
#         self.__setitem__("rhoCp", self.rhoCp)
#         self.__setitem__("volumeFractions", self.volumeFractions)
#         self.__setitem__("powerFractions", self.powerFractions)
#         self.__setitem__("T0", self.T0)

#         return super().__repr__(depth)


# TODO: Complete migration to attrs also, keeping the current checks on the attributes
# @offbeat_define
# class XYPosLattice(OpenFOAMDict):
#     """
#     Parameters
#     ----------
#     type : str
#         Lattice type. Options: `hexagonal`, `square`
#     nElements : List | list
#         Number of elements per rows and per column
#     pitch: float
#         Lattice element distance
#     lattice: List | list | str
#         Lattice element placement (`'e'` == no placement, else is placed). Format::

#             e e F F F       #  y -->
#              e F F F F      # x
#               F F F F F     # |
#                F F F F e    # v
#                 F F F e e

#     origin: List
#         Origin center of the lattice in absolute coordinates (default `[0, 0]`).
#     rotationAngle: float
#         Rotation angle of the lattice (default `0`).
#     """
#     type: str
#     nElements: List | list = field(factory=list)
#     pitch: float
#     lattice: List | list | str
#     origin: List | list = field(default=list[0, 0])
#     rotationAngle: int | float = 0

#     def __repr__(self, depth = 0):
#         self.__setitem__("lattice", self.lattice.__repr__(depth=depth+1, nCols=self.nElements[0]))
#         return super().__repr__(depth)

#     @property
#     def type(self):
#         return self._type_

#     @type.setter
#     def type(self, type) -> None:
#         check_type("type", type, str)
#         check_value("type", type, {"hexagonal", "square"})
#         self.__setitem__("type", type)
#         self._type_ = type

#     @property
#     def nElements(self):
#         return self._nElements

#     @nElements.setter
#     def nElements(self, nElements) -> None:
#         check_type("nElements", nElements, (List, list))
#         if (len(nElements) != 2):
#             msg = f"nElements must have a length of 2, found {len(nElements)} elements"
#             raise ValueError(msg)
#         if (isinstance(nElements, list)):
#             self._nElements = List(nElements)
#         else:
#             self._nElements = nElements
#         self.__setitem__("nElements", f"{self._nElements!r}")

#     @property
#     def pitch(self):
#         return self._pitch

#     @pitch.setter
#     def pitch(self, pitch) -> None:
#         check_type("pitch", pitch, (float, int))
#         check_positive("pitch", pitch, is_strict=True)
#         self._pitch = pitch
#         self.__setitem__("pitch", pitch)

#     @property
#     def origin(self):
#         return self._origin

#     @origin.setter
#     def origin(self, origin) -> None:
#         check_type("origin", origin, (List, list))
#         if (len(origin) != 2):
#             msg = f"origin must have a length of 2, found {len(origin)} elements"
#             raise ValueError(msg)
#         if (isinstance(origin, list)):
#             self._origin = List(origin)
#         else:
#             self._origin = origin
#         self.__setitem__("origin", f"{self._origin!r}")

#     @property
#     def lattice(self):
#         return self._lattice

#     @lattice.setter
#     def lattice(self, lattice) -> None:
#         check_type("lattice", lattice, (List, list, str))
#         if (isinstance(lattice, list)):
#             self._lattice = List(lattice)
#         elif (isinstance(lattice, str)):
#             self._lattice = List(' '.join(lattice.split('\n')).split())
#         else:
#             self._lattice = lattice

#     @property
#     def rotationAngle(self):
#         return self._rotationAngle

#     @rotationAngle.setter
#     def rotationAngle(self, rotationAngle) -> None:
#         check_type("rotationAngle", rotationAngle, (float, int))
#         self._rotationAngle = rotationAngle
#         self.__setitem__("rotationAngle", rotationAngle)


#     def computeCoordinates(self):
#         if (self.type == "hexagonal"):
#             return(self.computeCoordinatesHexagon())
#         elif (self.type == "square"):
#             return(self.computeCoordinatesSquare())

#     def computeCoordinatesSquare(self):
#         nx = self.nElements[0]
#         ny = self.nElements[1]
#         listX, listY = [], []
#         xmin = -self.pitch*nx/2
#         ymin = -self.pitch*ny/2
#         nullElement = "e"
#         for i in range(nx):
#             for j in range(ny):
#                 if (self.lattice[ny*i+j] != nullElement):
#                     x = xmin + i*self.pitch
#                     y = ymin + j*self.pitch

#                     listX.append(
#                         self.origin[0] + x*np.cos(self.rotationAngle) - y*np.sin(self.rotationAngle)
#                     )
#                     listY.append(
#                         self.origin[1] + x*np.sin(self.rotationAngle) + y*np.cos(self.rotationAngle)
#                     )

#         return(listX, listY)

#     def computeCoordinatesHexagon(self):
#         nx = self.nElements[0]
#         ny = self.nElements[1]
#         listX, listY = [], []
#         summitToSummit = self.pitch * 2.0/np.sqrt(3.0)
#         nullElement = "e"
#         for i in range(nx):
#             for j in range(ny):
#                 if (self.lattice[ny*i+j] != nullElement):
#                     x = 0.75 * summitToSummit * (i - (nx-1.0)/2.0)
#                     y = self.pitch * ((j - (ny-1.0)/2.0) + 0.5 * (i - (nx-1.0)/2.0))

#                     listX.append(
#                         self.origin[0] + x*np.cos(self.rotationAngle) - y*np.sin(self.rotationAngle)
#                     )
#                     listY.append(
#                         self.origin[1] + x*np.sin(self.rotationAngle) + y*np.cos(self.rotationAngle)
#                     )

#         return(listX, listY)

#     def plotLattice(self, ax=None):
#         """
#         Parameters
#         ----------
#         ax
#             Matplotlib axes object e.g `fig, ax = plt.subplots()`
#             (default `None`).
#         """
#         axNone = ax is None
#         if (axNone):
#             fig, ax = plt.subplots()
#         # circle1 = patches.Circle((0, 0), radius=1, color="lightgrey")
#         # ax.add_patch(circle1)
#         listX, listY = self.computeCoordinates()
#         nElements = len(listX)
#         ax.scatter(listX, listY, label=f"{nElements} elements")
#         ax.set_xlabel("X [m]")
#         ax.set_ylabel("Y [m]")
#         # ax.set_xticks([0], minor=False)
#         # ax.xaxis.grid(True, which='major')
#         # ax.set_yticks([0], minor=False)
#         # ax.yaxis.grid(True, which='major')
#         ax.set_aspect('equal', adjustable='box')
#         if (axNone):
#             ax.legend(
#                 bbox_to_anchor=(0., 1.02, 1., .102),
#                 loc='lower left',
#                 ncols=2,
#                 mode="expand",
#                 borderaxespad=0.
#             )
#             fig.savefig(f"fig_lattice_placement_{self.type}_{nElements}elements.png")
#             plt.close()


# class NuclearFuelFMU(PowerModel):
#     """
#     Parameters
#     ----------
#     axialPowerInterpolationMethod : str
#         Options: `linear`
#     radialBasisFunctionMethod : str
#         Options: `gaussian`, `polyharmonicSpline`, `kriging`
#     epsilon : List
#         Epsilon factors (one per direction) for Gaussion Radial Basis Function
#     """

#     def __init__(
#             self,
#             volumetricArea: float,
#             axialPowerInterpolationMethod: str,
#             axialLocationsNameToFMU: str,
#             axialLocations: List,
#             fuelLength: float,
#             radialBasisFunctionMethod: str,
#             avgPowerDensityNameToFMU: List,
#             axialProfilePowerDensityNameToFMU: List,
#             initPowerDensity: float,
#             initTstruct: float,
#             initHeatFlux: float=None,
#             initRhoCpdTdt: float=None,
#             heatFluxNameFromFMU: List=None,
#             rhoCpdTdtNameFromFMU: List=None,
#             htcNameToFMU: List=None,
#             TfluidNameToFMU: List=None,
#             TFuelNameFromFMU: List=None,
#             TCladNameFromFMU: List=None,
#             TstructNameToFMU: List=None,
#             fuelFraction: float=1,
#             fractionOfPowerFromNeutronics: float=1,
#             xyPosLattice: XYPosLattice=None,
#             xPos: List=None,
#             yPos: List=None,
#             krigingOptions: OpenFOAMDict=None,
#             krigingOptionsFuel: OpenFOAMDict=None,
#             krigingOptionsClad: OpenFOAMDict=None,
#             epsilon: List=None,
#             relaxationFactorTfluid: float=1,
#             relaxationFactorTstruct: float=1,
#             relaxationFactorPowerDensity: float=1,
#             relaxationFactorHtc: float=1,
#         ):
#         super().__init__("nuclearFuelFMU")

#         self.volumetricArea = volumetricArea
#         self.axialPowerInterpolationMethod = axialPowerInterpolationMethod
#         self.axialLocationsNameToFMU = axialLocationsNameToFMU
#         self.axialLocations = axialLocations
#         self.fuelLength = fuelLength
#         self.fuelFraction = fuelFraction
#         self.fractionOfPowerFromNeutronics = fractionOfPowerFromNeutronics
#         self.radialBasisFunctionMethod = radialBasisFunctionMethod
#         self.initPowerDensity = initPowerDensity
#         self.initTstruct = initTstruct
#         self.initHeatFlux = initHeatFlux
#         self.initRhoCpdTdt = initRhoCpdTdt
#         self.avgPowerDensityNameToFMU = avgPowerDensityNameToFMU
#         self.axialProfilePowerDensityNameToFMU = axialProfilePowerDensityNameToFMU
#         self.htcNameToFMU = htcNameToFMU
#         self.heatFluxNameFromFMU = heatFluxNameFromFMU
#         self.rhoCpdTdtNameFromFMU = rhoCpdTdtNameFromFMU
#         self.TfluidNameToFMU = TfluidNameToFMU
#         self.TFuelNameFromFMU = TFuelNameFromFMU
#         self.TCladNameFromFMU = TCladNameFromFMU
#         self.TstructNameToFMU = TstructNameToFMU
#         self.xyPosLattice = xyPosLattice
#         self.xPos = xPos
#         self.yPos = yPos
#         self.krigingOptions = krigingOptions
#         self.krigingOptionsFuel = krigingOptionsFuel
#         self.krigingOptionsClad = krigingOptionsClad
#         self.epsilon = epsilon
#         self.relaxationFactorTfluid = relaxationFactorTfluid
#         self.relaxationFactorTstruct = relaxationFactorTstruct
#         self.relaxationFactorPowerDensity = relaxationFactorPowerDensity
#         self.relaxationFactorHtc = relaxationFactorHtc


#     def __repr__(self, depth=0):
#         self.__setitem__("axialLocations", self.axialLocations.__repr__(depth=depth+1))

#         if (self.xPos is None and self.yPos is None and self.xyPosLattice is None):
#             msg = "None of (xPos, yPos) or xyPosLattice have been provided. It is mandatory to add XY-locations samples."
#             raise ValueError(msg)
#         if (self.xyPosLattice is not None):
#             self.__setitem__("xyPosLattice", self.xyPosLattice)
#         if (self.xPos is not None and self.yPos is not None):
#             self.__setitem__("xPos", self.xPos)
#             self.__setitem__("yPos", self.yPos)

#         self.__setitem__("avgPowerDensityNameToFMU", self.avgPowerDensityNameToFMU.__repr__(depth=depth+1))
#         self.__setitem__("axialProfilePowerDensityNameToFMU", self.axialProfilePowerDensityNameToFMU.__repr__(depth=depth+1))

#         if (self.TfluidNameToFMU is not None):
#             self.__setitem__("TfluidNameToFMU", self.TfluidNameToFMU.__repr__(depth=depth+1))
#         if (self.TFuelNameFromFMU is not None):
#             self.__setitem__("TFuelNameFromFMU", self.TFuelNameFromFMU.__repr__(depth=depth+1))
#         if (self.TCladNameFromFMU is not None):
#             self.__setitem__("TCladNameFromFMU", self.TCladNameFromFMU.__repr__(depth=depth+1))
#         if (self.TstructNameToFMU is not None):
#             self.__setitem__("TstructNameToFMU", self.TstructNameToFMU.__repr__(depth=depth+1))
#         if (self.htcNameToFMU is not None):
#             self.__setitem__("htcNameToFMU", self.htcNameToFMU.__repr__(depth=depth+1))
#         if (self.heatFluxNameFromFMU is not None):
#             self.__setitem__("heatFluxNameFromFMU", self.heatFluxNameFromFMU.__repr__(depth=depth+1))
#             self.__setitem__("initHeatFlux", self.initHeatFlux)
#         if (self.rhoCpdTdtNameFromFMU is not None):
#             self.__setitem__("rhoCpdTdtNameFromFMU", self.rhoCpdTdtNameFromFMU.__repr__(depth=depth+1))
#             self.__setitem__("initRhoCpdTdt", self.initRhoCpdTdt)

#         if (self.radialBasisFunctionMethod == 'kriging'):
#             if (self.krigingOptions is not None and self.krigingOptionsFuel is not None and self.krigingOptionsClad is not None):
#                 self.__setitem__("krigingOptions", self.krigingOptions)
#                 self.__setitem__("krigingOptionsFuel", self.krigingOptionsFuel)
#                 self.__setitem__("krigingOptionsClad", self.krigingOptionsClad)
#             else:
#                 msg = "Kriging needs 'krigingOptions', 'krigingOptionsFuel' and 'krigingOptionsClad' parameters"
#                 raise ValueError(msg)

#         if (self.radialBasisFunctionMethod == 'gaussian'):
#             if (self.epsilon is not None):
#                 self.__setitem__("epsilon", self.epsilon)
#             else:
#                 msg = "Gaussian RBF needs 'epsilon' parameters"
#                 raise ValueError(msg)

#         return super().__repr__(depth)

#     @property
#     def volumetricArea(self):
#         return self._volumetricArea

#     @volumetricArea.setter
#     def volumetricArea(self, volumetricArea) -> None:
#         check_type("volumetricArea", volumetricArea, (float, int))
#         check_positive("volumetricArea", volumetricArea)
#         self.__setitem__("volumetricArea", volumetricArea)
#         self._volumetricArea = volumetricArea

#     @property
#     def axialPowerInterpolationMethod(self):
#         return self._axialPowerInterpolationMethod

#     @axialPowerInterpolationMethod.setter
#     def axialPowerInterpolationMethod(self, axialPowerInterpolationMethod) -> None:
#         check_type("axialPowerInterpolationMethod", axialPowerInterpolationMethod, str)
#         check_value("axialPowerInterpolationMethod", axialPowerInterpolationMethod, {"linear"})
#         self.__setitem__("axialPowerInterpolationMethod", axialPowerInterpolationMethod)
#         self._axialPowerInterpolationMethod = axialPowerInterpolationMethod

#     @property
#     def axialLocationsNameToFMU(self):
#         return self._axialLocationsNameToFMU

#     @axialLocationsNameToFMU.setter
#     def axialLocationsNameToFMU(self, axialLocationsNameToFMU) -> None:
#         check_type("axialLocationsNameToFMU", axialLocationsNameToFMU, str)
#         self.__setitem__("axialLocationsNameToFMU", axialLocationsNameToFMU)
#         self._axialLocationsNameToFMU = axialLocationsNameToFMU

#     @property
#     def axialLocations(self):
#         return self._axialLocations

#     @axialLocations.setter
#     def axialLocations(self, axialLocations) -> None:
#         check_type("axialLocations", axialLocations, (List, list))
#         if (isinstance(axialLocations, list)):
#             self._axialLocations = List(axialLocations)
#         else:
#             self._axialLocations = axialLocations

#     @property
#     def fuelLength(self):
#         return self._fuelLength

#     @fuelLength.setter
#     def fuelLength(self, fuelLength) -> None:
#         check_type("fuelLength", fuelLength, (float, int))
#         check_positive("fuelLength", fuelLength)
#         self.__setitem__("fuelLength", fuelLength)
#         self._fuelLength = fuelLength

#     @property
#     def fuelFraction(self):
#         return self._fuelFraction

#     @fuelFraction.setter
#     def fuelFraction(self, fuelFraction) -> None:
#         check_type("fuelFraction", fuelFraction, (float, int))
#         check_positive("fuelFraction", fuelFraction)
#         self.__setitem__("fuelFraction", fuelFraction)
#         self._fuelFraction = fuelFraction

#     @property
#     def fractionOfPowerFromNeutronics(self):
#         return self._fractionOfPowerFromNeutronics

#     @fractionOfPowerFromNeutronics.setter
#     def fractionOfPowerFromNeutronics(self, fractionOfPowerFromNeutronics) -> None:
#         check_type("fractionOfPowerFromNeutronics", fractionOfPowerFromNeutronics, (float, int))
#         check_positive("fractionOfPowerFromNeutronics", fractionOfPowerFromNeutronics)
#         self.__setitem__("fractionOfPowerFromNeutronics", fractionOfPowerFromNeutronics)
#         self._fractionOfPowerFromNeutronics = fractionOfPowerFromNeutronics

#     @property
#     def initRhoCpdTdt(self):
#         return self._initRhoCpdTdt

#     @initRhoCpdTdt.setter
#     def initRhoCpdTdt(self, initRhoCpdTdt) -> None:
#         check_type("initRhoCpdTdt", initRhoCpdTdt, (float, int), none_ok=True)
#         self._initRhoCpdTdt = initRhoCpdTdt

#     @property
#     def initHeatFlux(self):
#         return self._initHeatFlux

#     @initHeatFlux.setter
#     def initHeatFlux(self, initHeatFlux) -> None:
#         check_type("initHeatFlux", initHeatFlux, (float, int), none_ok=True)
#         self._initHeatFlux = initHeatFlux

#     @property
#     def initTstruct(self):
#         return self._initTstruct

#     @initTstruct.setter
#     def initTstruct(self, initTstruct) -> None:
#         check_type("initTstruct", initTstruct, (float, int))
#         self.__setitem__("initTstruct", initTstruct)
#         self._initTstruct = initTstruct

#     @property
#     def initPowerDensity(self):
#         return self._initPowerDensity

#     @initPowerDensity.setter
#     def initPowerDensity(self, initPowerDensity) -> None:
#         check_type("initPowerDensity", initPowerDensity, (float, int))
#         self.__setitem__("initPowerDensity", initPowerDensity)
#         self._initPowerDensity = initPowerDensity

#     @property
#     def radialBasisFunctionMethod(self):
#         return self._radialBasisFunctionMethod

#     @radialBasisFunctionMethod.setter
#     def radialBasisFunctionMethod(self, radialBasisFunctionMethod) -> None:
#         check_type("radialBasisFunctionMethod", radialBasisFunctionMethod, str)
#         check_value("radialBasisFunctionMethod", radialBasisFunctionMethod, {"gaussian", "polyharmonicSpline", "kriging"})
#         self.__setitem__("radialBasisFunctionMethod", radialBasisFunctionMethod)
#         self._radialBasisFunctionMethod = radialBasisFunctionMethod

#     @property
#     def avgPowerDensityNameToFMU(self):
#         return self._avgPowerDensityNameToFMU

#     @avgPowerDensityNameToFMU.setter
#     def avgPowerDensityNameToFMU(self, avgPowerDensityNameToFMU) -> None:
#         check_type("avgPowerDensityNameToFMU", avgPowerDensityNameToFMU, (List, list))
#         if (isinstance(avgPowerDensityNameToFMU, list)):
#             self._avgPowerDensityNameToFMU = List(avgPowerDensityNameToFMU)
#         else:
#             self._avgPowerDensityNameToFMU = avgPowerDensityNameToFMU

#     @property
#     def axialProfilePowerDensityNameToFMU(self):
#         return self._axialProfilePowerDensityNameToFMU

#     @axialProfilePowerDensityNameToFMU.setter
#     def axialProfilePowerDensityNameToFMU(self, axialProfilePowerDensityNameToFMU) -> None:
#         check_type("axialProfilePowerDensityNameToFMU", axialProfilePowerDensityNameToFMU, (List, list))
#         if (isinstance(axialProfilePowerDensityNameToFMU, list)):
#             self._axialProfilePowerDensityNameToFMU = List(axialProfilePowerDensityNameToFMU)
#         else:
#             self._axialProfilePowerDensityNameToFMU = axialProfilePowerDensityNameToFMU

#     @property
#     def TstructNameToFMU(self):
#         return self._TstructNameToFMU

#     @TstructNameToFMU.setter
#     def TstructNameToFMU(self, TstructNameToFMU) -> None:
#         check_type("TstructNameToFMU", TstructNameToFMU, (List, list), none_ok=True)
#         if (isinstance(TstructNameToFMU, list)):
#             self._TstructNameToFMU = List(TstructNameToFMU)
#         else:
#             self._TstructNameToFMU = TstructNameToFMU

#     @property
#     def TfluidNameToFMU(self):
#         return self._TfluidNameToFMU

#     @TfluidNameToFMU.setter
#     def TfluidNameToFMU(self, TfluidNameToFMU) -> None:
#         check_type("TfluidNameToFMU", TfluidNameToFMU, (List, list), none_ok=True)
#         if (isinstance(TfluidNameToFMU, list)):
#             self._TfluidNameToFMU = List(TfluidNameToFMU)
#         else:
#             self._TfluidNameToFMU = TfluidNameToFMU

#     @property
#     def TFuelNameFromFMU(self):
#         return self._TFuelNameFromFMU

#     @TFuelNameFromFMU.setter
#     def TFuelNameFromFMU(self, TFuelNameFromFMU) -> None:
#         check_type("TFuelNameFromFMU", TFuelNameFromFMU, (List, list), none_ok=True)
#         if (isinstance(TFuelNameFromFMU, list)):
#             self._TFuelNameFromFMU = List(TFuelNameFromFMU)
#         else:
#             self._TFuelNameFromFMU = TFuelNameFromFMU

#     @property
#     def TCladNameFromFMU(self):
#         return self._TCladNameFromFMU

#     @TCladNameFromFMU.setter
#     def TCladNameFromFMU(self, TCladNameFromFMU) -> None:
#         check_type("TCladNameFromFMU", TCladNameFromFMU, (List, list), none_ok=True)
#         if (isinstance(TCladNameFromFMU, list)):
#             self._TCladNameFromFMU = List(TCladNameFromFMU)
#         else:
#             self._TCladNameFromFMU = TCladNameFromFMU

#     @property
#     def htcNameToFMU(self):
#         return self._htcNameToFMU

#     @htcNameToFMU.setter
#     def htcNameToFMU(self, htcNameToFMU) -> None:
#         check_type("htcNameToFMU", htcNameToFMU, (List, list), none_ok=True)
#         if (isinstance(htcNameToFMU, list)):
#             self._htcNameToFMU = List(htcNameToFMU)
#         else:
#             self._htcNameToFMU = htcNameToFMU

#     @property
#     def heatFluxNameFromFMU(self):
#         return self._heatFluxNameFromFMU

#     @heatFluxNameFromFMU.setter
#     def heatFluxNameFromFMU(self, heatFluxNameFromFMU) -> None:
#         check_type("heatFluxNameFromFMU", heatFluxNameFromFMU, (List, list), none_ok=True)
#         if (isinstance(heatFluxNameFromFMU, list)):
#             self._heatFluxNameFromFMU = List(heatFluxNameFromFMU)
#         else:
#             self._heatFluxNameFromFMU = heatFluxNameFromFMU

#     @property
#     def rhoCpdTdtNameFromFMU(self):
#         return self._rhoCpdTdtNameFromFMU

#     @rhoCpdTdtNameFromFMU.setter
#     def rhoCpdTdtNameFromFMU(self, rhoCpdTdtNameFromFMU) -> None:
#         check_type("rhoCpdTdtNameFromFMU", rhoCpdTdtNameFromFMU, (List, list), none_ok=True)
#         if (isinstance(rhoCpdTdtNameFromFMU, list)):
#             self._rhoCpdTdtNameFromFMU = List(rhoCpdTdtNameFromFMU)
#         else:
#             self._rhoCpdTdtNameFromFMU = rhoCpdTdtNameFromFMU
