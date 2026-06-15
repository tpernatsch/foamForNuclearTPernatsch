import numpy as np
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import OpenFOAMDict, Vector


_PATCH_TYPES = {
    # Basic
    "calculated", "fixedValue", "fixedGradient", "zeroGradient", "mixed",
    "directionMixed", "extrapolatedCalculated",
    # Constraint
    "cyclic", "cyclicACMI", "cyclicAMI", "cyclicSlip", "empty", "jumpCyclic",
    "jumpCyclicAMI", "nonuniformTransformCyclic", "processor", "processorCyclic",
    "symmetry", "symmetryPlane", "wedge",
    # Inlet
    "cylindricalInletVelocity", "fanPressure", "fixedFluxExtrapolatedPressure",
    "fixedFluxPressure", "fixedMean", "fixedMeanOutletInlet",
    "fixedNormalInletOutletVelocity", "fixedPressureCompressibleDensity",
    "flowRateInletVelocity", "freestream", "freestreamPressure",
    "freestreamVelocity", "mappedFlowRate", "mappedVelocityFluxFixedValue",
    "outletInlet", "outletMappedUniformInlet", "plenumPressure",
    "pressureDirectedInletOutletVelocity", "pressureDirectedInletVelocity",
    "pressureInletOutletParSlipVelocity", "pressureInletOutletVelocity",
    "pressureInletUniformVelocity", "pressureInletVelocity",
    "pressureNormalInletOutletVelocity", "pressurePIDControlInletVelocity",
    "rotatingPressureInletOutletVelocity", "rotatingTotalPressure",
    "supersonicFreestream", "surfaceNormalFixedValue",
    "swirlFlowRateInletVelocity", "swirlInletVelocity", "syringePressure",
    "timeVaryingMappedFixedValue", "totalPressure", "totalTemperature",
    "turbulentDFSEMInlet", "turbulentDigitalFilterInlet", "turbulentInlet",
    "turbulentIntensityKineticEnergyInlet", "uniformNormalFixedValue",
    "uniformTotalPressure", "variableHeightFlowRateInletVelocity",
    "variableHeightFlowRate", "waveSurfacePressure",
    # Outlet
    "advective", "fixedNormalInletOutletVelocity", "flowRateOutletVelocity",
    "fluxCorrectedVelocity", "inletOutlet", "inletOutletTotalTemperature",
    "matchedFlowRateOutletVelocity", "outletPhaseMeanVelocity",
    "uniformInletOutlet", "waveTransmissive",
    # Wall
    "fixedNormalSlip", "movingWallVelocity", "noSlip", "partialSlip",
    "rotatingWallVelocity", "slip", "translatingWallVelocity",
    "externalWallHeatFluxTemperature", "lumpedMassWallTemperature",
    "wallHeatTransfer",
    # Coupled
    "activeBaffleVelocity", "activePressureForceBaffleVelocity", "fan",
    "fixedJumpAMI", "fixedJump", "mappedField", "mappedFixedInternalValue",
    "mappedFixedPushedInternalValue", "mappedFixedValue", "mappedFlowRate",
    "mappedVelocityFluxFixedValue", "swirlFanVelocity", "temperatureCoupledBase",
    "timeVaryingMappedFixedValue", "uniformJumpAMI", "uniformJump",
    # Generic
    "mapped", "mappedVelocityFlux",
    "codedFixedValue", "codedMixed", "fixedInternalValueFvPatchField",
    "fixedProfile", "interfaceCompression", "phaseHydrostaticPressure",
    "prghPressure", "prghTotalHydrostaticPressure", "prghTotalPressure",
    "scaledFixedValue", "uniformDensityHydrostaticPressure",
    "uniformFixedGradient", "uniformFixedValue",
    # GeN-Foam
    "albedoSP3", "NusseltThermalBaffle1D",
    "epsilonWallFunction", "kqRWallFunction", "nutkWallFunction",
    "fixedMassFlowRate",
    # OFFBEAT
    "resistiveGap", "fuelRodGap", "gapContact", "tractionDisplacement", "coolantPressure",
    "fixedDisplacementZeroShear", "plenumSpringPressure", "topCladRingPressure",
    "implicitGapContact", "gapPressure",
    "unilateralContact", 
    "coolantChannel", "coolantChannelRIA", 
    # FMU4FOAM
    "coupledUniformExternalValue", "coupledFlowRateInletVelocity",
    "coupledFlowRateOutletVelocity", "coupledUniformExternalValue",
    # Custom
    "customPatch"
}


class Patch(OpenFOAMDict):
    """
    Patch base class.

    Parameters
    ----------
    type : str
        Patch type
    value : float | int | Vector
        Value at the patch
    """
    def __init__(
            self,
            type: str="empty",
            value: float | int | Vector=None
        ):
        super().__init__()
        self.type = type
        self.value = value


    def __repr__(self, depth = 0):
        if (not self.is_type_no_value() and self.value is not None):
            if (type(self.value) in [float, int, Vector]):
                self.__setitem__('value', f"uniform {self.value}")

        return super().__repr__(depth)


    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _PATCH_TYPES)
        self._type = type
        if (type != "customPatch"):
            self.__setitem__('type', self.type)

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value) -> None:
        if (value is not None):
            check_type("value", value, (int, float, list, Vector, np.ndarray))
            if (isinstance(value, (list, np.ndarray)) and len(value) == 3):
                self._value = Vector(value[0], value[1], value[2])
            else:
                self._value = value

            if (not self.is_type_no_value()):
                if (type(self.value) in [float, int, Vector]):
                    self.__setitem__('value', f"uniform {self.value}")
        else:
            self._value = None

    def is_type_no_value(self):
        return(self.type in [
            'zeroGradient', 'empty', 'slip', 'noSlip', 'wedge',
            'uniformFixedValue', 'cyclic', 'cyclicAMI'
        ])
