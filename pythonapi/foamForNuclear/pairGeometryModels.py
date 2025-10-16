from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import OpenFOAMDict, tab


class DispersionModel(OpenFOAMDict):
    def __init__(self, type: str):
        super().__init__(name="dispersionModel")

        self.type = type

    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        # check_value("type", type, _DISPERSED_DIAMETER_MODEL_TYPES)
        self._type_ = type
        self.__setitem__('type', type)


class ByRegimeDispersionModel(DispersionModel):
    def __init__(
            self,
            regimeMap: str,
            regimes: dict[OpenFOAMDict]={}
        ):
        super().__init__(type="byRegime")

        self.regimeMap = regimeMap
        self.regimes: dict[OpenFOAMDict] = regimes

    def __repr__(self, depth: int=0):
        for regimeName, regime in self.regimes.items():
            self.__setitem__(regimeName, regime)

        return super().__repr__(depth)

    @property
    def regimeMap(self):
        return self._regimeMap

    @regimeMap.setter
    def regimeMap(self, regimeMap) -> None:
        check_type("regimeMap", regimeMap, str)
        self._regimeMap = regimeMap
        self.__setitem__('regimeMap', regimeMap)

    def add_regime(
            self,
            regimeName: str,
            regimeDict: OpenFOAMDict
        ) -> None:
        check_type("regimeName", regimeName, str)
        check_type("regimeDict", regimeDict, OpenFOAMDict)
        self.regimes[regimeName] = regimeDict

    def add_constant_regime(
            self,
            regimeName: str,
            dispersedPhaseName: str
        ) -> None:
        check_type("regimeName", regimeName, str)
        check_type("dispersedPhaseName", dispersedPhaseName, str)
        self.add_regime(regimeName, OpenFOAMDict({
            "type": "constant",
            "dispersedPhase": dispersedPhaseName
        }))


class ConstantDispersionModel(DispersionModel):
    def __init__(
            self,
            dispersedPhase: str
        ):
        super().__init__(type="constant")

        self.dispersedPhase = dispersedPhase

    @property
    def dispersedPhase(self):
        return self._dispersedPhase

    @dispersedPhase.setter
    def dispersedPhase(self, dispersedPhase) -> None:
        check_type("dispersedPhase", dispersedPhase, str)
        self._dispersedPhase = dispersedPhase
        self.__setitem__('dispersedPhase', dispersedPhase)


class InterfacialAreaDensityModel(OpenFOAMDict):
    def __init__(self, type: str):
        super().__init__(name="interfacialAreaDensityModel")

        self.type = type

    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        # check_value("type", type, _DISPERSED_DIAMETER_MODEL_TYPES)
        self._type_ = type
        self.__setitem__('type', type)


class SphericalInterfacialAreaDensityModel(InterfacialAreaDensityModel):
    def __init__(self):
        super().__init__("spherical")


class ContactPartitionModel(OpenFOAMDict):
    def __init__(self, type: str):
        super().__init__(name="contactPartitionModel")

        self.type = type

    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        # check_value("type", type, _DISPERSED_DIAMETER_MODEL_TYPES)
        self._type_ = type
        self.__setitem__('type', type)


class ByRegimeContactPartitionModel(ContactPartitionModel):
    def __init__(
            self,
            regimeMap: str,
            regimes: dict[OpenFOAMDict]={}
        ):
        super().__init__(type="byRegime")

        self.regimeMap = regimeMap
        self.regimes: dict[OpenFOAMDict] = regimes

    def __repr__(self, depth: int=0):
        for regimeName, regime in self.regimes.items():
            self.__setitem__(regimeName, regime)

        return super().__repr__(depth)

    @property
    def regimeMap(self):
        return self._regimeMap

    @regimeMap.setter
    def regimeMap(self, regimeMap) -> None:
        check_type("regimeMap", regimeMap, str)
        self._regimeMap = regimeMap
        self.__setitem__('regimeMap', regimeMap)

    def add_regime(
            self,
            regimeName: str,
            regimeDict: OpenFOAMDict
        ) -> None:
        check_type("regimeName", regimeName, str)
        check_type("regimeDict", regimeDict, OpenFOAMDict)
        self.regimes[regimeName] = regimeDict

    def add_constant_regime(
            self,
            regimeName: str,
            value: float
        ) -> None:
        check_type("regimeName", regimeName, str)
        check_type("value", value, (float, int))
        self.add_regime(regimeName, OpenFOAMDict({
            "type": "constant",
            "value": value
        }))


class ConstantContactPartitionModel(ContactPartitionModel):
    def __init__(
            self,
            value: float
        ):
        super().__init__(type="constant")

        self.value = value

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value) -> None:
        check_type("value", value, (float, int))
        self._value = value
        self.__setitem__('value', value)


class PairGeometryModel:
    def __init__(
            self,
            dispersionModel: DispersionModel=None,
            interfacialAreaDensityModel: InterfacialAreaDensityModel=None,
            contactPartitionModel: ContactPartitionModel=None,
        ):
        self.dispersionModel = dispersionModel
        self.interfacialAreaDensityModel = interfacialAreaDensityModel
        self.contactPartitionModel = contactPartitionModel

    def __repr__(self, depth: int=0):
        text = ""
        if (self.dispersionModel is not None):
            text += self.dispersionModel.__repr__(depth=depth)
            text += depth*tab
        if (self.interfacialAreaDensityModel is not None):
            text += self.interfacialAreaDensityModel.__repr__(depth=depth)
            text += depth*tab
        if (self.contactPartitionModel is not None):
            text += self.contactPartitionModel.__repr__(depth=depth)
            text += depth*tab

        return(text)

    @property
    def dispersionModel(self):
        return self._dispersionModel

    @dispersionModel.setter
    def dispersionModel(self, dispersionModel) -> None:
        check_type("dispersionModel", dispersionModel, DispersionModel, none_ok=True)
        self._dispersionModel = dispersionModel

    @property
    def interfacialAreaDensityModel(self):
        return self._interfacialAreaDensityModel

    @interfacialAreaDensityModel.setter
    def interfacialAreaDensityModel(self, interfacialAreaDensityModel) -> None:
        check_type("interfacialAreaDensityModel", interfacialAreaDensityModel, InterfacialAreaDensityModel, none_ok=True)
        self._interfacialAreaDensityModel = interfacialAreaDensityModel

    @property
    def contactPartitionModel(self):
        return self._contactPartitionModel

    @contactPartitionModel.setter
    def contactPartitionModel(self, contactPartitionModel) -> None:
        check_type("contactPartitionModel", contactPartitionModel, ContactPartitionModel, none_ok=True)
        self._contactPartitionModel = contactPartitionModel