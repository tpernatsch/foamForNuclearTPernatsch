from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import OpenFOAMDict


_FIELD_OPERATION_TYPES = {"min", "max"}
_CRITERION_TYPES = {"valueAboveThreshold", "valueBelowThresholdAreAvailable"}


class PowerOffCriterionModel(OpenFOAMDict):
    """
    Used to turn off structure power sources (of all powerModels in all
    cellZones) if a certain criterion is met. Power sources can't currently be
    turned on again after they are turned off. This feature was added
    specifically to model pin bundle boiling experiments
    """
    def __init__(
            self,
            type
        ):
        super().__init__(name='powerOffCriterionModel')
        self.type = type

    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        self._type_ = type
        self.__setitem__("type", type)


class TimerPowerOffCriterionModel(PowerOffCriterionModel):
    """
    Used to turn off structure power sources (of all powerModels in all
    cellZones) if a certain criterion is met. Power sources can't currently be
    turned on again after they are turned off. This feature was added
    specifically to model pin bundle boiling experiments
    """
    def __init__(
            self,
            time: float
        ):
        super().__init__("timer")
        self.time = time

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, time) -> None:
        check_type("time", time, (float, int))
        self._time = time
        self.__setitem__("time", time)


class FieldValuePowerOffCriterionModel(PowerOffCriterionModel):
    """
    Used to turn off structure power sources (of all powerModels in all
    cellZones) if a certain criterion is met. Power sources can't currently be
    turned on again after they are turned off. This feature was added
    specifically to model pin bundle boiling experiments

    Parameters
    ----------
    fieldName : str
        Name of the volScalarField used for the criterion evaluation.
    fieldOperation : str
        Field operation performed to obtain value for comparison against
        a threshold. Currently, only `min` and `max` are available
    criterion : str
        Criterion to be applied. Currently, only `valueAboveThreshold` and
        `valueBelowThresholdAreAvailable`
    threshold : float
    timeDelay : float
        Time delay between the satifaction of the criterion and the moment
        power is actually turned off (default `0`).
    """
    def __init__(
            self,
            fieldName: str,
            fieldOperation: str,
            criterion: str,
            threshold: float,
            timeDelay: float=0
        ):
        super().__init__("fieldValue")
        self.fieldName = fieldName
        self.fieldOperation = fieldOperation
        self.criterion = criterion
        self.threshold = threshold
        self.timeDelay = timeDelay

    @property
    def fieldName(self):
        return self._fieldName

    @fieldName.setter
    def fieldName(self, fieldName) -> None:
        check_type("fieldName", fieldName, str)
        self._fieldName = fieldName
        self.__setitem__("fieldName", fieldName)

    @property
    def fieldOperation(self):
        return self._fieldOperation

    @fieldOperation.setter
    def fieldOperation(self, fieldOperation) -> None:
        check_type("fieldOperation", fieldOperation, str)
        check_value("fieldOperation", fieldOperation, _FIELD_OPERATION_TYPES)
        self._fieldOperation = fieldOperation
        self.__setitem__("fieldOperation", fieldOperation)

    @property
    def criterion(self):
        return self._criterion

    @criterion.setter
    def criterion(self, criterion) -> None:
        check_type("criterion", criterion, str)
        check_value("criterion", criterion, _CRITERION_TYPES)
        self._criterion = criterion
        self.__setitem__("criterion", criterion)

    @property
    def threshold(self):
        return self._threshold

    @threshold.setter
    def threshold(self, threshold) -> None:
        check_type("threshold", threshold, (float, int))
        self._threshold = threshold
        self.__setitem__("threshold", threshold)

    @property
    def timeDelay(self):
        return self._timeDelay

    @timeDelay.setter
    def timeDelay(self, timeDelay) -> None:
        check_type("timeDelay", timeDelay, (float, int))
        self._timeDelay = timeDelay
        self.__setitem__("timeDelay", timeDelay)
