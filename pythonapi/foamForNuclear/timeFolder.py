import os

from foamForNuclear.checkvalue import CheckedList, check_type
from foamForNuclear.field import Field
from foamForNuclear.openfoamFile import OpenFOAMFile


class TimeFolder(CheckedList):
    """
    Object containing all :class:`Field` objects. This is used as initial
    conditions.

    Parameters
    ----------
    time : str | float | int
        Time value.
    fields : list[Field | OpenFOAMFile]
        List of :class:`Field` objects.
    """
    def __init__(
            self,
            time: str | float | int=0,
            fields: list[Field | OpenFOAMFile]=None
        ):
        super().__init__((Field | OpenFOAMFile), 'Fields collection', fields)

        self.time = time


    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, time) -> None:
        check_type("time", time, (str, float, int))
        self._time = f"{time}"

    def append(self, item):
        """
        Append item to list

        Parameters
        ----------
        item : object
            Item to append
        """
        item.folder = self.time
        super().append(item)


    def export_to_openfoam(self):
        if (not os.path.exists(f"{self.time}")):
            os.mkdir(f"{self.time}")

        for field in self:
            field.folder = self.time

            if (field.region is not None and not os.path.exists(f"{self.time}/{field.region}")):
                os.mkdir(f"{self.time}/{field.region}")

            field.export_to_openfoam()


    def import_from_openfoam(self):
        for field in self:
            field.folder = self.time
            field.import_from_openfoam()
