import os

from foamForNuclear.checkvalue import CheckedList, check_type
from foamForNuclear.fields import Field
from foamForNuclear.openfoamFile import OpenFOAMFile
from foamForNuclear.mesh import Mesh


class TimeFolder(CheckedList):
    """
    Object containing all Field objects associated with a given time.
    Mainly used internally to group fields for OpenFOAM export/import.

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

        if isinstance(time, str):
            self._time = time
        elif isinstance(time, int):
            self._time = str(time)
        elif isinstance(time, float):
            self._time = str(int(time)) if time.is_integer() else f"{time:g}"

    def add_field(self, item: Field | list[Field]):
        """
        Append item to list of fields

        Parameters
        ----------
        item : object
            Item to append
        """
        if(isinstance(item, list)):
            for item_i in item:
                item_i.folder = self.time
                super().append(item_i)
        else:
            item.folder = self.time
            super().append(item)


    def export_to_openfoam(
        self,
        region_meshes: dict[str | None, Mesh] | None = None,
    ):
        if not os.path.exists(f"{self.time}"):
            os.mkdir(f"{self.time}")

        for field in self:
            field.folder = self.time

            # --- auto-prepare boundary field if we know the mesh ---
            if region_meshes is not None:
                region_name = field.region
                mesh = region_meshes.get(region_name)

                if mesh is not None and hasattr(field, "prepare_boundary_field"):
                    field.prepare_boundary_field(mesh)

            if (field.region is not None and not os.path.exists(f"{self.time}/{field.region}")):
                os.mkdir(f"{self.time}/{field.region}")

            field.export_to_openfoam()

    def import_from_openfoam(self):
        for field in self:
            field.folder = self.time
            field.import_from_openfoam()
