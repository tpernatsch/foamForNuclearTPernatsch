from abc import abstractmethod
import os

import foamlib

from foamForNuclear.common import *
from foamForNuclear.checkvalue import check_type

import attrs as attr
from attrs import define, field
import attrs.validators as v
from ._attrs_tools import auto_type_validator


@define(
    slots=True,
    on_setattr=[attr.setters.convert, attr.setters.validate],
    field_transformer=auto_type_validator,
    repr=False,
)
class OpenFOAMFile:
    """
    Base class for OpenFOAM file

    Parameters
    ----------
    name : str
        Name of the file
    folder : str
        Name of the folder to place the file (default: `''`)
    region : str
        Name of the region (default: `''`)
    ext : str
        Extension at the end of the file name, e.g `".liquid"` (default: `''`)
    """
    name: str = ""
    folder: str = ""
    region: str = ""
    ext: str = ""

    @property
    def path(self):
        return(f"{self.folder}/{self.region}/{self.name}{'.'+self.ext if self.ext != '' else ''}")
        return(f"{OpenFOAMFile.caseFolder}/{self.folder}/{self.region}/{self.name}{'.'+self.ext if self.ext != '' else ''}")


    def get_attributes_as_list(self):
        return([
            a for a in dir(self)
            if not a.startswith('__') and not callable(getattr(self, a))
        ])


    def _write_to_file(func) -> str:
        def wrapper(self, *args, **kwargs) -> str:
            if (not os.path.exists(self.folder)):
                os.mkdir(self.folder)
            if (self.region != '' and not os.path.exists(f"{self.folder}/{self.region}")):
                os.mkdir(f"{self.folder}/{self.region}")
            # if (not os.path.exists(f"{OpenFOAMFile.caseFolder}/{self.folder}")):
            #     os.mkdir(f"{OpenFOAMFile.caseFolder}/{self.folder}")
            # if (self.region != '' and not os.path.exists(f"{OpenFOAMFile.caseFolder}/{self.folder}/{self.region}")):
            #     os.mkdir(f"{OpenFOAMFile.caseFolder}/{self.folder}/{self.region}")

            text  = openfoamHeader
            text += openfoamFileHeader(self.name)
            text += func(self, *args, **kwargs)

            text += openfoamFooterLine

            with open(self.path, 'w') as f:
                f.write(text)

            return(text)

        return wrapper

    @abstractmethod
    def import_from_openfoam(self):
        """
        Read OpenFOAM file.

        Returns
        -------
            Return the file read through `foamLib`.
        """
        return foamlib.FoamFile(self.path)


    @abstractmethod
    @_write_to_file
    def export_to_openfoam(self):
        """
        Write current file into OpenFOAM format.
        """
