from abc import abstractmethod
import os

import foamlib

from foamForNuclear.common import *
from foamForNuclear.checkvalue import check_type


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
    def __init__(self, name, folder: str="", region: str="", ext: str=""):
        self.name = name
        self.folder = folder
        self.region = region
        self.ext = ext

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, name: str) -> None:
        check_type("name", name, str, none_ok=True)
        self._name = name

    @property
    def folder(self) -> str:
        return self._folder

    @folder.setter
    def folder(self, folder: str) -> None:
        if folder is not None:
            check_type("folder", folder, str)
            self._folder = folder
        else:
            self._folder = ''

    @property
    def region(self) -> str:
        return self._region

    @region.setter
    def region(self, region: str) -> None:
        if region is not None:
            check_type("region", region, str)
            self._region = region
        else:
            self._region = ''

    @property
    def ext(self) -> str:
        return self._ext

    @ext.setter
    def ext(self, ext: str) -> None:
        if ext is not None:
            check_type("ext", ext, str)
            self._ext = ext
        else:
            self._ext = ''

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
