from foamForNuclear.checkvalue import check_value
from foamForNuclear.common import *
from foamForNuclear.openfoamFile import OpenFOAMFile


_TRANSPORT_MODEL_TYPES = {"Newtonian"}


class TransportProperties(OpenFOAMFile):
    """
    Transport properties object.

    Parameters
    ----------
    region : str
        Name of the region (default `""`).
    diffusionCoefficient : float
        Diffusion coefficient in m2/s (default `None`).
    transportModel : str {"Newtonian"}
        Transport model (default `None`).
    nu : float
    ext : str
        To be used if multiple phase.
    """
    def __init__(
            self,
            region: str="",
            diffusionCoefficient: float=None,
            transportModel: str=None,
            nu: float=None,
            ext: str=""
        ):
        super().__init__("transportProperties", "constant", region, ext)

        self.diffusionCoefficient = diffusionCoefficient
        self.transportModel = transportModel
        self.nu = nu


    @property
    def diffusionCoefficient(self):
        return self._diffusionCoefficient

    @diffusionCoefficient.setter
    def diffusionCoefficient(self, diffusionCoefficient) -> None:
        check_type("diffusionCoefficient", diffusionCoefficient, (float, int), none_ok=True)
        self._diffusionCoefficient = diffusionCoefficient

    @property
    def transportModel(self):
        return self._transportModel

    @transportModel.setter
    def transportModel(self, transportModel) -> None:
        check_type("transportModel", transportModel, str, none_ok=True)
        if (transportModel is not None):
            check_value("transportModel", transportModel, _TRANSPORT_MODEL_TYPES)
        self._transportModel = transportModel

    @property
    def nu(self):
        return self._nu

    @nu.setter
    def nu(self, nu) -> None:
        check_type("nu", nu, (float, int), none_ok=True)
        self._nu = nu


    def import_from_openfoam(self):
        foamFile = super().import_from_openfoam()
        attributes = self.get_attributes_as_list()
        for attr in attributes:
            if (attr in foamFile.keys()):
                setattr(self, attr, foamFile[attr])


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""

        if (self.diffusionCoefficient is not None):
            text += addParameter("DT", self.diffusionCoefficient, isAddExtraLine=True)
        if (self.transportModel is not None):
            text += addParameter("transportModel", self.transportModel, isAddExtraLine=True)
        if (self.nu is not None):
            text += addParameter("nu", self.nu, isAddExtraLine=True)

        return(text)
