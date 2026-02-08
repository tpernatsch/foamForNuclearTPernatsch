import foamlib
from foamForNuclear.common import *
from foamForNuclear.openfoamFile import OpenFOAMFile


class fvSchemes(OpenFOAMFile):
    """
    Scheme selection. More info in the OpenFOAM website
    https://www.openfoam.com/documentation/user-guide/6-solving/6.2-numerical-schemes

    Attributes
    ----------
    d2dt2Schemes : OpenFOAMDict
    ddtSchemes : OpenFOAMDict
    gradSchemes : OpenFOAMDict
    divSchemes : OpenFOAMDict
    laplacianSchemes : OpenFOAMDict
    interpolationSchemes : OpenFOAMDict
    snGradSchemes : OpenFOAMDict
    fluxRequired : OpenFOAMDict
    """

    def __init__(self, region=""):
        super().__init__(name="fvSchemes", folder="system", region=region)

        self.d2dt2Schemes: OpenFOAMDict = OpenFOAMDict({'default': 'Euler'}, name="d2dt2Schemes")
        self.ddtSchemes: OpenFOAMDict = OpenFOAMDict({'default': 'Euler'}, name="ddtSchemes")
        self.gradSchemes: OpenFOAMDict = OpenFOAMDict({'default': 'Gauss linear'}, name="gradSchemes")
        self.divSchemes: OpenFOAMDict = OpenFOAMDict({'default': 'Gauss linear'}, name="divSchemes")
        self.laplacianSchemes: OpenFOAMDict = OpenFOAMDict({'default': 'Gauss linear corrected'}, name="laplacianSchemes")
        self.interpolationSchemes: OpenFOAMDict = OpenFOAMDict({'default': 'linear'}, name="interpolationSchemes")
        self.snGradSchemes: OpenFOAMDict = OpenFOAMDict({'default': 'corrected'}, name="snGradSchemes")
        self.fluxRequired: OpenFOAMDict = OpenFOAMDict({'default': 'false'}, name="fluxRequired")


    def import_from_openfoam(self):
        foamFile = super().import_from_openfoam()
        attributes = self.get_attributes_as_list()
        for attr in attributes:
            if (
                attr in foamFile.keys()
                and isinstance(foamFile[attr], foamlib._files._files.FoamFile.SubDict)
            ):
                setattr(self, attr, OpenFOAMDict(foamFile[attr], name=attr))
                

    def _format_scheme_dict(self, name: str, mapping: dict) -> str:
        lines: list[str] = []
        lines.append(name)
        lines.append("{")

        for key, val in mapping.items():
            # If the value is a tuple/list, join elements with spaces
            if isinstance(val, (tuple, list)):
                val_str = " ".join(str(x) for x in val)
            else:
                val_str = str(val)
            if( key == "default"):
                lines.append(f'    {key}'.ljust(30) + f" {val_str};")
            else:
                lines.append(f'    "{key}"'.ljust(30) + f" {val_str};")


        lines.append("}")
        lines.append("")  # blank line after block
        return "\n".join(lines)


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        parts = [
           self._format_scheme_dict("d2dt2Schemes",        self.d2dt2Schemes),
           self._format_scheme_dict("ddtSchemes",          self.ddtSchemes),
           self._format_scheme_dict("gradSchemes",         self.gradSchemes),
           self._format_scheme_dict("divSchemes",          self.divSchemes),
           self._format_scheme_dict("laplacianSchemes",    self.laplacianSchemes),
           self._format_scheme_dict("interpolationSchemes", self.interpolationSchemes),
           self._format_scheme_dict("snGradSchemes",       self.snGradSchemes),
           self._format_scheme_dict("fluxRequired",        self.fluxRequired),
        ]
        return "\n".join(parts)

