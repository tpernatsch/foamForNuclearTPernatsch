from foamForNuclear.common import copyFolder
from .mesh import Mesh


class UnvMesh(Mesh):
    """
    Copy a `.unv` file in the current case. Overwrite if a file already exist in
    the location.

    Parameters
    ----------
    srcpath : str
        Path of the `.unv` file
    """

    def __init__(self, region="", srcpath: str=""):
        super().__init__(region=region)

        self.srcpath = srcpath


    def export_to_openfoam(self):
        Mesh.export_to_openfoam(self)

        copyFolder(self.srcpath, f"constant/{self.region}")
