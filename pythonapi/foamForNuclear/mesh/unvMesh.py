from foamForNuclear.common import copyFolder
from .mesh import Mesh

from pathlib import Path
from .blockMesh import Face

import os

class UnvMesh(Mesh):
    """
    Import a `.unv` mesh into the case and allow declaring patch types
    (wedge/empty/symmetry/cyclic...) for auto BC defaults.
    """

    def __init__(
        self,
        region: str = "",
        srcpath: str = "",
        *,
        special_patches: dict[str, object] | None = None,
    ):
        super().__init__(region=region)

        # resolve path early (your current approach)
        creation_cwd = Path(os.getcwd())
        p = Path(srcpath)
        if not p.is_absolute():
            p = (creation_cwd / p).resolve()
        self.srcpath = str(p)

        self.faces: list[Face] = []
        self.special_patches = special_patches or {}

        # store as Face objects (subset only)
        for name, spec in self.special_patches.items():
            if isinstance(spec, tuple):
                btype, extra = spec
                self.add_boundary(Face(name, boundaryType=btype, extraParameters=dict(extra)))
            else:
                self.add_boundary(Face(name, boundaryType=str(spec)))

    def add_boundary(self, face: Face):
        # mirror BlockMesh API
        self.faces.append(face)

    def export_to_openfoam(self):
        Mesh.export_to_openfoam(self)
        copyFolder(self.srcpath, f"constant/{self.region}")

