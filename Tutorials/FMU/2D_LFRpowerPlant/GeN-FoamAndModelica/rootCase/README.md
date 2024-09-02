# Root Case of the ALFRED_Critical_Wedge

---

## Salome

- Launch Salome
- Go to: `File/Load Script`
- Select `geometry/ALFRED_2D.py"`
- Remove for each mesh "Partition_<region>" the "Groups of Edges"
- Export each mesh "Partition_<region>" in UNV format with
    `File/Export/UNV file`. The name of the exported file must match the name in
    the `mesh/mesh.sh`.

## Build the geometry from extracted UNV files

### Shortcut

Execute in command line in the `rootCase` folder:
```bash
# Automatically (re)generate the polyMesh folders from the UNV files in mesh/
./mesh/mesh.sh
```

### Explanations

#### Main geometry

The geometry must be included in the `mesh` folder. To build it, the following
command must be executed in the openfoam shell.

```bash
ideasUnvToFoam Mesh_<region>.unv
```

Then apply the mesh merging through the OpenFOAM routine `createPatch`
```bash
createPatch -overwrite -region <region>
```

#### Baffles

To merge the faceZones that compose the baffles (internal walls), run:
```bash
topoSet -region <region>
```

Then create the baffles:
```bash
createBaffles -overwrite -region <region>
```
