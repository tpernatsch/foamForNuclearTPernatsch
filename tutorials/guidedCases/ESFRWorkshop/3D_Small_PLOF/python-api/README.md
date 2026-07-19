# 3D Small ESFR — Python API Workshop Cases

This folder contains Python API versions of the ESFR workshop tutorial cases.
Each case produces the same OpenFOAM simulation as its counterpart in the
parent directory, but is set up entirely through Python — no manual editing
of OpenFOAM dictionaries required.

---

## Workshop Progression

| Case | Physics | Description |
|------|---------|-------------|
| `steadyState_neutronics/` | Neutronics only | Eigenvalue calculation — find k_eff and flux shape |
| `steadyState_TH/` | TH only | Sodium flow and heat transfer without nuclear power |
| `steadyState_THAndNeutronics/` | TH + Neutronics | Coupled: power feeds TH, temperatures feed neutronics |
| *(TM cases — coming soon)* | TH + N + TM | Full coupling including core mechanical deformations |

Each case has a **completed** version and a **`_toComplete`** version (for
students) with key parameters left as `???` to fill in.

---

## Case Structure

Each case folder contains four Python files:

```
case_name/
├── mesh.py               # Load or generate the computational mesh(es)
├── boundaryConditions.py # Define initial conditions and boundary conditions
├── caseSetup.py          # Configure solvers, coupling, and run settings
└── case.py               # Main script — runs everything end-to-end
```

Run a case with:
```bash
cd steadyState_neutronics    # (or whichever case)
python case.py
```

---

## Mesh Options

Each `mesh.py` supports two approaches (controlled by `USE_PREMADE_MESH`):

**Option A — Pre-generated mesh** (`USE_PREMADE_MESH = True`, *recommended*)
> Loads existing polyMesh files from the `meshes/` folder.  Fast — no extra
> dependencies.  Use this for the workshop.

**Option B — Generate mesh** (`USE_PREMADE_MESH = False`)
> Runs the `blockMesh` meshing algorithm from scratch.  Takes ~5–10 minutes
> and requires the `honeycomb` meshing tool.  Use this if you want to change
> the core geometry.

Pre-generated meshes are in:
```
meshes/
├── polyMeshNeutro/    neutronics mesh (one cell per assembly per axial layer)
├── polyMeshFluid/     fluid mesh (with inter-assembly baffles)
└── polyMeshMeca/      thermomechanics mesh (available for future TM cases)
```

---

## Existing Full Case

A complete combined script for all physics (including thermomechanics) is
available in `Allrun.py` at this directory level.  The workshop cases in the
subfolders split that script into modular, pedagogically structured files.
