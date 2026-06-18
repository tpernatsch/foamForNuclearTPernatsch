# OpenFOAM patches

Fixes that must be applied to the **OpenFOAM source tree itself** (not foamForNuclear), and so
cannot be carried by this repository. Run once per OpenFOAM installation.

## `patch_OpenFOAM_mappedPatchBase`

Fixes an OpenFOAM-v2512 bug where parallel FSI runs crash with `MPI_ERR_TRUNCATE`:
`mappedPatchBase::upToDate()` uses per-processor local booleans as conditions guarding
collective MPI calls, so with an unequal decomposition (zero-face patches on some ranks) the
collectives misalign. The patch adds `Pstream::reduceAnd()` synchronisation before those
conditions and rebuilds `libmeshTools.so`. Idempotent (no-op on an already-patched tree).

```sh
./patch_OpenFOAM_mappedPatchBase
```

Full write-up: [`OpenFOAM-v2512_MPIbug.md`](OpenFOAM-v2512_MPIbug.md). This fix is also a
candidate for submission to ESI OpenFOAM upstream.
