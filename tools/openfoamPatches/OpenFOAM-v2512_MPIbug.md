# OpenFOAM v2512 Bug: MPI_ERR_TRUNCATE in Parallel FSI Simulations

## Summary

A bug in `mappedPatchBase::upToDate()` causes `MPI_ERR_TRUNCATE` crashes during parallel FSI
(Fluid-Structure Interaction) simulations when using `defaultCommsType = nonBlocking`. The root
cause is that per-processor local booleans are used as conditions for entering blocks that contain
collective MPI operations, leading to misaligned `MPI_Allreduce` calls across processors.

---

## Environment

- **OpenFOAM version**: v2512
- **Solver**: GeN-Foam (FSI)
- **MPI comms type**: `nonBlocking`
- **Parallel decomposition**: 8 processors
- **Crash point**: `Time = 6.25e-07` (reproducible)

---

## Symptoms

The simulation crashes with the following MPI error on one or more ranks:

```
[rankN] Fatal error in PMPI_Allreduce: MPI_ERR_TRUNCATE: message truncated
```

The crash occurs at a specific timestep and is fully reproducible. Different ranks may report the
error depending on the processor decomposition (in the observed case, ranks 1 and 7).

---

## Root Cause

**File**: `OpenFOAM-v2512/src/meshTools/mappedPatches/mappedPolyPatch/mappedPatchBase.C`  
**Function**: `mappedPatchBase::upToDate()`

The function computes two local (per-processor) booleans:

```cpp
bool thisUpToDate   = thisMesh.upToDatePoints(updateMeshTime());
bool sampleUpToDate = sampleMesh().upToDatePoints(updateSampleMeshTime());
```

It then conditionally calls a lambda `checkPointMovement()`, which internally calls
`Pstream::reduceAnd()` — a **collective MPI operation** that requires all processors to
participate simultaneously:

```cpp
if (!thisUpToDate && thisMesh.moving())
    thisUpToDate = checkPointMovement(thisMesh, patch_, updateMeshTime());

if (!sampleUpToDate && sampleMesh().moving() && ...)
    sampleUpToDate = checkPointMovement(sampleMesh(), samplePolyPatch(), updateSampleMeshTime());
```

Because `thisUpToDate` and `sampleUpToDate` are evaluated **locally** on each processor, they can
differ across ranks. This happens in FSI cases where the fluid patch (`outerWall`) has **zero
faces** on some processors (due to unequal decomposition). Those processors never called
`setUpToDate()` via `calcMapping()` / `calcAMI()` for the sample mesh, so their
`updateSampleMeshTime_` event number lags behind the solid mesh's `points_` event number.

**Outcome**: Some processors enter the `checkPointMovement()` block (calling `reduceAnd()`), while
others skip it. This misaligns the collective operations:

| Processor | `sampleUpToDate` | `reduceAnd` calls inside `upToDate()` |
|-----------|-----------------|---------------------------------------|
| 0, 2, 3, 4, 5, 6 | `true`  | 1 (only `thisUpToDate` path)          |
| 1, 7              | `false` | 2 (both `thisUpToDate` and `sampleUpToDate` paths) |

The MPI library detects the mismatch in collective call counts → `MPI_ERR_TRUNCATE`.

### Why `sampleUpToDate` diverges

`calcMapping()` calls `setUpToDate()` on both `updateMeshTime_` and `updateSampleMeshTime_`,
keeping them fresh. However, `calcAMI()` does **not** call `setUpToDate()` on either. After the
first timestep where only `calcAMI()` is invoked (no full remapping), processors that happen to
have zero `outerWall` faces never refresh `updateSampleMeshTime_`, causing it to fall behind the
solid mesh's moving-points event counter.

---

## Fix

The fix ensures that all processors agree on `thisUpToDate` and `sampleUpToDate` **before** any
processor decides whether to enter the collective `checkPointMovement()` block. Two
`Pstream::reduceAnd()` synchronization calls are added as guards:

```cpp
// SYNCHRONIZATION FIX: ensure all procs enter/skip each collective together
// (prevents MPI_ERR_TRUNCATE with nonBlocking comms when sampleUpToDate
//  diverges locally across procs with unequal patch face distributions)
if (thisMesh.moving())
    Pstream::reduceAnd(thisUpToDate);
if (!thisUpToDate && thisMesh.moving())
    thisUpToDate = checkPointMovement(thisMesh, patch_, updateMeshTime());

if (sameWorld() && sampleMesh().moving()
 && (mode_==NEARESTPATCHFACE || mode_==NEARESTPATCHFACEAMI || mode_==NEARESTPATCHPOINT))
    Pstream::reduceAnd(sampleUpToDate);
if (!sampleUpToDate && sampleMesh().moving()
 && (mode_==NEARESTPATCHFACE || mode_==NEARESTPATCHFACEAMI || mode_==NEARESTPATCHPOINT))
    sampleUpToDate = checkPointMovement(sampleMesh(), samplePolyPatch(), updateSampleMeshTime());
```

The guard `reduceAnd` uses the minimum of the local boolean across all processors (logical AND),
so if **any** processor considers the mesh out-of-date, all processors agree it is out-of-date and
all enter the `checkPointMovement()` block together. This keeps collective calls symmetric.

### Complete fixed `upToDate()` function

```cpp
bool Foam::mappedPatchBase::upToDate() const
{
    const polyMesh& thisMesh = patch_.boundaryMesh().mesh();

    bool thisUpToDate = thisMesh.upToDatePoints(updateMeshTime());
    bool sampleUpToDate =
    (
        sameWorld()
      ? sampleMesh().upToDatePoints(updateSampleMeshTime())
      : true
    );

    auto checkPointMovement = []
    (
        const polyMesh& mesh,
        const polyPatch& patch,
        regIOobject& state
    ) -> bool
    {
        bool upToDate = true;
        const auto& oldPoints = mesh.oldPoints();
        const auto& points = mesh.points();
        for (const label pointi : patch.meshPoints())
        {
            if (mag(oldPoints[pointi] - points[pointi]) > SMALL)
            { upToDate = false; break; }
        }
        Pstream::reduceAnd(upToDate);
        if (upToDate) { state.setUpToDate(); }
        return upToDate;
    };

    // SYNCHRONIZATION FIX: ensure all procs enter/skip each collective together
    if (thisMesh.moving())
        Pstream::reduceAnd(thisUpToDate);
    if (!thisUpToDate && thisMesh.moving())
        thisUpToDate = checkPointMovement(thisMesh, patch_, updateMeshTime());

    if (sameWorld() && sampleMesh().moving()
     && (mode_==NEARESTPATCHFACE || mode_==NEARESTPATCHFACEAMI || mode_==NEARESTPATCHPOINT))
        Pstream::reduceAnd(sampleUpToDate);
    if (!sampleUpToDate && sampleMesh().moving()
     && (mode_==NEARESTPATCHFACE || mode_==NEARESTPATCHFACEAMI || mode_==NEARESTPATCHPOINT))
        sampleUpToDate = checkPointMovement(sampleMesh(), samplePolyPatch(), updateSampleMeshTime());

    return (thisUpToDate && sampleUpToDate);
}
```

---

## Verification

After rebuilding `libmeshTools.so` with the fix (`wmake OpenFOAM-v2512/src/meshTools`), the
simulation:

- Passed `Time = 6.25e-07` — the previously reliable crash point — without any MPI error
- Continued cleanly to `Time = 1.55e-05` (25× beyond the crash point, ~120 timesteps) with zero
  `MPI_ERR`, `TRUNCATE`, or `Fatal` messages

---

## Affected Configurations

This bug is most likely to manifest when **all** of the following conditions are true:

1. Parallel run with `defaultCommsType = nonBlocking`
2. FSI simulation with a moving mesh (solid and fluid regions)
3. The fluid patch mapped to the solid patch has **zero faces** on some processors (unequal
   decomposition)
4. `calcAMI()` is called without a subsequent `calcMapping()` after the first timestep

---

## Files Modified

| File | Change |
|------|--------|
| `OpenFOAM-v2512/src/meshTools/mappedPatches/mappedPolyPatch/mappedPatchBase.C` | Added two `Pstream::reduceAnd()` guard calls in `upToDate()` |
