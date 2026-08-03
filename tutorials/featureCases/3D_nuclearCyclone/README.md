# 3D Nuclear Cyclone

Tags: [![badge](https://img.shields.io/badge/openFoamImported-nuclearParcelFoam-blue.svg)]()

## Description

This case exercises the `nuclearParcelFoam` solver (`applications/modules/openFoamImported/nuclearParcelFoam`)
on the geometry of the stock OpenFOAM `lagrangian/MPPICFoam/cyclone` tutorial.

`nuclearParcelFoam` is built on a plain incompressible `KinematicCloud` (as used by
`kinematicParcelFoam`), not on the multiphase/MPPIC machinery of `MPPICFoam` — so
this case is *not* a like-for-like reproduction of the stock cyclone tutorial. It
reuses the same cyclone geometry (mesh, inlet/outlet patches) purely as a
convenient, non-trivial 3D flow domain to check that the solver, its buoyant
Boussinesq energy equation, and the `FPCloud` fission-product cloud (drag,
gravity, thermophoresis, wall rebound/escape, convective heat transfer and decay
heat) all run together correctly through foamForNuclear's `GeN-Foam` driver.

The continuous phase is laminar air entering through the `inlet` patch and
leaving through `outlet`. Nuclear particles ("fission product" parcels) are
injected at the inlet, tracked through the swirling cyclone flow, decay
(releasing heat back to the gas via `constDecay`/`RanzMarshall`), and either
escape through the outlet or rebound off the walls.

The domain starts hot: both the `walls` patch and the initial internal field
are held at 330 K (`0.orig/fluidRegion/T`), well above the 293 K inlet flow.
The cyclone is then progressively cooled from the inlet jet inwards, giving a
real, evolving temperature gradient for the `Thermophoretic` particle force to
respond to: particles should be visibly pushed away from the hotter
(as-yet-uncooled) regions towards the cooler inlet-fed core as the flow
develops.

## Running

```
./Allrun
```

This fetches the cyclone STL geometry from `$FOAM_TUTORIALS/resources/geometry`,
runs `blockMesh`/`snappyHexMesh` for the `fluidRegion`, then runs `GeN-Foam`.
