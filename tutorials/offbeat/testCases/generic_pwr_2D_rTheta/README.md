# **generic_pwr_2D_rTheta — OFFBEAT/OpenFOAM case**

???+ warning

    This case is the 2D r–θ counterpart of the `generic_pwr_1D` tutorial.
    Most of the physical models, boundary conditions, and operating
    assumptions are identical. Only the geometry, mesh, and mechanical
    formulation differ and are highlighted explicitly below.

## **Overview**

This case represents a generic 2D PWR fuel rod cross-section simulation using an r–θ formulation. It is intended as a reference example to demonstrate the same coupled thermal, mechanics, gap, neutronics/burnup, and fission gas release models used in the 1D case, but applied to a two-dimensional pellet disc.

The main purpose of this case is to illustrate:
- the use of externally generated meshes,
- and the impact of a different mechanical assumption compared to the 1D formulation.

Unless otherwise stated, the physical models and time histories are the same as in the `generic_pwr_1D` case.

---

## **Geometry and mesh**

The geometry consists of a single **pellet disc** representing a transverse cross-section of the fuel rod.

Key characteristics:

- 2D r–θ geometry (no axial resolution)
- Single fuel disc, without explicit axial extent
- Circumferential discretization enabled

Mesh generation:

- The mesh is generated externally using **Salome**
- Mesh file: `PelletDiscCoarser.unv`
- Imported using `ideasUnvToFoam`

Mesh preprocessing steps:

- Scaling from meters to millimeters
- Translation to align the disc with the mid-height of the rod

Compared to the 1D case, the fuel and cladding radii are slightly different. 
This is reflected directly in the parameters used by relocation model (`GapCold` and `DiamCold` in the `relocation` subdictionary)

---

## **Physics included**

The set of physics models is identical to the `generic_pwr_1D` case unless explicitly noted below.

### **Thermal**

Thermal conduction is solved in the fuel and cladding over the 2D r–θ domain.

Boundary conditions and thermal histories are the same as in the 1D case.

### **Mechanics**

Solid mechanics are solved in plane stress over the pellet disc.

Key differences with respect to the 1D case:

- Plane stress formulation is used instead of modified plane strain
- Circumferential stress variations are explicitly resolved

Fuel–clad mechanical interaction is still handled through a gap/contact interface.

Fuel models:

- Relocation model identical in form to the 1D case, but with updated geometric parameters:

    - `GapCold` reflects the slightly different rod radius
    - `DiamCold` reflects the updated fuel diameter

- Densification model with reduced resintering density change:

    - resintering density change: 0.3 (instead of 1.0)

Cladding models:

- von Mises plasticity with irradiation creep, as in the 1D case

### **Gap gas**

Gap gas pressure and evolution are treated identically to the 1D case.

### **Neutronics and burnup**

Neutronics and burnup models are unchanged with respect to the 1D case and are used consistently to drive power and burnup-dependent quantities.

### **Heat source (LHGR)**

The same time-dependent LHGR history and axial power normalization used in the 1D case are applied here.

In this 2D configuration, the power is applied uniformly over the pellet disc area, with radial shaping still derived from burnup.

### **Fast flux**

The same time-dependent fast flux history is applied to fuel and cladding and drives irradiation-dependent material behavior.

### **Fission gas release (FGR)**

Fission gas release is computed through the same SCIANTIX interface and configuration used in the 1D case.

---

## **Materials**

### **Fuel: UO₂**

- Material name: `fuel`
- Grain radius: same as in the 1D case
- Enrichment: identical to the 1D case

Differences with respect to the 1D case:

- Reduced densification amplitude
- Updated geometric parameters used by relocation due to the different disc radius

### **Cladding: Zircaloy**

- Material name: `cladding`
- Constitutive behavior identical to the 1D case

---

## **Time control and numerical settings**

Time integration settings, time-step control, and write intervals are the same as in the `generic_pwr_1D` case.

This allows a direct comparison between 1D and 2D results under identical operating conditions.

---

## **Boundary patches used**

Boundary patches follow the same naming convention as the 1D case, adapted to the 2D mesh topology.

Thermal and mechanical boundary conditions are applied consistently with the pellet disc geometry.

---

## **Function objects and expected outputs**

The same class of diagnostic quantities as in the 1D case is available, adapted to the 2D geometry.

Typical outputs include:

- Temperature and stress fields over the pellet disc
- Circumferential stress variations
- Gap-related quantities when contact occurs
- Burnup and fission gas release indicators

This case is particularly useful to assess azimuthal variations that cannot be captured in a 1D formulation.

---

## **What this case is for (and what it is not)**

Good for:

- Demonstrating 2D r–θ thermo-mechanics in OFFBEAT
- Assessing circumferential stress and temperature variations
- Comparing plane stress and modified plane strain assumptions
- Verifying consistency between 1D and 2D formulations

Not intended for:

- Full 3D fuel rod simulations
- Licensing-grade predictions
- Detailed axial behavior or end effects

---

## **How to run**

From the case directory, using a standard OpenFOAM workflow:

1. `ideasUnvToFoam PelletDiscCoarser.unv`
2. `transformPoints -scale "(0.001 0.001 0.001)"` if using OpenFOAM-ESI, alternatively
    `transformPoints "scale=(0.001 0.001 0.001)"` if using OpenFOAM-Foundation
3. `transformPoints -translate "(0.0 0.0 1.5)"` if using OpenFOAM-ESI, alternatively
    `transformPoints "translate=(0.0 0.0 1.5)"` if using OpenFOAM-Foundation    
4. `changeDictionary` as all patches created by `ideasUnvToFoam` are of the `patch` type
5. `offbeat` (or `offbeat | tee log.offbeat`)
6. Post-process fields and function object outputs as in the 1D case