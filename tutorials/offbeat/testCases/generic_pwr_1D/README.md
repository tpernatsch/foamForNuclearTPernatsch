# **generic_pwr_1D — OFFBEAT/OpenFOAM case**

## **Overview**

This case represents a generic 1D PWR fuel rod simulation (axisymmetric r–z formulation with a coarse axial discretization). It is intended as a reference case to exercise the coupled thermal, mechanics, gap, neutronics/burnup, and fission gas release models in a typical LWR operating scenario with a time-dependent power history.

The model includes:

- UO₂ fuel and Zircaloy cladding
- a fuel–clad gap with gap heat transfer, gap gas pressure, and contact
- time-dependent linear heat generation rate (LHGR) with an axial shape
- fast flux history
- burnup evolution and radial power shaping from burnup
- fission gas release via an external model (SCIANTIX interface)

This is not a plant-specific case: geometry, histories, and models are representative and chosen for robustness and demonstration purposes.

---

## **Geometry and mesh**

The rod is modeled as two concentric solid regions:

- Fuel pellet: inner radius `fuel_ri = 0` to outer radius `fuel_ro = 4.5 mm`
- Cladding: inner radius `clad_ri = 4.565 mm` to outer radius `clad_ro = 5.315 mm`

Axial extent:

- Active fuel length: `fuel_length = 3.0 m`
- Plenum length: `plenum_length = 0.2 m`

Discretization:

- Fuel radial cells: `fuel_nr = 30`
- Cladding radial cells: `clad_nr = 10`
- Axial slices: `slices = 10` (uniform segmentation of the rod length)

???+ note

    The axial discretization is intentionally coarse to mirror traditional 1.5D fuel
    performance codes and to keep runtime low. The model is axisymmetric; only the
    radial coordinate is resolved within each axial slice.

---

## **Physics included**

### **Thermal**

- Solid heat conduction in fuel and cladding.
- Boundary conditions:

    - Fuel outer surface: gap heat transfer condition (fuel side)
    - Clad inner surface: gap heat transfer condition (clad side)
    - Clad outer surface: imposed cladding temperature history (tabulated in time)

The initial temperature is uniform:

- `T0 = 300 K`

The cladding outer temperature is ramped to hot conditions:

- `(t = 0 h → 300 K)`, `(t = 1 h → 600 K)`, then held at `600 K` until end time

### **Mechanics**

- Small-strain solid mechanics in fuel and cladding.
- Fuel–clad mechanical interaction handled via a gap/contact interface on:

    - `fuelOuter`
    - `cladInner`

Clad outer mechanical loading is applied via coolant pressure:

- a tabulated pressure history from near-atmospheric to high pressure, then held

Fuel models:

- Relocation (FRAPCON-type) based on cold gap and cold fuel diameter, with a recovery fraction
- Densification (FRAPCON-type) with a 1% resintering density change
- Isotropic cracking damage model

Cladding models:

- von Mises plasticity with constant yield stress and a creep model
- Yield stress is set large relative to this irradiation scenario, while creep provides the dominant time-dependent deformation mechanism

### **Gap gas**

- Gap gas pressure is initialized to a single value and evolves through the gap gas model
- Gap patches: `fuelOuter`, `cladInner`

### **Neutronics and burnup**

- A simplified diffusion-based neutronics model is enabled to support burnup- and power-related quantities
- Burnup is computed using a Lassmann- or TUBRNP-inspired model with a specified convergence tolerance

### **Heat source (LHGR)**

The volumetric heat source is driven by a time-dependent LHGR history:

- Time points (hours): `0 → 1 → 35000`
- LHGR (W/m): `0 → 20 kW/m → 10 kW/m`

Spatial shaping:

- Axial profile: time-dependent tabulated shape, normalized and applied along the axial direction
- Radial profile: derived from burnup (`FromBurnup`), enabling burnup-dependent radial peaking

### **Fast flux**

A time-dependent fast flux history is applied to fuel and cladding:

- `0 → 2e13 → 1e13` (units as defined by the case setup)

This drives irradiation-dependent models such as cladding creep.

### **Fission gas release (FGR)**

Fission gas behavior is computed through the SCIANTIX interface with a selected set of enabled sub-models (grain growth, diffusion, bubble evolution, etc.). The FGR output is collected through a dedicated function object.

---

## **Materials**

### **Fuel: UO₂**

- Material name: `fuel`
- Grain radius: `rGrain = 5 µm`
- Enrichment (weight fraction):

    - U-235: 4.5%
    - U-238: 95.5%

- Included models:

    - isotropic cracking / damage
    - densification
    - relocation

### **Cladding: Zircaloy**

- Material name: `cladding`
- Constitutive behavior:

    - von Mises plasticity with irradiation creep
    - creep relaxation factor defined in the model configuration

---

## **Time control and numerical settings**

- End time: `35,000 h`
- User time unit: `hours`
- Initial `deltaT`: `0.001 h`
- Adjustable time stepping enabled with:

    - `maxDeltaT = 100 h`
    - `minDeltaT = 0.01 h`

Write-out:

- `writeControl = timeStep`
- `writeInterval = 100`

Stress analysis:

- `maxOuterIter = 100`

> The time-step controls are intentionally permissive to allow aggressive growth after the initial transient. For production-oriented simulations, tighter limits are recommended.

---

## **Boundary patches used**

Thermal:

- `fuelOuter`: gap heat transfer (fuel side)
- `cladInner`: gap heat transfer (clad side)
- `cladOuter`: prescribed temperature vs time

Mechanics:

- `fuelOuter`: gap/contact (fuel side)
- `cladInner`: gap/contact (clad side)
- `cladOuter`: coolant pressure loading vs time

---

## **Function objects and expected outputs**

The case writes a limited set of high-value diagnostics:

1. **Fuel centerline probes**

    - Fields: `T`, `Bu`
    - Location: mid-plane, centerline

2. **Fuel outer surface patch probes**

    - Fields: `T`, `Bu`, `gapWidth`, `interfaceP`, `hGap`
    - Patch: `fuelOuter`
    - Location: mid-plane

3. **Radial profile graph**

    - Line from `r = fuel_ri` to `r = fuel_ro` at a selected axial slice
    - Fields: `T`, `sigma`
    - Used to assess temperature gradients and stress distribution in the fuel

4. **Fission gas release output**

    - Collected through the FGR function object (SCIANTIX-driven)

Typical post-processing checks:

- Temperature ramp and approach to quasi-steady radial profiles under power
- Evolution of gap width and contact pressure during power increase
- Burnup accumulation and correlation with radial power shaping
- Cladding response under pressure and creep
- FGR trends with increasing burnup and temperature

---

## **What this case is for (and what it is not)**

Good for:

- Regression testing of coupled physics in a robust, low-cost configuration
- Sanity checks of time histories and boundary condition propagation
- Demonstration of typical rod-level quantities (temperature, gap conductance, burnup)

Not intended for:

- Licensing-grade predictions
- Detailed axial resolution (only 10 axial slices)
- High-fidelity coolant–cladding coupling (clad outer temperature is prescribed)

---

## **How to run**

From the case directory, using a standard OpenFOAM workflow:

1. `blockMesh`
2. `changeDictionary` (optional, to modify dictionaries)
3. `offbeat` (or `offbeat | tee log.offbeat`)
4. `gnuplot Residuals.gp` to visualize residual evolution
5. A Python post-processing script (`plot.py`) is provided as an example
