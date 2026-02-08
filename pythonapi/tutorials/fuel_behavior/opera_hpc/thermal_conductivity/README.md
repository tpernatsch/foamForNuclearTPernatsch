# Exercise 1 — OFFBEAT (Thermal only)

## Description

In the next few exercises, we study a short rodlet (10 pellets) over a simplified **one-year** irradiation. We start with basics and **thermal analysis only**: no relocation, densification, or creep-down. Think of the fuel as a hot cylinder that releases heat to the Zircaloy cladding across a **fixed-width gap** with **fixed thermal contact conductance**.

---

## Data

### Geometry

* No central hole
* Fuel outer radius: **4.5 mm**
* Gap: **60 µm = 0.06 mm**
* Cladding radial width: **700 µm = 0.7 mm**
* Fuel column: **10** pellets × **11 mm** = **110 mm**
* Top plenum: **10 mm**
* Ignore end caps

### Discretization

* Fuel: **30** cells radially, **1** axially
* Cladding: **10** cells radially, **1** axially

### Temperature / BCs

* Initial temperature: ( T_0 = 293,\mathrm{K} )
* Cladding outer temperature: ( T_{c0} = 573,\mathrm{K} )
* Fixed gap conductance: ( h_{\mathrm{gap}} = 5000,\mathrm{W/(m^2,K)} )

### Linear heat rate (LHGR), 1 year

* Ramp-up: 0 → **40 kW/m** in **0.04 days** (~1 h)
* Hold: **40 kW/m** until **day 364**
* Ramp-down: **40 → 0 kW/m** in **0.04 days**
* Cool-down to **day 365**

### Materials (constant k)

* Fuel conductivity: ( $k_f$ = 3 W/(mK) )
* Clad conductivity: ( $k_c$ = 21 W/(mK) )

---

## Tasks

1. **Copy the case folder.**
   Duplicate `Case1/` into your working directory (keep `exercise01.py` and this `README.md`).

2. **Set up geometry & mesh in Python.**
   Open `exercise01.py` and edit the **Parameters** cell at the top (all lengths in **SI**—use `mm = 1e-3`).
   Create the 1D rod mesh with:

   * `mesh.rod_1d(...)` using the values above
   * Fuel/clad radial cell counts: 30 / 10
   * Axial slices: 1

3. **Define fields and boundary conditions.**

   * Create `TimeFolder(0)`
   * Create Temperature `Field` with internal field **293 K**
   * Apply BCs:

     * `fuelOuter` and `cladInner`: **ResistiveGap** with ( $h_{gap}$=5000 )
     * `cladOuter`: **FixedValue (@573 K)**

4. **Assign material properties.**

   * Use `materials.Constant.preset(name="fuel", k=3.0)`
   * Use `materials.Constant.preset(name="cladding", k=21.0)`
   * Ensure the `name` matches mesh zone names (e.g., `"fuel"`, `"cladding"`).

5. **Configure the OFFBEAT thermal solver.**

   * Thermal only: `SolidConductionThermalSolver()`
   * Heat source: `TimeDependentLhgrHeatSource(...)` with **days** and **W/m** values as specified
   * Flat axial/radial profiles; apply to **fuel** only

6. **Case settings, probe, export, and run.**

   * Create `Case()`, add the solver and the time folder
   * Use **days** as `userTime`; set `endTime=365` and a reasonable starting `deltaT` (e.g., `0.01`)
   * Add a **probe** at **(0, 0, 0.055 m)** named `fuelCenterline` for field `T`
   * With `model.export_to_openfoam()` the API writes the case files
   * Run with `ffn.run(model=model)`

7. **Post-process (optional).**

   * Read the probe data (centerline temperature vs time) via API utilities
   * Plot/save a simple PNG (e.g., `fig_centerline_T.png`)
   * (Optional) Sample a **radial line** at z = 0.055 m if you want a radial profile

---

## How to run

### Option A — VS Code “percent” notebook

1. Open `exercise01.py`
2. Run cell-by-cell (`# %%` separators)

### Option B — CLI

```bash
python exercise01.py
```

---

## What to submit / check

* A **successful run** with stable temperature evolution
* A **centerline temperature plot** over time (PNG)
* (Optional) A **radial temperature profile** at a chosen time

---
