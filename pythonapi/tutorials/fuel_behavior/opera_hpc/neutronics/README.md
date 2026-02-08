# **Exercise 1 — Neutronics: Radial Power & Fuel Composition Evolution**

## **Objective**

In this exercise we analyse a **1-D UO₂ fuel rod** under a constant LHGR and study how neutronics and fuel composition evolve over several hundred days of irradiation.

You will:

* Build the case from the Python API (`case.py`)
* Run a diffusion neutronics + burnup calculation
* Track **radial profiles** of burnup & nuclide densities
* Post-process via `post.py` to generate plots automatically

This is a *neutronics-only* exercise:
**no mechanics, gap models, creep, relocation, or thermal feedback**.

---

## **What the script sets up**

### **Geometry & Mesh**

* Cylindrical 1-D fuel rod
* Radius: **4.65 mm**
* Axial thickness: **1 mm** (slice)
* Radial discretization: **1000 cells**
* Single axial slice → focus is on **radial effects only**

Configured in `mesh.rod_1d(...)` in **`case.py`**:


### **Physics**

| Case      | Object                                   |
| ---------- | ---------------------------------------- |
| Neutronics | Diffusion solver                         |
| Burnup     | Lassmann model                           |
| Fuel       | UO₂ enriched to 4.6%                     |
| Power      | Constant LHGR = 60 kW/m for 1096 days    |
| Time       | End time 600 days, output every 100 days |

---

## **Data extraction configured**

From `case.py`:
✔ centerline probe (Bu & T)
✔ radial line profile with ~1000 equally spaced points
✔ volume averages of isotopes
✔ integral HM inventory (U + Pu)


From `post.py` (automatic plots):

* Burnup vs time at centerline
* Normalized radial profiles (Bu, U-235, Pu-239)
* Heavy-metal inventory evolution


---

## **How to run**

### **1️⃣ Build and export the case**

```bash
python3 case.py
```

This generates the OpenFOAM case structure and writes mesh, fields, and function objects.

### **2️⃣ Run**

```bash
python3 run.py
```

### **3️⃣ Post-process**

```bash
python3 post.py
```

This will save:

* `fig_center.png` → Burnup evolution at centerline
* `fig_radial_Bu.png` → Radial burnup + isotope profiles
* `fig_N_evolution.png` → HM isotopic inventory evolution vs burnup

---

## **What you should observe**

* Burnup grows smoothly at the rod centre
* Radial burnup profile becomes **rim-peaked**
* U-235 fraction decreases; Pu-239/240/241 builds up
* Ratio Pu/HM increases with burnup

---

## **Questions to answer**

1. Why does burnup peak near the outer radius?
2. Which Pu isotope grows fastest, and why?
3. Does the U-235 radial profile flatten or sharpen with time?
4. How would a **non-flat** LHGR profile modify results?
5. If we added temperatures + Doppler, what trend would you expect?

---

## **What to submit**

* `fig_center.png` (burnup vs time)
* `fig_radial_Bu.png` (radial normalized profile)
* Short paragraph answering the 5 questions above