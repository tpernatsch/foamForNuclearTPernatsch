## Description

In the next few exercises, we study a short rodlet (10 pellets) over a simplified **one-year** irradiation. We start with the basics and perform **thermal analysis only**: no relocation, densification, or creep-down. Think of the fuel as a hot cylinder that releases heat to the Zircaloy cladding across a **fixed-width gap** with **fixed thermal contact conductance**.

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

* Initial temperature: $T_0 = 293,\mathrm{K}$
* Cladding outer temperature: $T_{c0} = 573,\mathrm{K}$
* Fixed gap conductance: $h_{gap} = 5000,\mathrm{W/(m^2,K)}$

### Linear heat rate (LHGR), 1 year

* Ramp-up: $0 \rightarrow 40,\mathrm{kW/m}$ in $0.04$ days (~1 h)
* Hold: $40,\mathrm{kW/m}$ until day $364$
* Ramp-down: $40 \rightarrow 0,\mathrm{kW/m}$ in $0.04$ days
* Cool-down to day $365$

### Materials (constant $k$)

* Fuel conductivity: $k_f = 3,\mathrm{W/(m,K)}$
* Fuel grain radius: $rGrain = 28e-6$
* Clad conductivity: $k_c = 21,\mathrm{W/(m,K)}$

---

## Tasks

1. **Copy `case01/` in your working directory.**

2. **Update `case01/case.py`.**
   Use the data above to complete the simulation setup. In particular:

   * Complete the **geometry and mesh settings** (rod dimensions, number of radial cells).
   * Complete the **initial and boundary conditions** for the temperature field.
   * Complete the **material properties** and the **LHGR time history** (ramp-up / hold / ramp-down / cool-down).

3. **Run OFFBEAT.**
   From the `case01` directory, run:

   ```
   python run.py
   ```

   or

   ```
   python3 run.py
   ```

   This will generate the OpenFOAM case, create the mesh, and run OFFBEAT.

4. **Post-process.**
   From the same directory, plot the results with:

   ```
   python post.py
   ```

   or

   ```
   python3 post.py
   ```

5. **Extract data for the MATLAB grader.**
   Extract:

   * `case/postProcessing/fuelCenterline/0/T`
   * `case/postProcessing/radialProfile/##/line_T.xy`

   You can do this manually, or run from the case root:

   ```
   python extract.py
   ```

   or

   ```
   python3 extract.py
   ```

6. **Copy the extracted data into the MATLAB grader.**

---

## Notes

* The `run.py` script handles case generation, meshing, and execution automatically. If you need to restart a case, run:

   ```
   python clean.py
   ```

   or

   ```
   python3 clean.py
   ```

  before re-running.