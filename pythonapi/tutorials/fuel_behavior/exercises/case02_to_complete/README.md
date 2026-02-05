# Description

In Case 2, we start incorporating nuclear-specific phenomena into our thermal analysis. While maintaining the same settings as in Case1, a key difference here is that thermal conductivity is no longer constant but varies as a function of temperature. Additionally, as for the UO2 fuel conductivity is also significantly influenced by burnup, we activate the burnup module which integrates the volumetric power density (times the density) over time.

---

# Tasks

1. **Create the folder.**
   Either copy `case02/` in your working directory, or duplicate `case01/` and rename it to `case02`.
   Then clean the folder using:

   ```
   python3 clean.py
   ```

2. **Update `case02/case.py`.**
   Modify the case to activate temperature- and burnup-dependent material behaviour:

   * **Use default conductivity correlations.**
     Remove any explicit constant conductivity model so that the built-in temperature-dependent correlations are used for fuel and cladding.

   * **Activate the burnup solver.**
     Enable a burnup model based on power integration (the `fromPower` model), which computes burnup from the volumetric power density over time.

   * **Track burnup at the fuel centerline.**
     Add the burnup field `Bu` to the `fuelCenterline` probe so that burnup is recorded together with temperature.

3. **Run OFFBEAT.**
   From the `case02` directory, run:

   ```
   python run.py
   ```

   or

   ```
   python3 run.py
   ```

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
   * `case/postProcessing/fuelCenterline/0/Bu`
   * `case/postProcessing/radialProfile/##/line_T.csv`

   using:

   ```
   python extractData.py
   ```

   or

   ```
   python3 extractData.py
   ```

6. **Copy the extracted data into the MATLAB grader.**

---

## Notes

* The **geometry, mesh, boundary conditions, and power history** are the same as in Case 01; only the **material models and burnup physics** change.

* The `run.py` script automatically generates the OpenFOAM case, creates the mesh, and runs OFFBEAT.If you need to restart the simulation from scratch, run:

  ```
  python3 clean.py
  ```
  before running again.