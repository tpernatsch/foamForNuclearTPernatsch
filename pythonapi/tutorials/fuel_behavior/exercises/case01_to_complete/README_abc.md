# Description

This exercise involves three sub-cases, each a variation of the original Case1. These variations are designed to explore different aspects of thermal behavior in a nuclear fuel rod simulation without the complexity of mechanics and nuclear-specific phenomena.

## Sub-case Description

* **Case1a**: Similar to Case1, but with a different constant value for the fuel thermal conductivity of $k_f = 2,\mathrm{W/(m\cdot K)}$
* **Case1b**: Similar to Case1, but with a different gap conductance, equal to $2000,\mathrm{W/(m^2 \cdot K)}$
* **Case1c**: Similar to Case1, but the fuel in this scenario has a central hole with a radius of $0.75,\mathrm{mm}$, altering its thermal profile.

---

# Tasks (for each sub-case)

1. **Create the folder.**
   Copy the completed `case01` folder, rename it (e.g. `case01a`, `case01b`, `case01c`), and clean it using:

   ```
   python clean.py
   ```

2. **Update `case.py`.**
   Modify `case0Xa/case.py` using the specifications above:

   * Case1a: change the fuel thermal conductivity.
   * Case1b: change the gap conductance.
   * Case1c: introduce a fuel central hole of radius $0.75,\mathrm{mm}$ and update the probe and radial profile definitions accordingly.

3. **Run OFFBEAT.**
   From the sub-case folder (`case01a`, `case01b`, or `case01c`), run:

   ```
   python run.py
   ```

   or

   ```
   python3 run.py
   ```

4. **Post-process.**
   Plot the results from the same folder using:

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

* In **Case1c**, the presence of a central hole means the probe and radial line must start from the **fuel inner radius** instead of zero. This must be reflected in `case.py` when defining `fuelCenterline` and `radialProfile`.

* The `run.py` script handles case generation, meshing, and execution automatically. If you need to restart a case, run:

  ```
  python clean.py
  ```

  before re-running.