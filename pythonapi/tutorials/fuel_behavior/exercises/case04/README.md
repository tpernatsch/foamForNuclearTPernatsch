# Description

Case 4 is identical to Case 3 (`case03_a`), but it introduces two crucial nuclear fuel phenomena: densification and relocation. Densification refers to the closing of porosity in nuclear fuel leading to pellet shrinkage and gap opening. Relocation is the outward movement of fuel fragments post-cracking. This phenomenon is considered the major contributor to gap closing.

---

# Tasks

1. **Create the folder.**
   Either copy `case04/` in your working directory, or duplicate `case03_a/` and rename it to `case04`.
   Then clean the folder using:

   ```
   python3 clean.py
   ```

2. **Update `case04/case.py`.**
   Introduce densification and relocation in the fuel material model:

   * **Introduce a densification model.**
     In the fuel material settings, enable `densificationModel` and select the **empirical** model.

   * **Complete the densification model settings.**
     Set:

     * `densificationDensityChange` = 1 (percent)
     * `densificationTimeConstant` = 1000 $\mathrm{MWd/tU}$`

   * **Introduce a relocation model.**
     In the fuel material settings, enable `relocationModel` and select the **UO2Frapcon** model.

   * **Complete the relocation model settings.**
     Provide:

     * `GapCold` (in $\mathrm{m}$)
     * `DiamCold` (in $\mathrm{m}$)
     * the name of the **fuel outer patch** (outerPatch)

3. **Run OFFBEAT.**
   From the `case04` directory, run:

   ```
   python run.py
   ```

   or

   ```
   python3 run.py
   ```

4. **Post-process.**
   Plot the results with:

   ```
   python post.py
   ```

5. **Extract data for the MATLAB grader.**
   Extract:

   * `case/postProcessing/fuelCenterline/0/T`
   * `case/postProcessing/fuelOuter/0/gapWidth`

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

* Case 04 is identical to **Case 03** (thermo-mechanics + gap model), but adds **densification** and **relocation** in the fuel behaviour.

* The `run.py` script handles case generation, meshing, and execution.
  To restart from scratch, run:

  ```
  python clean.py
  ```

  before running again.
