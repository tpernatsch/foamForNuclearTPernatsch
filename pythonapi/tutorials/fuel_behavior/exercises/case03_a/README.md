# Description

In Case 3, we increase the complexity of our simulation and move to thermo-mechanics by activating the mechanical solver. The code now solves also for the radial displacement, following the 1.5D approximation in small strain, and the gap size changes over time due to thermal expansion. For this reason, we also include a gap model and related BC to include the effect of the changing gap size on the temperature distribution. Then, in Case 3a, we apply a single but important modification: starting from Case 3, we apply the isotropic cracking model so that the presence of cracks in the fuel pellet which affects the fuel stress state can be taken into account.

---

# Tasks — Case 03

1. **Create the folder.**
   Either copy `case03/` in your working directory, or duplicate `case02/` and rename it to `case03`.
   Then clean the folder using:

   ```
   python3 clean.py
   ```

2. **Update `case03/case.py`.**
   Modify the case to activate thermo-mechanical coupling and a dynamic gap:

   * **Use a thermo-mechanical gap boundary condition.**
     Replace the temperature boundary conditions on `fuelOuter` and `cladInner` with a `fuelRodGap` boundary condition so that heat transfer depends on the evolving gap size.

   * **Activate mechanical analysis.**
     Enable a small-strain mechanical solver (`smallStrain`) and set the rheology to use the **1.5D (modified plane strain)** approximation.

   * **Apply system pressure.**
     Apply a system pressure of **10 MPa** on the **outer cladding surface** in the displacement boundary conditions.

   * **Activate a gap gas model.**
     Enable the `Frapcon` gap-gas model.

   * **Initialize the gap gas.**
     Set the gap to contain only **He** with an initial pressure of **20 bar**.

   * **Monitor the evolving gap width.**
     Add a probe on the `fuelOuter` surface to record `gapWidth` during the simulation.

3. **Run OFFBEAT.**
   From the `case03` directory, run:

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

   or

   ```
   python3 post.py
   ```

5. **Extract data for the MATLAB grader.**
   Extract:

   * `case/postProcessing/fuelCenterline/0/T`
   * `case/postProcessing/fuelOuter/0/gapWidth`
   * `case/postProcessing/radialProfile/##/line_T.csv`
   * `case/postProcessing/radialProfile/##/line_sigma.csv`

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

# Tasks — Case 03a (Fuel cracking)

1. **Create the folder.**
   Duplicate `case03` and rename it to `case03a`.
   Then clean the new folder using:

   ```
   python3 clean.py
   ```

2. **Update `case03a/case.py`.**

   * **Activate isotropic cracking.**
     Enable the `isotropicCracking` model for the fuel material and set
     `nCracksMax = 12`.

3. **Run OFFBEAT.**
   From the `case03a` directory, run:

   ```
   python run.py
   ```

   or

   ```
   python3 run.py
   ```

4. **Post-process.**

   ```
   python post.py
   ```

5. **Extract data for the MATLAB grader.**
   Extract:

   * `case/postProcessing/radialProfile/168.04/line_sigma.csv`

   using:

   ```
   python extractData.py
   ```

6. **Copy only the **radial** and **hoop stress** at **168.04 days** into the MATLAB grader.**

---

## Notes

* Geometry, mesh, power history, and thermal models are inherited from **Case 02**.
* The `run.py` script performs case generation, meshing, and execution.
* To restart any case from scratch, run:

  ```
  python clean.py
  ```
  
  before running again.