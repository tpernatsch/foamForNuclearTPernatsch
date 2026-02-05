# Description

In Case 5, we increase further the complexity of our simulation by adding the phenomenon of cladding creep, moving beyond the purely elastic behavior observed in the earlier cases. Due to creep, the Zircaloy cladding accumulates plastic (or permanent) strain when exposed to moderately high temperatures and stresses over extended periods of time. The impact of neutron fast flux, which is known to enhance the creep phenomenon, is now a crucial factor to consider.

Due to differential pressure exerted on the cladding's inner and outer surfaces - typically higher on the outside, especially at lower fission gas release fractions - cladding creep results in the material "creeping" down towards the fuel pellets. This movement, known as cladding creep-down, is a key contributor to the closure of the gap between the fuel and the cladding, thus influencing the overall behavior of the fuel rod under operational conditions.

---

# Tasks

1. **Create the folder.**
   Either copy `case05/` in your working directory, or duplicate `case04/` and rename it to `case05`.
   Then clean the folder using:

   ```
   python3 clean.py
   ```

2. **Update `case05/case.py`.**
   Extend the thermo-mechanical model by activating cladding creep and fast-flux–driven creep enhancement:

   * **Change the cladding rheology model.**
     Set the cladding `constitutiveLaw` to **`misesPlasticCreep`**.

   * **Configure the creep law.**
     Add a `creep` and `yieldStress` model options for the cladding:

     * Select **`ZircaloyLimback`** as the `creepModel`.
     * Set a **relaxation factor of 0.9** for the creep.
     * Define a **constant yield stress of $400,\mathrm{MPa}$**.

   * **Activate a fast-flux model.**
     Enable a `fastFlux` model of type **`timeDependentAxialProfile`**.

   * **Define the fast-flux history.**
     Set the fast flux proportional to power:

     * Ramp from $0$ to $1\times10^{14},\mathrm{n/cm^2/s}$ in $0.04$ days
     * Hold constant until day $364$
     * Ramp down to $0$ in $0.04$ days
     * Remain at $0$ until day $365$
       Use a **flat axial profile** and apply the model to both **fuel and cladding**.

3. **Run OFFBEAT.**
   From the `case05` directory, run:

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

* Case 05 builds on **Case 04** (thermo-mechanics + gap + densification + relocation) and adds **cladding creep** driven by **stress, temperature, and fast neutron flux**.

* The `run.py` script performs case generation, meshing, and execution.
  To restart from scratch, run:

  ```
  python clean.py
  ```

  before running again.