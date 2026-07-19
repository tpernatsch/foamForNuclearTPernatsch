# Stand-alone Steady-State Neutronics — Python API  [TO COMPLETE]

This is the student version of the neutronics tutorial.  Complete the tasks
below in `boundaryConditions.py` and `caseSetup.py`, then run `case.py`.

The completed reference is in `../steadyState_neutronics/`.

---

## Tasks

1. **`boundaryConditions.py`** — Set the neutron flux boundary conditions.
   Apply `bc.FixedValue(0)` on all three outer surfaces (`wall`, `top`, `bottom`).

2. **`caseSetup.py`** — Configure the solver:
   - Choose `solver = "diffusionNeutronics"`
   - Set `REACTOR_POWER = 8e8` W and `INITIAL_KEFF = 0.9388902`
   - Set `maxNeutronIterations = 50`

3. **Run the case** once all TODOs are filled:
   ```bash
   python case.py
   ```

4. **Compare k_eff** with the result from the OpenFOAM case in
   `../../../steadyState_neutronics/`.

5. **(Optional)** Change `solver` to `"SP3Neutronics"` and re-run.
   Does k_eff change?  Why?
