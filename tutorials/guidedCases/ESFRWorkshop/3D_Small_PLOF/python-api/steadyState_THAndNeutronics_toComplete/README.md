# Steady-State Coupled TH + Neutronics — Python API  [TO COMPLETE]

This is the student version of the coupled TH + neutronics tutorial.
Complete the tasks below in `caseSetup.py`, then run `case.py`.

The completed reference is in `../steadyState_THAndNeutronics/`.

---

## Tasks

1. **`caseSetup.py`** — Add the coupling field transfers.

   The coupling block looks like this:
   ```python
   coupling.add_field_transfer(from_solver, to_solver, "source_field", "target_field")
   ```

   You need to define:
   - **Neutronics → TH** (2 transfers): power density to the fluid mesh
   - **TH → Neutronics** (5 transfers): coolant and fuel temperatures,
     coolant density, passive structure temperature

   Use the field name reference table at:
   https://foam-for-nuclear.gitlab.io/GeN-Foam/user_guide/coupling/

2. **`caseSetup.py`** — Set `settings.endTime = 200`.

3. **Run the case**:
   ```bash
   python case.py
   ```

4. **Compare with stand-alone cases**: Does adding feedback change k_eff?
   Why does the outlet temperature now vary across assemblies?
