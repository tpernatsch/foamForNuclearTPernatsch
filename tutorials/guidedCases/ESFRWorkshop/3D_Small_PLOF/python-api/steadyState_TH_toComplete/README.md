# Stand-alone Steady-State Thermal-Hydraulics — Python API  [TO COMPLETE]

This is the student version of the TH tutorial.  Complete the tasks below
in `boundaryConditions.py` and `caseSetup.py`, then run `case.py`.

The completed reference is in `../steadyState_TH/`.

---

## Tasks

1. **`boundaryConditions.py`** — Set boundary conditions for temperature `T`:
   - `FixedValue(INLET_TEMPERATURE)` at the inlet (`bottom`)
   - `ZeroGradient()` at all other boundaries

2. **`boundaryConditions.py`** — Set boundary conditions for velocity `U`:
   - `FixedValue(U.internalField)` at the inlet (`bottom`)
   - `ZeroGradient()` at the outlet (`top`)
   - `Slip()` at walls and baffles

3. **`caseSetup.py`** — Set `solver = "onePhase"` in the
   `ThermalHydraulicsSolver` definition.

4. **Run the case** once all TODOs are filled:
   ```bash
   python case.py
   ```

5. **Check the result**: with zero power, the coolant temperature should
   remain at 668 K everywhere.  Verify the velocity distribution is
   symmetric and physically sensible.
