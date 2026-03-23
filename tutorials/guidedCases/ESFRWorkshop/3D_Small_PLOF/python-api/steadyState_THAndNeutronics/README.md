# Steady-State Coupled TH + Neutronics — Python API

Fully coupled steady-state simulation: neutronics computes the power shape,
thermal-hydraulics computes the temperatures, and the two exchange feedback
at every iteration.  This is the Python-API equivalent of the
`../../../steadyState_THAndNeutronics` OpenFOAM case.

---

## Tasks

1. **Review `caseSetup.py`**: study the `coupling.add_field_transfer` calls
   and understand which field goes from which solver to which.

2. **Run the case**:
   ```bash
   python case.py
   ```
   In parallel (4 cores):
   ```bash
   decomposePar -allRegions
   mpirun -np 4 GeN-Foam -parallel | tee log.GeN-Foam
   ```

3. **Compare with the stand-alone results**:
   - How does k_eff differ from the pure neutronics case?
   - How does the coolant temperature profile differ from the pure TH case?
   - What is driving those differences?

4. **(Optional) Tight coupling**: edit `caseSetup.py` to wrap the two
   solvers in a Picard loop (see the OpenFOAM README for the syntax).

---

## File Structure

| File | Purpose |
|------|---------|
| `mesh.py` | Load neutronics + fluid meshes |
| `boundaryConditions.py` | All initial and boundary conditions |
| `caseSetup.py` | Both solvers + coupling field transfers |
| `case.py` | Main script: export, run, post-process |
