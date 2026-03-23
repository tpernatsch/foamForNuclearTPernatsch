# Stand-alone Steady-State Thermal-Hydraulics — Python API

Simulates the sodium coolant flow through the ESFR core without any nuclear
power source.  This is the Python-API equivalent of the
`../../../steadyState_TH` OpenFOAM case.

---

## Tasks

1. **Review `boundaryConditions.py`**: understand the temperature and
   velocity boundary conditions at the inlet, outlet, walls and baffles.

2. **Review `caseSetup.py`**: understand how the porous-medium model,
   fuel pin model, drag and heat transfer correlations are set up.

3. **Run the case**:
   ```bash
   python case.py
   ```

4. **Examine the results**: coolant velocity distribution, temperature
   (should be uniform at 668 K — no power source).

---

## File Structure

| File | Purpose |
|------|---------|
| `mesh.py` | Load or generate the fluid-region mesh |
| `boundaryConditions.py` | Temperature, velocity, pressure IC and BC |
| `caseSetup.py` | TH solver, sodium properties, porous-medium models |
| `case.py` | Main script: export, run, post-process |
| `../meshes/polyMeshFluid/` | Pre-generated fluid mesh (with baffles) |
