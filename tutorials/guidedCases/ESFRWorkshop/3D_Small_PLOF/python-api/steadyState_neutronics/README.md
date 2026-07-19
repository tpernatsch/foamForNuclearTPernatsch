# Stand-alone Steady-State Neutronics — Python API

Performs an eigenvalue neutronics calculation for the 3D Small ESFR core using
the GeN-Foam Python API.  This is the Python-API equivalent of the
`../../../steadyState_neutronics` OpenFOAM case.

---

## Tasks

1. **Set neutron flux boundary conditions** in `boundaryConditions.py`.
   Apply a black (zero-flux) boundary condition on all outer surfaces.

2. **Configure the neutronics solver** in `caseSetup.py`:
   - Choose the solver type (`"diffusionNeutronics"` or `"SP3Neutronics"`)
   - Set the total reactor power and the initial k_eff guess
   - Set the number of neutron transport inner iterations

3. **Run the case**:
   ```bash
   python case.py
   ```
   Or in parallel (4 cores):
   ```bash
   decomposePar -allRegions
   mpirun -np 4 GeN-Foam -parallel | tee log.GeN-Foam
   ```

4. **Examine the results**: k_eff convergence plot, neutron flux shape,
   power density distribution.

5. **(Optional)** Switch to `"SP3Neutronics"` and compare k_eff values.

---

## File Structure

| File | Purpose |
|------|---------|
| `mesh.py` | Load or generate the neutronics mesh |
| `boundaryConditions.py` | Neutron flux initial and boundary conditions |
| `caseSetup.py` | Neutronics solver setup and simulation settings |
| `case.py` | Main script: export, run, post-process |
| `../meshes/polyMeshNeutro/` | Pre-generated neutronics mesh |
