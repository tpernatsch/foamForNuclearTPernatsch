"""
case.py — Stand-alone steady-state neutronics of the 3D Small ESFR
===================================================================

WHAT THIS CASE DOES
-------------------
Solves the multi-group neutron diffusion equations for the ESFR core in
eigenvalue mode.  The result is:
  - k_eff: the effective multiplication factor
            (k_eff < 1 → subcritical, k_eff = 1 → critical, k_eff > 1 → supercritical)
  - Neutron flux shape: the spatial distribution of neutrons in the core,
            which is proportional to the local power density

Temperatures and coolant density are held fixed at their hot-zero-power
values (no thermal-hydraulics or thermomechanics in this step).

HOW TO RUN
----------
    python case.py

For parallel execution (4 cores):
    decomposePar -allRegions
    mpirun -np 4 GeN-Foam -parallel | tee log.GeN-Foam
"""

import foamForNuclear as ffn
import matplotlib.pyplot as plt

# Import all objects defined in the other modules of this case
from mesh               import nMesh
from boundaryConditions import timeFolder0
from caseSetup          import model, neutronicsSolver, solvers, coupling

# ---------------------------------------------------------------------------
# Nuclear data: cross-section libraries
# ---------------------------------------------------------------------------
# Cross-section tables describe how neutrons interact with each material
# as a function of temperature and density.  These files were pre-computed
# with a lattice physics code (Serpent) and are shared with the OpenFOAM cases.
XS_SOURCE = "../../steadyState_THNeutronicsAndTM/constant/neutroRegion"
XS_FILES  = [
    'nuclearData',       # group-wise cross sections for all materials
    'XSref',             # reference-state cross sections
    'XSaxialExpansion',  # XS sensitivity to axial core expansion
    'XSradialExpansion', # XS sensitivity to radial core expansion
    'XSrhoCool500kgm3',  # XS sensitivity to coolant density
    'XSTClad1950K',      # XS sensitivity to cladding temperature
    'XSTFuel1200K',      # XS sensitivity to fuel temperature
]

if __name__ == "__main__":

    # -----------------------------------------------------------------------
    # 1. Export Python objects → OpenFOAM input files
    # -----------------------------------------------------------------------
    ffn.allclean()
    model.export_to_openfoam()
    print(model)

    # Show which solvers talk to which (useful for multi-physics cases)
    coupling.plot_coupling_graph()
    coupling.plot_solving_graph()

    # -----------------------------------------------------------------------
    # 2. Copy cross-section data into the case folder
    # -----------------------------------------------------------------------
    for xs_file in XS_FILES:
        ffn.copyFolder(f"{XS_SOURCE}/{xs_file}", f"constant/{nMesh.region}")

    # -----------------------------------------------------------------------
    # 3. Pre-processing (mesh check, renumbering, etc.)
    # -----------------------------------------------------------------------
    ffn.run_preprocessing(model=model)

    # Visualise the mesh (opens an interactive 3D window)
    model.plot_mesh(region=nMesh, show_edges=True)

    # -----------------------------------------------------------------------
    # 4. Run the simulation
    # -----------------------------------------------------------------------
    ffn.run(model=model)

    # -----------------------------------------------------------------------
    # 5. Reconstruct parallel results (if run with mpirun)
    # -----------------------------------------------------------------------
    if model.is_parallel:
        ffn.run_reconstruction(model, isLatestTime=True)

    # -----------------------------------------------------------------------
    # 6. Post-processing
    # -----------------------------------------------------------------------
    lastTime = model.get_time_steps()[-1]
    print(f"\nFinal k_eff = {model.keff(time=lastTime):.6f}")

    # Convergence history of k_eff over the outer iterations
    res = model.get_keff_from_log()
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(res['time'], res['keff'])
    ax.axhline(1.0, color='r', linestyle='--', label='critical (k=1)')
    ax.set_xlabel("Outer iteration")
    ax.set_ylabel(r"$k_\mathrm{eff}$")
    ax.set_title("Eigenvalue convergence")
    ax.legend()
    fig.tight_layout()
    fig.savefig("keff_convergence.png")
    print("Saved: keff_convergence.png")

    # Flux distribution along the core centreline
    model.plot_residuals(parameters=['fluxStar0'], title="Neutronics residuals")

    model.plot_slice(
        region=nMesh, time=lastTime,
        fieldName='flux0', normal='z',
        cmap='Blues', unit='neutron/m²/s', show_edges=True
    )
    model.plot_slice(
        region=nMesh, time=lastTime,
        fieldName='powerDensity', normal='z',
        cmap='inferno', unit='W/m³', show_edges=True
    )
