"""
case.py — Steady-state coupled thermal-hydraulics + neutronics (3D Small ESFR)
===============================================================================

WHAT THIS CASE DOES
-------------------
Solves the coupled neutronics and thermal-hydraulics of the ESFR core.
The two physics solvers iterate together (Picard iteration) until a
self-consistent solution is reached:

  1. Neutronics computes the power distribution (flux shape × total power)
  2. TH uses this power to heat the coolant → temperatures and density
  3. Neutronics updates cross sections with new temperatures/density → new flux
  4. Repeat until both converge

This coupling captures the most important reactor physics feedbacks:
  Doppler effect         : higher fuel temperature → broader absorption resonances
                           → fewer fast neutrons survive → lower reactivity
  Coolant void effect    : lower Na density → less moderation and absorption
                           → changes reactivity (negative in the ESFR)
  Structural temperature : thermal expansion of passive structures
                           → small reactivity effect

HOW TO RUN
----------
    python case.py

For parallel execution (4 cores):
    decomposePar -allRegions
    mpirun -np 4 GeN-Foam -parallel | tee log.GeN-Foam
"""

import foamForNuclear as ffn
import matplotlib.pyplot as plt

from mesh               import nMesh, thMesh
from boundaryConditions import timeFolder0
from caseSetup          import model, neutronicsSolver, thSolver, solvers, coupling

# ---------------------------------------------------------------------------
XS_SOURCE = "../../steadyState_THNeutronicsAndTM/constant/neutroRegion"
XS_FILES  = [
    'nuclearData', 'XSref',
    'XSaxialExpansion', 'XSradialExpansion',
    'XSrhoCool500kgm3', 'XSTClad1950K', 'XSTFuel1200K',
]

if __name__ == "__main__":

    # -----------------------------------------------------------------------
    # 1. Export → OpenFOAM input files
    # -----------------------------------------------------------------------
    ffn.allclean()
    model.export_to_openfoam()
    print(model)

    # Visualise how the two physics solvers exchange data
    coupling.plot_coupling_graph()
    coupling.plot_solving_graph()

    # -----------------------------------------------------------------------
    # 2. Copy cross-section data
    # -----------------------------------------------------------------------
    for xs_file in XS_FILES:
        ffn.copyFolder(f"{XS_SOURCE}/{xs_file}", f"constant/{nMesh.region}")

    # -----------------------------------------------------------------------
    # 3. Pre-processing
    # -----------------------------------------------------------------------
    ffn.run_preprocessing(model=model)

    model.plot_mesh(region=nMesh,  show_edges=True)
    model.plot_mesh(region=thMesh, show_edges=True)

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

    model.plot_residuals(parameters=['fluxStar0'],   title="Neutronics residuals")
    model.plot_residuals(parameters=['p_rgh', 'h'],  title="TH residuals")

    # k_eff convergence
    res = model.get_keff_from_log()
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(res['time'], res['keff'])
    ax.axhline(1.0, color='r', linestyle='--', label='critical')
    ax.set_xlabel("Outer iteration"); ax.set_ylabel(r"$k_\mathrm{eff}$")
    ax.legend(); fig.tight_layout(); fig.savefig("keff_convergence.png")

    # Temperature fields on the neutronics mesh (averaged over the volume)
    for field in ['TCool', 'TFuel', 'TClad', 'TStructMech']:
        model.plot_slice(
            region=nMesh, time=lastTime,
            fieldName=field, normal='x',
            cmap='RdBu_r', unit='K', show_edges=True
        )

    # Axial power profile along the core centreline
    model.plot_slice(
        region=nMesh, time=lastTime,
        fieldName='powerDensity', normal='z',
        cmap='inferno', unit='W/m³', show_edges=True
    )

    # Coolant temperature distribution
    model.plot_mesh(
        region=thMesh, time=lastTime,
        fieldName='T', cmap='RdBu_r', unit='K'
    )
