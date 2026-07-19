"""
case.py — Stand-alone steady-state thermal-hydraulics of the 3D Small ESFR
===========================================================================

WHAT THIS CASE DOES
-------------------
Solves the sodium coolant flow and heat transfer through the ESFR core at
nominal operating conditions, but without any nuclear power source.

Because there is no neutronics solver in this step, the power density is
zero throughout the core.  The coolant flows in at 668 K and exits at the
same temperature.  This is useful to:
  - Verify that the flow distribution is correct (all assemblies receive flow)
  - Check that pressure drops and velocities are physically reasonable
  - Establish the baseline for the coupled neutronics step

HOW TO RUN
----------
    python case.py

For parallel execution (4 cores):
    decomposePar -allRegions
    mpirun -np 4 GeN-Foam -parallel | tee log.GeN-Foam
"""

import shutil
import foamForNuclear as ffn
import matplotlib.pyplot as plt

from mesh               import thMesh
from boundaryConditions import timeFolder0
from caseSetup          import model, thSolver, solvers, coupling

# ---------------------------------------------------------------------------
# Pre-computed nuclear power map
# ---------------------------------------------------------------------------
# The stand-alone TH case needs a power distribution to actually heat the
# coolant.  We reuse the pre-computed field from the reference OpenFOAM case
# (originally produced by a coupled neutronics+TH run).

if __name__ == "__main__":

    # -----------------------------------------------------------------------
    # 1. Export Python objects → OpenFOAM input files
    # -----------------------------------------------------------------------
    ffn.allclean()
    model.export_to_openfoam()
    print(model)

    # (no coupling in stand-alone TH — coupling graph is skipped)

    # -----------------------------------------------------------------------
    # 2. Pre-processing
    # -----------------------------------------------------------------------
    ffn.run_preprocessing(model=model)

    model.plot_mesh(region=thMesh, show_edges=True)
    model.plot_boundary(region=thMesh, boundaryName='baffle0')

    # -----------------------------------------------------------------------
    # 3. Run the simulation
    # -----------------------------------------------------------------------
    ffn.run(model=model)

    # -----------------------------------------------------------------------
    # 4. Reconstruct parallel results (if run with mpirun)
    # -----------------------------------------------------------------------
    if model.is_parallel:
        ffn.run_reconstruction(model, isLatestTime=True)

    # -----------------------------------------------------------------------
    # 5. Post-processing
    # -----------------------------------------------------------------------
    lastTime = model.get_time_steps()[-1]

    model.plot_residuals(parameters=['p_rgh', 'h'], title="TH residuals")

    # Temperature distribution through the core
    model.plot_mesh(
        region=thMesh, time=lastTime,
        fieldName='T', cmap='RdBu_r', unit='K'
    )

    # Coolant velocity distribution
    model.plot_slice(
        region=thMesh, time=lastTime,
        fieldName='U', normal='z', cmap='viridis',
        unit='m/s', show_edges=True
    )
