"""
caseSetup.py — Neutronics solver configuration  [TO COMPLETE]
=============================================================

TASK 2 — Choose the neutronics solver
--------------------------------------
GeN-Foam offers two neutron transport approximations:
  "diffusionNeutronics" — multi-group diffusion (faster, standard choice)
  "SP3Neutronics"       — simplified P3 transport (more accurate, slower)

Fill in the solver name in the NeutronicsSolver definition below.

TASK 3 — Set the reactor power and initial k_eff guess
-------------------------------------------------------
The solver needs to know the total reactor power so it can normalise the
flux to physical units.  It also needs a starting estimate for k_eff
(the effective multiplication factor) to begin the eigenvalue iteration.

  REACTOR_POWER  : total thermal power of the reactor [W]  → use 8×10⁸ W
  INITIAL_KEFF   : starting guess for k_eff [-]            → use 0.9388902

TASK 4 — Set the number of neutron transport iterations
--------------------------------------------------------
More inner iterations per outer step gives a more accurate flux at each
step but takes longer.  Fill in the number of iterations (use 50).
"""

import foamForNuclear as ffn

from mesh               import nMesh
from boundaryConditions import timeFolder0

# ---------------------------------------------------------------------------
# TODO: Set the reactor power [W] and the initial k_eff guess [-]
# ---------------------------------------------------------------------------
REACTOR_POWER = ???    # total thermal power [W]
INITIAL_KEFF  = ???    # initial eigenvalue guess [-]

# ---------------------------------------------------------------------------
# Neutronics solver
# ---------------------------------------------------------------------------
neutronicsSolver = ffn.NeutronicsSolver(
    # TODO: choose the solver — "diffusionNeutronics" or "SP3Neutronics"
    solver  = "???",
    mesh    = nMesh,
    region  = nMesh.region,
    power   = REACTOR_POWER,
    keff    = INITIAL_KEFF,
    isMeshDeformation = False,
)

# TODO: Set the maximum number of neutron transport inner iterations per step
neutronicsSolver.neutronTransportOptions.maxNeutronIterations = ???

# ---------------------------------------------------------------------------
# Solver list, coupling (empty for stand-alone neutronics), and model
# ---------------------------------------------------------------------------
solvers  = ffn.Solvers([neutronicsSolver])

for solver in solvers:
    solver.decomposeParDict.numberOfSubdomains = 4
    solver.decomposeParDict.method             = 'simple'
    solver.decomposeParDict.simpleCoeffs['n']  = ffn.Vector(1, 1, 4)

coupling = ffn.Coupling(solvers=solvers)

model    = ffn.Model(
    solvers     = solvers,
    coupling    = coupling,
    timeFolders = [timeFolder0]
)

settings = model.settings
settings.application       = 'GeN-Foam'
settings.endTime           = 1
settings.deltaT            = 0.01
settings.writeControl      = 'runTime'
settings.writeInterval     = 1
settings.writePrecision    = 7
settings.runTimeModifiable = False
settings.adjustTimeStep    = False
