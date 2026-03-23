"""
caseSetup.py — Neutronics solver configuration
===============================================

This file configures the neutronics solver and the simulation settings.

NEUTRON DIFFUSION SOLVER
-------------------------
GeN-Foam solves the multi-group neutron diffusion equation to find:
  - The effective multiplication factor k_eff (criticality indicator)
  - The neutron flux spatial distribution (and therefore the power shape)

We run in EIGENVALUE mode: the solver iterates over the flux shape and k_eff
simultaneously until both converge.  The result tells us how far the reactor
is from criticality (k_eff = 1 means exactly critical) and where the power
is being generated inside the core.

k_eff initial guess: The solver needs a starting estimate.  The value here is
close to the expected answer, which speeds up convergence.  A rough guess
(e.g. 1.0) also works but requires more iterations.
"""

import foamForNuclear as ffn

from mesh               import nMesh
from boundaryConditions import timeFolder0

# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------
REACTOR_POWER  = 8e8        # Total reactor thermal power [W]
INITIAL_KEFF   = 1  # Initial guess for k_eff [-]  (converges ~10–50 iters)

# ---------------------------------------------------------------------------
# Neutronics solver
# ---------------------------------------------------------------------------
neutronicsSolver = ffn.NeutronicsSolver(
    solver  = "diffusionNeutronics",   # multi-group diffusion (fastest option)
                                       # alternative: "SP3Neutronics" (more accurate)
    mesh    = nMesh,
    region  = nMesh.region,
    power   = REACTOR_POWER,
    keff    = INITIAL_KEFF,
    isMeshDeformation = False,         # no structural deformation in this stand-alone case
)

# Number of inner neutron transport iterations per outer time step.
# More iterations → more accurate flux per step, but slower.
# 50 is a good balance for eigenvalue calculations on this core size.
neutronicsSolver.neutronTransportOptions.maxNeutronIterations = 50

# ---------------------------------------------------------------------------
# Solver list  (only one solver for this stand-alone case)
# ---------------------------------------------------------------------------
solvers = ffn.Solvers([neutronicsSolver])

# Parallel decomposition: split the domain across 4 CPU cores axially
for solver in solvers:
    solver.decomposeParDict.numberOfSubdomains = 1
    solver.decomposeParDict.method             = 'simple'
    solver.decomposeParDict.simpleCoeffs['n']  = ffn.Vector(1, 1, 4)

# ---------------------------------------------------------------------------
# Coupling  (not needed for stand-alone neutronics — no other physics)
# ---------------------------------------------------------------------------
coupling = ffn.Coupling(solvers=solvers)

# ---------------------------------------------------------------------------
# Simulation model
# ---------------------------------------------------------------------------
model = ffn.Model(
    solvers     = solvers,
    coupling    = coupling,
    timeFolders = [timeFolder0]
)

# ---------------------------------------------------------------------------
# Time and output settings
# ---------------------------------------------------------------------------
# In eigenvalue mode the "time" axis represents successive global iterations,
# not physical time.  endTime = 50 means up to 50 outer iterations.
settings = model.settings
settings.application        = 'GeN-Foam'
settings.endTime            = 1       # simulation "time" in eigenvalue mode (not physical time)
settings.deltaT             = 0.01   # each step = one outer neutronics iteration
settings.writeControl       = 'runTime'
settings.writeInterval      = 1      # write every step
settings.writePrecision     = 7
settings.runTimeModifiable  = False
settings.adjustTimeStep     = False
