"""
Run ESFR transient
"""
# Import model from steady-state
from Allrun import *

# Modify parameters
neutronicsSolver.eigenvalueNeutronics = False

fluxSolver = ffn.fvSolutionSolver(
    solver="GAMG",
    tolerance=1e-6,
    relTol=1e-2,
    smoother="DICGaussSeidel",
    nPostSweeps=1,
    nPreSweeps=1,
    nFinestSweeps=1,
    cacheAgglomeration=True,
    nCellsInCoarsestLevel=20,
    agglomerator="faceAreaPair",
    mergeLevels=1,
    processorAgglomerator="masterCoarsest"
)
# precSolver = ffn.fvSolutionSolver(
#     solver="PBiCG",
#     preconditioner="DILU",
#     tolerance=1e-6,
#     relTol=1e-3,
# )

neutronicsSolver.fvSolution.solvers['"flux.*"'] = fluxSolver
neutronicsSolver.fvSolution.solvers['"adjoint_flux.*"'] = fluxSolver

model.settings.deltaT = 0.001
model.settings.endTime = 500
model.settings.writeControl = "runTime"
model.settings.writeInterval = 10
model.settings.adjustTimeStep = True
model.settings.runTimeModifiable = True
model.settings.maxDeltaT = 0.1
model.settings.maxCo = 10
model.settings.maxPowerVariation = 0.025

# Re-export to update the case
model.export_to_openfoam()
print(model)

# Duplicate and overwrite files
for filename in [
    'nuclearData', 'XSaxialExpansion', 'XSradialExpansion', 'XSref',
    'XSrhoCool500kgm3', 'XSTClad1950K', 'XSTFuel1200K'
]:
    ffn.copyFolder(f"../steadyState_THNeutronicsAndTM/constant/neutroRegion/{filename}", f"constant/{nMesh.region}")


# Run the case
ffn.run(model=model)
