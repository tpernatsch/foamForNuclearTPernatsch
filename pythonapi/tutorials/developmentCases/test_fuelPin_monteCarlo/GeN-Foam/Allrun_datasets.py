#==============================================================================*
# Imports

from Allrun import *

import tabulate


#==============================================================================*
# Small test loop over solvers

keff = []

for solver in ["diffusionNeutronics", "SP3Neutronics"]:
    neutronicsSolver.solver = solver

    temp = [solver]

    for nuclearDataSet in [serpentState, openmcState]:

        # Empty nuclear data
        neutronicsSolver.nuclearData.states.pop(0)
        neutronicsSolver.nuclearData.add_state(nuclearDataSet)

        # Clean and export
        ffn.cleanCase0()
        model.export_to_openfoam()

        # Run
        ffn.run(model, is_preprocessing=True)

        # Extract keff
        temp.append(model.keff())

    keff.append(temp)


print("Disclaimer: MGXS might not be well converged")
print(tabulate.tabulate(
    keff,
    headers=['Solver', 'Serpent', 'OpenMC']
))


#==============================================================================*
