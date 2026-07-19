############################################################
# %% Imports
############################################################
import foamForNuclear as ffn
import matplotlib.pyplot as plt
import numpy as np

solver = ffn.solvers.offbeat.Offbeat()
solver.import_from_openfoam()

solver.thermalSolver = ffn.offbeat_lib.misc.SolidConductionThermalSolver()

solver.export_to_openfoam()






# model = ffn.Case()
# model.import_from_openfoam()

# print(model)

# fuelCenterline = ffn.functions.Probes(
#     name="fuelCenterline",
#     fields=['Bu', 'T'],
#     probeLocations=[ffn.Vector(0.0, 0.0, fuel_length/2)],  # NOTE: in meters;
#     enabled=True
# )

# radialProfile = ffn.functions.Graph(
#     name="radialProfile",
#     start=ffn.Vector(fuel_ri, 0.0, fuel_length/2),
#     end=ffn.Vector(fuel_ro, 0.0,  fuel_length/2),
#     fields=["D", "Bu", "T", "N_U", "N_U235", "N_U238", "N_Pu238", "N_Pu239", "N_Pu240", "N_Pu241", "N_Pu242"],
#     nPoints=fuel_nr,
#     writeControl = "timeStep"
# )

# volAverage = ffn.functions.VolFieldValue(
#     name="volAverage",
#     fields=["D", "Bu", "T", "N_U", "N_U235", "N_U238", "N_Pu238", "N_Pu239", "N_Pu240", "N_Pu241", "N_Pu242"],
#     operation="volAverage",
#     regionName="fuel"
# )

# sumFields = ffn.functions.VolFieldValue(
#     name="sum",
#     fields=["N_U", "N_U235", "N_U238", "N_Pu238", "N_Pu239", "N_Pu240", "N_Pu241", "N_Pu242"],
#     operation="sum",
#     regionName="fuel"
# )