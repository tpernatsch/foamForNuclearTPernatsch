"""
Script to test the coupling between 2 FMUs from Modelica and Simulink. 
Works perfectly.

Author: Thomas Guilbaud, EPFL/Transmutex, 09/05/2023
"""

# Imports
from pyfmi import load_fmu
from pyfmi.master import Master
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the FMUs
system = load_fmu("lumpedSystem.fmu")
controller = load_fmu("pidControllers.fmu")

models = [system, controller]

# Connections from output to input
connections = [ 
    (system, "y", controller, "measFrequency"),
    (controller, "pidFrequencyResponse", system, "u")
]

# controller.set('measFrequency', np.float64(1.0))

# system.set("port",6006)
# system.set("outputPath", "TempControl")

# Create the Master simulator
master_simulator = Master(models, connections)
opts = master_simulator.simulate_options()
opts["step_size"] = 0.005
opts["result_handling"] = "file"

# Run the simulation
res = master_simulator.simulate(
    start_time = 0.0,
    final_time = 100.0,
    options = opts
)

# Extract the results
results = pd.DataFrame()
results["time"] = res[controller]["time"]
# results["cmd"] = res[system]["step.y"]
results["y"] = res[system]["y"]
results["measFrequency"] = res[controller]["measFrequency"]
results["pidFrequencyResponse"] = res[controller]["pidFrequencyResponse"]
results.to_csv("results.csv",index=False)

# Plot the results
plt.plot(results["time"], [f/50 for f in results["measFrequency"]], label="measFrequency")
# plt.plot(results["time"], results["cmd"], label="Cmd")
plt.plot(results["time"], [f/50 for f in results["y"]], label="y")
plt.plot(results["time"], results["pidFrequencyResponse"], label="pidFrequencyResponse")
plt.grid()
plt.legend()
plt.show()