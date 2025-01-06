# Imports
from pyfmi import load_fmu
from pyfmi.master import Master
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the FMUs
controller = load_fmu("pidControllers.fmu")
controller.initialize()

models = [controller]

# Connections from output to input
connections = []

# controller.set('measFrequency', np.float64(1.0))

# system.set("port",6006)
# system.set("outputPath", "TempControl")

# Create the Master simulator
master_simulator = Master(models, connections)

controller.set("measFrequency", 50)

def getValue(model, name: str):
    valueRef = model.get_variable_valueref(name)
    return(model.get_real([valueRef])[0])


# Run the simulation
t, dt = 0, 0.005
while (t < 5):
    print(t)
    for model in master_simulator.models:
        status = model.do_step(t, dt, False)
        print(
            getValue(model, "measFrequency"), 
            " -> ", 
            getValue(model, "pidFrequencyResponse")
        )

        
    t += dt

    if (t > 2.5):
        controller.set("measFrequency", 25)