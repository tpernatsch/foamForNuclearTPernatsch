import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# File path
file_path = 'postProcessing/probes/fluid/0/T'

# Load file, skipping comments
with open(file_path, 'r') as f:
    lines = [line.strip() for line in f if not line.strip().startswith('#')]

# Convert data to NumPy array
data = np.array([list(map(float, line.split())) for line in lines])

# Find the row where time == 2.0
time_index = np.where(np.isclose(data[:, 0], 8))[0]

if time_index.size == 0:
    raise ValueError("Time not found in the data.")
else:
    time_row = data[time_index[0], 1:]  # skip time column

# Normalize temperatures
T_norm = (time_row - 300) / 10

# Get expected values for comparison

numerical = pd.read_csv('benchmark/expected_numerical.csv', header=None, names=['x', 'y'])
analytical = pd.read_csv('benchmark/expected_analytical.csv', header=None, names=['x', 'y'])


# Plot
plt.figure(figsize=(8, 4))
plt.plot(np.linspace(0.01, 0.96, 20), T_norm, marker='o', label = 'foamForNuclear')
plt.plot(numerical['x'], numerical['y'], label='Expected Numerical', linestyle='--')
plt.plot(analytical['x'], analytical['y'], label='Expected Analytical', linestyle='-.')
plt.xlabel('Probe Index')
plt.ylabel('Normalized Temperature (T - 300) / (T_s - 300))')
plt.title('Normalized Temperature at t = 8.0s')
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig("temperaturePlot.png")
