

import matplotlib.pyplot as plt
import pandas as pd
import sys

# Check input command line
if (len(sys.argv) != 2):
    print(f"Usage: python3 {sys.argv[0]} path/to/GeN-Foam.dat")
    sys.exit()

# Extract filename
filename = sys.argv[1]

# Extract data
data = pd.read_csv(filename, sep=';')

# Plot
plt.plot(data['time(s)'], data['power(W)'])
plt.xlabel("Time [s]")
plt.ylabel("Power [W]")
plt.grid(True)

plt.savefig("powerEvolution.png")

plt.show()
