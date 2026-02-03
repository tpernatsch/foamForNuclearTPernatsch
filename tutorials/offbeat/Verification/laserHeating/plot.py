import matplotlib.pyplot as plt
import numpy as np

def analytical(I0, alpha, dist):
    return I0 * np.exp(-alpha * dist)

# Load data
path = "./postProcessing/axialIntensity/1/axial_I.csv"
data = np.genfromtxt(path, skip_header=1, delimiter=",")

# Extract distances and numerical intensities
dist = 1e-3 - data[:, 0]  # Ensure correct sign convention
I_numerical = data[:, 1]

# Compute analytical values at numerical points
I_analytical = analytical(I0=100, alpha=1e4, dist=dist)

# Compute Error (Scalar)
mae = np.mean(np.abs(I_numerical - I_analytical))  # Mean Absolute Error (MAE)
rmse = np.sqrt(np.mean((I_numerical - I_analytical) ** 2))  # Root Mean Square Error (RMSE)

# print(f"Mean Absolute Error (MAE): {mae:.6f}")
# print(f"Root Mean Square Error (RMSE): {rmse:.6f}")

# Plotting
z = np.linspace(0, max(dist), 100)
I_analytical_full = analytical(I0=100, alpha=1e4, dist=z)

plt.figure()
plt.plot(z, I_analytical_full, "-r", label=r"$I(z)=I_0 \cdot \exp(-\alpha z)$")
plt.scatter(dist, I_numerical, marker="o", c="k", label="Numerical")

plt.xlabel("z (m)")
plt.ylabel("I (1/m² s)")
plt.legend(frameon=False)
plt.title("Photon Attenuation Solver")

plt.savefig("results.png")

# Save RMSE to a file
with open("results", "w") as f:
    f.write(f"RMSE wrt analytical solution: {rmse*100:.2f} %\n")
    if rmse*100 < 0.5:
        f.write(f"Test Passed")
    else:
        f.write("Test Failed")


