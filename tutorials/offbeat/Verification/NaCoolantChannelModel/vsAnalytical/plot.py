import matplotlib.pyplot as plt
import numpy as np

# Radius of the heater
ro = 5e-3

# Mass flowrate
m_dot = 0.05

# Sodium cp (assumed constant) 
cp = 1275

# Na coolant inlet temperature 
Tin = 600

# Power factor, q'''(z) = P*sin(pi*z)
P = 3e8

def analytical(z):
    # Enthalpy rise deltaH = 1/m_dot * integr(q'' * pi * heatedDiam * dz)
    # q'' = ro/2 * q'''
    dH = (ro**2*P)/m_dot * (1-np.cos(np.pi*z))
    return Tin + dH/cp

# Get Temperature OFFBEAT surface
data = np.genfromtxt("postProcessing/patchSurface/10/T_patch.raw", comments="#")
z_OFFBEAT = data[:,2]
T_surface = data[:,-1]

# Get Temperature OFFBEAT coolant
dataNa = np.genfromtxt("postProcessing/patchSurface/10/Tcoolant_patch.raw", comments="#")
z_Na = dataNa[:,2]
T_Na = dataNa[:,-1]

# Compute RMSE
T_analytical = analytical(z_Na)
rmse = np.sqrt(np.mean((T_Na - T_analytical) ** 2))

# Write RMSE to file
with open("results", "w") as f:
    f.write(f"RMSE (OFFBEAT vs Analytical) = {rmse:.2f} K\n")

# Plotting
zeta = np.linspace(0, 1)
cm = 1/2.54
fig, ax = plt.subplots(1, figsize=(8*cm, 8*cm))
ax2 = ax.twinx()
ax2.plot(zeta, np.sin(np.pi*zeta), color="b", zorder=3)
ax.plot(zeta, analytical(zeta), color ="k", label="T Na, analytical", zorder=3)
ax.scatter(z_Na, T_Na, s=15, marker="o", 
           facecolors="none", edgecolors="g",
           label="T Na, OFFBEAT", zorder=2)

ax2.set_ylabel("Norm. Power Density (-)")
ax.legend(frameon = False, bbox_to_anchor=(0.5, 1.05), loc="center", ncol=2, fontsize=7)
ax.set_xlabel("Axial Coordinate (m)")
ax.set_ylabel("Temperature (K)")
ax.grid(which='major', color='#DDDDDD', linewidth=0.8, zorder=1)
ax.grid(which='minor', color='#DDDDDD', linestyle=':', linewidth=0.5, zorder=1)
ax.minorticks_on()
ax.tick_params(axis='both', which='both', direction='in')

plt.tight_layout()
plt.savefig("plt_NaCoolantVsAnalytical.png")