#==============================================================================*
# Imports

import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate


#==============================================================================*

m_dry = 20e3 # kg
T_exh = 2200 # K
massFlow = 30 # kg/s
# LEO 250km to LMO 100km
deltaV_tot = 2*(2440 + 679 + 145 + 676) + 120 # m/s

M_H2 = 2.018e-3 # kg/mol
gamma = 1.4
R = 8.314 # J/K/mol
g = 9.8 # m/s2

v_exh = np.sqrt(2*gamma*R*T_exh / ((gamma-1) * M_H2))

m_prop = m_dry * (np.exp(deltaV_tot/v_exh) - 1)

m_wet = m_prop + m_dry

thrust = massFlow * v_exh

Isp = v_exh / g

print(tabulate(
    [
        ["deltaV total", deltaV_tot, "m/s"],
        ["Dry mass", m_dry, "kg"],
        ["Propellant mass", m_prop, "kg"],
        ["Wet mass", m_wet, "kg"],
        ["Exhaust temperature", T_exh, "K"],
        ["Exhaust velocity", v_exh, "m/s"],
        ["Thrust", thrust/1000, "kN"],
        ["Isp", Isp, "s"],
    ],
    ["Parameter", "Value", "Unit"],
    tablefmt="simple_grid"
))


# Moon transfer
deltaV_tot = 2440 + 679
m_prop = m_wet * (1 - np.exp(-deltaV_tot/v_exh))
burnTime = m_prop/massFlow

print(tabulate(
    [
        ["Propellant mass for Moon transfer", m_prop, "kg"],
        ["Estimated time of burn", f"{burnTime:.3f}", "s"],
    ],
    ["Parameter", "Value", "Unit"],
    tablefmt="simple_grid"
))


deltaV = lambda t: v_exh*np.log(m_wet/(m_wet - massFlow*t))

tRange = np.linspace(0, burnTime, 100)
fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
ax.plot(tRange, deltaV(tRange))
ax.set_xlabel("Time [s]")
ax.set_ylabel("$\Delta v$ [m/s]")
ax.set_xlim(0)
ax.set_ylim(0)
ax.set_xticks([burnTime], minor=True)
ax.xaxis.grid(True, which='minor')
ax.set_yticks([deltaV_tot], minor=True)
ax.yaxis.grid(True, which='minor')
fig.tight_layout()
fig.savefig("fig_compute_deltaV.png")


#==============================================================================*
