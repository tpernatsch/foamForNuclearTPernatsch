#==============================================================================*
# Imports

import matplotlib.pyplot as plt
import pandas as pd


#==============================================================================*
# Extract data

data = pd.read_csv("Controller.csv")

time = data['time']
massFlowCmd = data['pid.u_s']
massFlowDamped = data['pid.u_m']
massFlow = data['massFlow']
pumpMomentum = data['pumpMomentum']


#==============================================================================*
# Plot

fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
ax.plot(time, massFlowCmd, label="Cmd", color="black")
ax.plot(time, massFlow, label="Measure", color="tab:gray", alpha=0.5, ls='--')
ax.plot(time, massFlowDamped, label="Measure damped", color="black", ls='-.')
ax.set_xlim(0)
ax.set_ylim(0)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Mass flow rate [kg/s]')
ax.legend()
fig.tight_layout()
fig.savefig("fig_results_massFlow.png")
plt.close()

fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
ax.plot(time, pumpMomentum, color="black")
ax.set_xlim(0)
# ax.set_ylim(0)
ax.set_xlabel('Time [s]')
fig.tight_layout()
fig.savefig("fig_results_pumpMomentum.png")
plt.close()


#==============================================================================*
