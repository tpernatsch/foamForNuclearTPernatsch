import sys
import matplotlib.pyplot as plt
import numpy as np
# import scienceplots
# plt.style.use("ieee")

# Mode: ESI (default) oppure FOUNDATION
mode = "ESI"
if len(sys.argv) >= 2:
    m = sys.argv[1].strip().upper()
    if m in ("ESI", "FOUNDATION"):
        mode = m
    else:
        print(f"[WARN] Argument not recognized '{sys.argv[1]}', using default ESI.")

# Se FOUNDATION -> heatFlux da 10/, altrimenti da 0/
hf_folder = "10" if mode == "FOUNDATION" else "0"

# OFFBEAT data
dataT = np.genfromtxt("./postProcessing/probes/0/T", comments="#")
dataTcool = np.genfromtxt("./postProcessing/probes/0/Tcoolant", comments="#")
dataHF = np.genfromtxt(f"./postProcessing/probes/{hf_folder}/heatFlux", comments="#")

dataT_TRACE = np.genfromtxt("./TRACE_Twall", skip_header=3)
dataTcool_TRACE = np.genfromtxt("./TRACE_T0", skip_header=3)
dataHF_TRACE = np.genfromtxt("./TRACE_qwall", skip_header=3)

i_node_TRACE = 5

# Interpolate TRACE data to match OFFBEAT time points
interp_T_TRACE = np.interp(dataT[:,0], dataT_TRACE[:,0], dataT_TRACE[:,i_node_TRACE])
interp_Tcool_TRACE = np.interp(dataTcool[:,0], dataTcool_TRACE[:,0], dataTcool_TRACE[:,i_node_TRACE])

# Compute RMSE for Temperature
rmse_T_wall = np.sqrt(np.mean((dataT[:,1] - interp_T_TRACE) ** 2))
rmse_T_coolant = np.sqrt(np.mean((dataTcool[:,1] - interp_Tcool_TRACE) ** 2))

# Write RMSE to file
with open("results", "w") as f:
    f.write(f"RMSE (OFFBEAT vs TRACE) - Wall Temperature = {rmse_T_wall:.2f} K\n")
#     f.write(f"RMSE (OFFBEAT vs TRACE) - Coolant Temperature = {rmse_T_coolant:.2f} K\n")

cm = 1/2.54
fig, ax = plt.subplots(1, figsize=(8*cm, 8*cm))

ax.plot(dataT_TRACE[:,0]/3600 , dataT_TRACE[:,i_node_TRACE], 
        label="TRACE - wall", color='lightgreen', ls='-', lw = 3, alpha=0.8)
ax.plot(dataT[:,0]/3600 , dataT[:,1],
        label="OFFBEAT - wall", color='darkgreen', ls='-', lw=1, zorder=3)

ax.plot(dataTcool_TRACE[:,0]/3600 , dataTcool_TRACE[:,i_node_TRACE], 
        label="TRACE - Na", color='lightblue', ls='-', lw = 3, alpha=0.8)
ax.plot(dataTcool[:,0]/3600 , dataTcool[:,1], 
        label="OFFBEAT - Na", color='b', ls='-', lw=1, zorder=3)

ax.minorticks_on()
ax.tick_params(axis='both', which='both', direction='in')
ax.legend(frameon=False, loc="upper right", fontsize=6, ncol=1)
ax.set_xlabel("Time (hours)")
ax.set_ylabel("Temperature (K)")
plt.tight_layout()
plt.savefig("plt_NaCoolantVsTRACE.png")

# ------------------- HEAT FLUX PLOT

cm = 1/2.54
fig, ax = plt.subplots(1, figsize=(8*cm, 8*cm))

ax.plot(dataHF_TRACE[:,0]/3600 , dataHF_TRACE[:,i_node_TRACE]/max(dataHF_TRACE[:,i_node_TRACE]), 
        label="TRACE", color='lightblue', ls='-', lw = 3, alpha=0.8)
ax.plot(dataHF[:,0]/3600 , dataHF[:,1]/max(dataHF[:,1]), 
        label="OFFBEAT", color='b', ls='-', lw=1)

ax.minorticks_on()
ax.tick_params(axis='both', which='both', direction='in')
ax.legend(frameon=False)
ax.set_xlabel("Time (hours)")
ax.set_ylabel(r"Heat Flux (W/m$^2$)")
plt.tight_layout()
plt.savefig("plt_NaCoolantVsTRACE_heatflux.png")
