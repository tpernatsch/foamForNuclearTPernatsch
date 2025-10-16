#=============================================================================*
# Imports

import sys
import matplotlib.pyplot as plt


#=============================================================================*

if (len(sys.argv) < 2):
    print("Usage: python plot.py folder/path ...")
    sys.exit(1)


#=============================================================================*
# Main

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

ax1.set_ylabel("power (W)")
ax2.set_ylabel("reactivity (pcm)")
ax3.set_ylabel("fuelTemperature (K)")
ax3.set_xlabel("time(s)")

linestyles = ["-", "--"]

folders = sys.argv[1:]

for filename, linestyle in zip(folders, linestyles):

    with open(filename, "r") as file:
        lines = file.readlines()
        powers = []
        times = []
        totRhos = []
        TFuels = []
        i = 0
        power = 0
        time = 0
        totRho = 0
        TFuel = 0
        for line in lines :
            if "totalPower" in line :
                power = float(line.split()[2])
            if "totalReactivity" in line :
                totRho = float(line.split()[2])
            if "TFuel = " in line :
                TFuel = float(line.split()[2])
            if "Time =" in line and not "ExecutionTime" in line :
                time = float(line.split()[2])
                powers.append(power)
                times.append(time)
                totRhos.append(totRho)
                TFuels.append(TFuel)
        del totRhos[0]
        del powers[0]
        del TFuels[0]
        t0 = times[0]
        for i in range(len(times)) :
            times[i] = times[i] - t0
        del times[-1]
        ax1.plot(times, powers, label=filename, linestyle=linestyle)
        ax2.plot(times, totRhos, label=filename, linestyle=linestyle)
        ax3.plot(times, TFuels, label=filename, linestyle=linestyle)

ax1.legend()
ax1.grid(True)
ax2.legend()
ax2.grid(True)
ax3.legend()
ax3.grid(True)

fig.tight_layout()
fig.savefig(f"fig_results_{'_'.join([folder.replace('/', '') for folder in folders])}.png")


#=============================================================================*
