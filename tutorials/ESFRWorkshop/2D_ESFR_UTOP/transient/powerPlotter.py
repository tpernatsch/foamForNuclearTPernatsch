import re
import matplotlib.pyplot as plt

# Storage
times = []
fuel_Tmax = []
clad_Tmax = []
fluid_Tmax = []
power_vals = []

# Temporary holders per time step
current_time = None
last_fuel_Tmax = None
last_clad_Tmax = None
last_fluid_Tmax = None
last_power = None

with open("log.GeN-Foam", "r") as log:
    for line in log:
        line = line.strip()

        # New time step
        time_match = re.match(r"Time\s*=\s*([0-9Ee\.\+-]+)", line)
        if time_match:
            if current_time is not None:
                # Append only if we have all values
                if last_fuel_Tmax is not None:
                    times.append(current_time)
                    fuel_Tmax.append(last_fuel_Tmax)
                    clad_Tmax.append(last_clad_Tmax)
                    fluid_Tmax.append(last_fluid_Tmax)
                    power_vals.append(last_power)
            current_time = float(time_match.group(1))
            last_fuel_Tmax = None
            last_clad_Tmax = None
            last_fluid_Tmax = None
            last_power = None

        # Temperature matches
        fuel_match = re.match(r"T\.nuclearFuelPin\.fuel\s+\(avg min max\)\s*=\s*\S+\s+\S+\s+(\S+)", line)
        clad_match = re.match(r"T\.nuclearFuelPin\.clad\s+\(avg min max\)\s*=\s*\S+\s+\S+\s+(\S+)", line)
        fluid_match = re.match(r"T\s+\(avg min max\)\s*=\s*\S+\s+\S+\s+(\S+)", line)

        if fuel_match:
            last_fuel_Tmax = float(fuel_match.group(1))
        if clad_match:
            last_clad_Tmax = float(clad_match.group(1))
        if fluid_match:
            last_fluid_Tmax = float(fluid_match.group(1))

        # Power line
        power_match = re.match(r"totalPower =\s*([0-9Ee\.\+-]+)", line)
        if power_match:
            last_power = float(power_match.group(1))

# Add last timestep
if current_time is not None and last_fuel_Tmax is not None:
    times.append(current_time)
    fuel_Tmax.append(last_fuel_Tmax)
    clad_Tmax.append(last_clad_Tmax)
    fluid_Tmax.append(last_fluid_Tmax)
    power_vals.append(last_power)

# --- Plotting ---

def save_plot(x, y, ylabel, title, filename, color):
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, marker='o', linestyle='-', color=color)
    plt.xlabel('Time [s]')
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"✅ Saved: {filename}")

# Save all plots
save_plot(times, fuel_Tmax, 'Max Fuel Temperature [K]', 'Maximum Fuel Temperature vs Time', 'fuel_Tmax_vs_time.png', 'darkorange')
save_plot(times, clad_Tmax, 'Max Cladding Temperature [K]', 'Maximum Cladding Temperature vs Time', 'clad_Tmax_vs_time.png', 'steelblue')
save_plot(times, fluid_Tmax, 'Max Fluid Temperature [K]', 'Maximum Fluid Temperature vs Time', 'fluid_Tmax_vs_time.png', 'crimson')
save_plot(times, power_vals, 'Power [W]', 'Reactor Power vs Time', 'power_vs_time.png', 'seagreen')
