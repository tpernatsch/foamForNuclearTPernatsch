import re
import matplotlib.pyplot as plt
import numpy as np

def parse_log(filepath):
    times = []
    fuel_Tmax = []
    clad_Tmax = []
    fluid_Tmax = []
    power_vals = []
    keff_vals = []

    current_time = None
    last_fuel_Tmax = None
    last_clad_Tmax = None
    last_fluid_Tmax = None
    last_power = None
    last_keff = None

    with open(filepath, "r") as log:
        for line in log:
            line = line.strip()

            time_match = re.match(r"Time\s*=\s*([0-9Ee\.\+-]+)", line)
            if time_match:
                if current_time is not None and last_fuel_Tmax is not None:
                    times.append(current_time)
                    fuel_Tmax.append(last_fuel_Tmax)
                    clad_Tmax.append(last_clad_Tmax)
                    fluid_Tmax.append(last_fluid_Tmax)
                    power_vals.append(last_power)
                    keff_vals.append(last_keff)
                current_time = float(time_match.group(1))
                last_fuel_Tmax = None
                last_clad_Tmax = None
                last_fluid_Tmax = None
                last_power = None
                last_keff = None

            fuel_match = re.match(r"T\.nuclearFuelPin\.fuel\s+\(avg min max\)\s*=\s*\S+\s+\S+\s+(\S+)", line)
            clad_match = re.match(r"T\.nuclearFuelPin\.clad\s+\(avg min max\)\s*=\s*\S+\s+\S+\s+(\S+)", line)
            fluid_match = re.match(r"T\s+\(avg min max\)\s*=\s*\S+\s+\S+\s+(\S+)", line)

            if fuel_match:
                last_fuel_Tmax = float(fuel_match.group(1))
            if clad_match:
                last_clad_Tmax = float(clad_match.group(1))
            if fluid_match:
                last_fluid_Tmax = float(fluid_match.group(1))

            power_match = re.match(r"power:\s*([0-9Ee\.\+-]+)", line)
            if power_match:
                last_power = float(power_match.group(1))

            keff_match = re.match(r"keff_\s*=\s*([0-9Ee\.\+-]+)", line)
            if keff_match:
                last_keff = float(keff_match.group(1))

    if current_time is not None and last_fuel_Tmax is not None:
        times.append(current_time)
        fuel_Tmax.append(last_fuel_Tmax)
        clad_Tmax.append(last_clad_Tmax)
        fluid_Tmax.append(last_fluid_Tmax)
        power_vals.append(last_power)
        keff_vals.append(last_keff)

    return {
        "time": np.array(times),
        "fuel": np.array(fuel_Tmax),
        "clad": np.array(clad_Tmax),
        "fluid": np.array(fluid_Tmax),
        "power": np.array(power_vals),
        "keff": np.array(keff_vals)
    }

# Interpolation helper
def interpolate_old_to_new(new_times, old_times, old_values):
    return np.interp(new_times, old_times, old_values)

# Plotting function
def plot_difference(time, diff, ylabel, title, filename, color):
    plt.figure(figsize=(10, 6))
    plt.plot(time, diff, marker='o', linestyle='-', color=color)
    plt.xlabel('Time [s]')
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"✅ Saved: {filename}")

# Load data
new_data = parse_log("log.GeN-Foam")
old_data = parse_log("log.oldResults")

# Time alignment
t_new = new_data["time"]

# Interpolate old data onto new time points
def rel_diff(new, old_interp):
    return (new - old_interp) / old_interp

old_interp = {
    key: interpolate_old_to_new(t_new, old_data["time"], old_data[key])
    for key in ["fuel", "clad", "fluid", "power", "keff"]
}

# Compute differences
fuel_diff = rel_diff(new_data["fuel"], old_interp["fuel"])
clad_diff = rel_diff(new_data["clad"], old_interp["clad"])
fluid_diff = rel_diff(new_data["fluid"], old_interp["fluid"])
power_diff = rel_diff(new_data["power"], old_interp["power"])
keff_diff_pcm = (new_data["keff"] - old_interp["keff"]) * 1e5  # in pcm

# Plot
plot_difference(t_new, fuel_diff*100, 'Relative Difference [%]', 'Fuel Tmax Relative Difference', 'fuel_diff.png', 'darkorange')
plot_difference(t_new, clad_diff*100, 'Relative Difference [%]', 'Cladding Tmax Relative Difference', 'clad_diff.png', 'steelblue')
plot_difference(t_new, fluid_diff*100, 'Relative Difference [%]', 'Fluid Tmax Relative Difference', 'fluid_diff.png', 'crimson')
plot_difference(t_new, keff_diff_pcm, 'Δk-eff [pcm]', 'Absolute Difference in k-eff [pcm]', 'keff_diff.png', 'purple')
