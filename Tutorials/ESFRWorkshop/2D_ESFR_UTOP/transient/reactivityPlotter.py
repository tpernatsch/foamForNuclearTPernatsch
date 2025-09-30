import re
import matplotlib.pyplot as plt

# Storage
times = []
total_vals = []
external_vals = []
doppler_vals = []
fuel_exp_vals = []
clad_exp_vals = []
coolant_vals = []

# Temporary holders per time step
current_time = None
last_total = None
last_external = None
last_doppler = None
last_fuel_exp = None
last_clad_exp = None
last_coolant = None

with open("log.GeN-Foam", "r") as log:
    for line in log:
        line = line.strip()

        # New time step
        time_match = re.match(r"Time\s*=\s*([0-9Ee\.\+-]+)", line)
        if time_match:
            if current_time is not None:
                if last_total is not None:
                    times.append(current_time)
                    total_vals.append(last_total)
                    external_vals.append(last_external)
                    doppler_vals.append(last_doppler)
                    fuel_exp_vals.append(last_fuel_exp)
                    clad_exp_vals.append(last_clad_exp)
                    coolant_vals.append(last_coolant)
            current_time = float(time_match.group(1))
            last_total = None
            last_external = None
            last_doppler = None
            last_fuel_exp = None
            last_clad_exp = None
            last_coolant = None

        # Reactivity components
        total_match = re.match(r"totalReactivity\s*=\s*([0-9Ee\.\+-]+)\s*pcm", line)
        external_match = re.match(r"->\s*extReactivity\s*=\s*([0-9Ee\.\+-]+)\s*pcm", line)
        doppler_match = re.match(r"->\s*Doppler.*=\s*([0-9Ee\.\+-]+)\s*pcm", line)
        fuel_match = re.match(r"->\s*TFuel\s*=\s*([0-9Ee\.\+-]+)\s*pcm", line)
        clad_match = re.match(r"->\s*TClad\s*=\s*([0-9Ee\.\+-]+)\s*pcm", line)
        coolant_match = re.match(r"->\s*rhoCool\s*=\s*([0-9Ee\.\+-]+)\s*pcm", line)

        if total_match:
            last_total = float(total_match.group(1))
        if external_match:
            last_external = float(external_match.group(1))
        if doppler_match:
            last_doppler = float(doppler_match.group(1))
        if fuel_match:
            last_fuel_exp = float(fuel_match.group(1))
        if clad_match:
            last_clad_exp = float(clad_match.group(1))
        if coolant_match:
            last_coolant = float(coolant_match.group(1))

# Add final timestep
if current_time is not None and last_total is not None:
    times.append(current_time)
    total_vals.append(last_total)
    external_vals.append(last_external)
    doppler_vals.append(last_doppler)
    fuel_exp_vals.append(last_fuel_exp)
    clad_exp_vals.append(last_clad_exp)
    coolant_vals.append(last_coolant)

# --- Plotting ---
plt.figure(figsize=(12, 7))

plt.plot(times, total_vals, label='total', color='black', linewidth=2)
plt.plot(times, external_vals, label='external', color='green')
plt.plot(times, doppler_vals, label='doppler', color='red')
plt.plot(times, fuel_exp_vals, label='fuel expansion', color='orange')
plt.plot(times, clad_exp_vals, label='cladding expansion', color='blue')
plt.plot(times, coolant_vals, label='coolant density', color='purple')

plt.xlabel('Time [s]')
plt.ylabel('Reactivity [pcm]')
plt.title('Reactivity Components vs Time')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("reactivity_components_vs_time.png", dpi=300)
print("✅ Saved: reactivity_components_vs_time.png")
