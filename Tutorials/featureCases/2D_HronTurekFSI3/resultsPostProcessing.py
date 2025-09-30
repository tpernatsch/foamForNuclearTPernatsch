import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def analyze_signal(time, signal, label='Signal', plot=True):
    # Only use time > 4s
    mask = time > 4.2
    time = time[mask]
    signal = signal[mask]

    # Find maxima and minima
    peaks_max, _ = find_peaks(signal)
    peaks_min, _ = find_peaks(-signal)

    if len(peaks_max) < 1 or len(peaks_min) < 1:
        print(f"{label} → Not enough peaks found.")
        return None, None, None

    max_vals = signal[peaks_max]
    min_vals = signal[peaks_min]

    # Amplitude and mean
    amp = (np.max(max_vals) - np.min(min_vals)) / 2
    mean = (np.max(max_vals) + np.min(min_vals)) / 2

    # Frequency
    if len(peaks_max) >= 2:
        periods = np.diff(time[peaks_max])
        freq = 1 / np.mean(periods)
    else:
        freq = None

    # Plot
    if plot:
        plt.figure()
        plt.plot(time, signal, label=label)
        plt.plot(time[peaks_max], signal[peaks_max], 'ro', label='Maxima')
        plt.plot(time[peaks_min], signal[peaks_min], 'bo', label='Minima')
        plt.axhline(mean, color='k', linestyle='--', label='Mean')
        plt.xlabel("Time (s)")
        plt.ylabel("Displacement")
        plt.title(f"{label} Analysis")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    if freq:
        print(f"{label} → Amplitude: {amp:.6f}, Mean: {mean:.6f}, Frequency: {freq:.3f} Hz")
    else:
        print(f"{label} → Amplitude: {amp:.6f}, Mean: {mean:.6f}, Frequency: N/A")

    return amp, mean, freq

# === Step 1: Load data from file ===
file_path1 = "displacementHistory.dat"
try:
    data1 = np.loadtxt(file_path1)
except Exception as e:
    print(f"Error reading file: {e}")
    exit()

time_all1 = data1[:, 0]
x_disp1 = data1[:, 1]
y_disp1 = data1[:, 2]

# === Step 2: Analyze ===
amp_x, mean_x, freq_x = analyze_signal(time_all1, x_disp1, label='X Displacement_ffn', plot = False)
amp_y, mean_y, freq_y = analyze_signal(time_all1, y_disp1, label='Y Displacement_ffn', plot = False)


# === Plot: X displacement comparison ===
plt.figure()
plt.plot(time_all1, x_disp1, label='ffn')

# Add horizontal dashed lines
x_base = -2.69e-3
x_upper = (-2.69 + 2.53) * 1e-3
x_lower = (-2.69 - 2.53) * 1e-3
for val in [x_base, x_upper, x_lower]:
    plt.axhline(y=val, color='k', linestyle='--', linewidth=1.5)

plt.xlabel("Time (s)")
plt.ylabel("Displacement")
plt.title("X displacement Analysis")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("xDispComp.png")


# === Plot: Y displacement comparison ===
plt.figure()
plt.plot(time_all1, y_disp1, label='f4n')

# Add horizontal dashed lines
y_base = 1.5e-3
y_upper = (1.5 + 34) * 1e-3
y_lower = (1.5 - 34) * 1e-3
for val in [y_base, y_upper, y_lower]:
    plt.axhline(y=val, color='k', linestyle='--', linewidth=1.5)

plt.xlabel("Time (s)")
plt.ylabel("Displacement")
plt.title("Y displacement Analysis")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("yDispComp.png")

