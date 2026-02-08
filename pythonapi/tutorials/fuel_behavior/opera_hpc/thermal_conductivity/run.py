import numpy as np
import matplotlib.pyplot as plt
from case import setup_base_case, run_point

T_values  = range(500, 3000, 100)
Bu_values = [0, 20e3, 40e3]
results = {Bu: [] for Bu in Bu_values}

model, _ = setup_base_case()

for Bu in Bu_values:
    for T in T_values:
        print(f"Running (reused case): T={T}, Bu={Bu}")
        k_value = run_point(model, T, Bu)
        print(f"  -> k = {k_value}")
        results[Bu].append((T, k_value))

# -------------------------------
# Plot k vs T for each burnup
# -------------------------------
plt.figure(figsize=(6, 4), dpi=150)

for Bu, pairs in results.items():
    pairs = sorted(pairs, key=lambda x: x[0])  # sort by T
    T_list = [p[0] for p in pairs]
    k_list = [p[1] for p in pairs]

    plt.plot(T_list, k_list, marker='o', label=f"Bu = {Bu/1000:.0f} MWd/kg")

plt.xlabel("Temperature [K]")
plt.ylabel("Thermal conductivity k [W/mK]")
plt.title("k(T) for various burnup levels")
plt.grid(True, alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig("k_vs_T.png")

print("\nPlot saved as k_vs_T.png")
