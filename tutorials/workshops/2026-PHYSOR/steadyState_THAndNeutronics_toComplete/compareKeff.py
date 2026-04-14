import re

NOMINAL_KEFF = 1.00153

print("Feedback effects considered:")
print("  - Fuel temperature")
print("  - Cladding temperature")
print("  - Coolant density")
print()

last_keff = None
with open("log.GeN-Foam", "r") as log:
    for line in log:
        m = re.match(r"keff_\s*=\s*([0-9Ee\.\+-]+)", line.strip())
        if m:
            last_keff = float(m.group(1))

if last_keff is None:
    print("ERROR: no k-eff value found in log.GeN-Foam")
else:
    diff_pcm = (last_keff - NOMINAL_KEFF) * 1e5
    if abs(diff_pcm) < 1:
        diff_pcm = 0.0
    print(f"Last k-eff       : {last_keff:.6f}")
    print(f"Nominal k-eff    : {NOMINAL_KEFF:.6f}")
    print(f"Difference       : {diff_pcm:+.2f} pcm")
