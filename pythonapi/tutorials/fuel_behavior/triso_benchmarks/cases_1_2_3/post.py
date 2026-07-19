"""
Post-processing: compare Offbeat simulation to analytical Lamé solution.

Usage
-----
    python post.py                    # all cases
    python post.py case_2             # single case
    python post.py case_1 case_3      # subset

The script reads the last available radialProfile snapshot from the case
directory and overlays the analytical solution from analytical.py.

CSV column mapping (line_T_sigma.csv):
  x        : radial coordinate [m]
  sigma_0  : σ_xx = σ_rr  (radial stress)      [Pa]
  sigma_3  : σ_yy = σ_θθ  (tangential stress)  [Pa]
"""

import sys
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from analytical import analytical
from _config import CASES as _CASES

um = 1e-6   # µm → m
MPa = 1e6   # Pa → MPa


def _last_csv(case_id: str):
    pattern = os.path.join(
        case_id, "postProcessing", "radialProfile", "*", "line_T_sigma.csv"
    )
    files = sorted(glob.glob(pattern), key=lambda p: float(os.path.basename(os.path.dirname(p))))
    if not files:
        raise FileNotFoundError(f"No radialProfile CSV found in {case_id}/")
    return files[-1]


def post(case_id: str = "case_1"):
    cfg = _CASES[case_id]
    layers_cfg = cfg["layers"]

    # ------------------------------------------------------------------ #
    # Simulation data
    # ------------------------------------------------------------------ #
    csv_path = _last_csv(case_id)
    df = pd.read_csv(csv_path)
    r_sim   = df["x"].to_numpy()          # [m]
    sr_sim  = df["sigma_0"].to_numpy()    # σ_rr  [Pa]
    st_sim  = df["sigma_3"].to_numpy()    # σ_θθ  [Pa]

    # ------------------------------------------------------------------ #
    # Analytical solution
    # ------------------------------------------------------------------ #
    ana = analytical(case_id)

    # ------------------------------------------------------------------ #
    # Plot
    # ------------------------------------------------------------------ #
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=False)
    fig.suptitle(f"TRISO benchmark — {case_id}", fontsize=12)

    colors_ana = ["tab:blue", "tab:orange"]

    for ax, (key, label) in zip(axes, [("sigma_r", "Radial stress σ_r"), ("sigma_t", "Tangential stress σ_θ")]):
        # Simulation — step plot reflects cell-constant stress field
        sim_values = sr_sim if key == "sigma_r" else st_sim
        ax.step(r_sim / um, sim_values / MPa, where="mid",
                color="black", lw=1.2, label="Offbeat (sim)", zorder=5)

        # Analytical (per layer)
        for i, layer in enumerate(ana):
            lbl = f"Analytical — {layer['name']}" if len(ana) > 1 else "Analytical"
            ax.plot(layer["r"] / um, layer[key] / MPa,
                    color=colors_ana[i % len(colors_ana)], lw=1.8, label=lbl)

        # Layer boundaries
        r_bounds = [cfg["inner_radius"]] + [lay[1] for lay in layers_cfg]
        for rb in r_bounds:
            ax.axvline(rb, color="gray", lw=0.8, ls="--")

        ax.set_xlabel("r  [μm]")
        ax.set_ylabel("Stress  [MPa]")
        ax.set_title(label)
        ax.legend(fontsize=8)
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    os.makedirs("figures", exist_ok=True)

    # Print reference-point comparison (inner face of each layer)
    hdr = (f"\n{'':>6}  {'r [μm]':>8}"
           f"  {'σ_r sim':>10}  {'σ_r ana':>10}  {'err_r %':>8}"
           f"  {'σ_t sim':>10}  {'σ_t ana':>10}  {'err_t %':>8}  [MPa]")
    print(hdr)
    for layer in ana:
        r_check = layer["r"][0]
        idx_s = np.argmin(np.abs(r_sim - r_check))
        sr_s = sr_sim[idx_s] / MPa;  sr_a = layer["sigma_r"][0] / MPa
        st_s = st_sim[idx_s] / MPa;  st_a = layer["sigma_t"][0] / MPa
        err_r = (sr_s - sr_a) / abs(sr_a) * 100 if sr_a != 0 else float("nan")
        err_t = (st_s - st_a) / abs(st_a) * 100 if st_a != 0 else float("nan")
        print(f"  {layer['name']:>6}  {r_check/um:8.1f}"
              f"  {sr_s:10.2f}  {sr_a:10.2f}  {err_r:+8.2f}"
              f"  {st_s:10.2f}  {st_a:10.2f}  {err_t:+8.2f}")

    if len(ana) > 1:
        Pi = ana[0].get("Pi")
        if Pi is not None:
            # Π is contact pressure (positive=compressive); σ_r at interface = −Π
            print(f"\n  σ_r at interface = {-Pi/MPa:.2f} MPa  (ref: −18.8 MPa)")

    out_path = os.path.join("figures", f"verification_{case_id}.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"\nFigure saved to {out_path}")


if __name__ == "__main__":
    from _config import CASES
    case_ids = sys.argv[1:] if len(sys.argv) > 1 else list(CASES.keys())
    for case_id in case_ids:
        post(case_id)
    plt.show()
