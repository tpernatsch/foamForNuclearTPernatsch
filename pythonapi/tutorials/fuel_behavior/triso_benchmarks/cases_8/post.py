"""
Post-processing for PyC TRISO benchmark case 8 (temperature cycling).

Plots σ_θ at the inner face of each layer vs fast fluence, showing the
sawtooth stress evolution over 10 irradiation cycles. If reference SVGs
exist in reference_figures/ a comparison plot against the digitised TECDOC
curves is also produced. Axis scales are read from reference_figures/axes.py.

Usage
-----
    python post.py
"""

import os
import sys
import glob
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from _config import CASES as _CASES, PYC_FAST_FLUX

HERE      = os.path.dirname(os.path.abspath(__file__))
REFDIR    = os.path.join(HERE, "..", "reference_figures")
EXTRACTOR = os.path.join(REFDIR, "extractSVGPlotData.py")

sys.path.insert(0, REFDIR)
from axes import AXES, get_series_style, OFFBEAT_STYLE

um        = 1e-6
MPa       = 1e6
fast_flux = PYC_FAST_FLUX["fastFlux"][0]   # n/cm²/s

# Maps layer_name → SVG stem in reference_figures/
REFS = {
    "IPyC": "case_8_IPyC_tangential_stress",
    "SiC":  "case_8_SiC_tangential_stress",
}


def _all_csvs(case_id: str):
    pattern = os.path.join(
        HERE, case_id, "postProcessing", "radialProfile", "*", "line_T_sigma.csv"
    )
    files = sorted(
        glob.glob(pattern),
        key=lambda p: float(os.path.basename(os.path.dirname(p)))
    )
    if not files:
        raise FileNotFoundError(f"No radialProfile CSV found in {case_id}/")
    return files


def post(case_id: str = "case_8"):
    cfg        = _CASES[case_id]
    layers_cfg = cfg["layers"]
    csvs       = _all_csvs(case_id)

    layer_names = [name for name, _, _ in layers_cfg]
    r = cfg["inner_radius"]
    layer_r_in = []
    for _, thickness, _ in layers_cfg:
        layer_r_in.append(r * um)
        r += thickness

    times     = []
    st_series = {name: [] for name in layer_names}

    for csv_path in csvs:
        t = float(os.path.basename(os.path.dirname(csv_path)))
        df = pd.read_csv(csv_path)
        r_sim  = df["x"].to_numpy()
        st_sim = df["sigma_3"].to_numpy()

        times.append(t)
        for name, r_in in zip(layer_names, layer_r_in):
            idx = np.argmin(np.abs(r_sim - r_in))
            st_series[name].append(st_sim[idx] / MPa)

    times   = np.array(times)
    fluence = fast_flux * times * 1e4   # n/m²

    cycle_fluences = [0.30, 0.60, 0.90, 1.20, 1.50, 1.80, 2.10, 2.40, 2.70, 3.00]

    fig, ax = plt.subplots(figsize=(9, 4))
    fig.suptitle(f"TRISO {case_id} — inner-face σ_θ vs fluence (cycling)", fontsize=12)

    for name in layer_names:
        ax.plot(fluence / 1e25, st_series[name], label=name, lw=1.5)

    for phi in cycle_fluences:
        ax.axvline(phi, color="gray", lw=0.6, ls="--")

    ax.set_xlabel("Fast fluence  [10²⁵ n/m²]")
    ax.set_ylabel("σ_θ at inner face  [MPa]")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    figures_dir = os.path.join(HERE, "figures")
    os.makedirs(figures_dir, exist_ok=True)

    out = os.path.join(figures_dir, f"{case_id}_fluence.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")

    # Comparison plots against digitised TECDOC reference curves
    for layer_name, svg_stem in REFS.items():
        if layer_name not in st_series:
            continue
        xlsx_path = os.path.join(REFDIR, svg_stem + ".xlsx")
        if not os.path.isfile(xlsx_path):
            svg_path = os.path.join(REFDIR, svg_stem + ".svg")
            if not os.path.isfile(svg_path):
                print(f"No Excel or SVG reference found, skipping {layer_name}")
                continue
            subprocess.run(["python", EXTRACTOR, svg_path, REFDIR], check=True)

        axes = AXES[svg_stem]
        x0, x1 = axes["xlim"]
        y0, y1 = axes["ylim"]
        xl = pd.read_excel(xlsx_path, sheet_name=None)

        fig_ref, ax_ref = plt.subplots(figsize=(8, 5))
        for curve_name, df in xl.items():
            x_real = df["x_norm"] * (x1 - x0) + x0
            y_real = df["y_norm"] * (y1 - y0) + y0
            ax_ref.plot(x_real, y_real, **get_series_style(curve_name),
                        label=curve_name)

        ax_ref.plot(fluence / 1e25, st_series[layer_name], **OFFBEAT_STYLE, label="OFFBEAT")
        ax_ref.set_xlim(axes["xlim"])
        ax_ref.set_ylim(axes["ylim"])
        if "ytick_step" in axes:
            ax_ref.set_yticks(np.arange(y0, y1 + 1, axes["ytick_step"]))
        ax_ref.set_xlabel("Fast fluence  [10²⁵ n/m²]")
        ax_ref.set_ylabel("σ_θ at inner face  [MPa]")
        ax_ref.set_title(f"TRISO {case_id} — {layer_name} inner-face σ_θ")
        ax_ref.legend(loc="lower left")
        ax_ref.grid(True, alpha=0.3)
        out_ref = os.path.join(figures_dir, f"{case_id}_{layer_name}_comparison.png")
        fig_ref.savefig(out_ref, dpi=150, bbox_inches="tight")
        print(f"Saved {out_ref}")

    t_last_days = times[-1] / 86400
    print(f"\n{case_id}  —  end of life  (t = {t_last_days:.0f} days,"
          f"  fluence = {fluence[-1]/1e25:.2f}×10²⁵ n/m²)")
    print(f"  {'layer':>6}  {'σ_t(r_in) [MPa]':>18}")
    for name in layer_names:
        print(f"  {name:>6}  {st_series[name][-1]:>18.1f}")


if __name__ == "__main__":
    post()
    plt.show()
