"""
Post-processing for PyC TRISO benchmark cases 5–7.

Cases 5–7 have no closed-form analytical solution; the IAEA TECDOC-1674
benchmark compares codes against each other (Figs 9.13–9.15).

Two plots are produced per case:
  1. Radial stress profile at the final time step.
  2. Max tangential stress at the inner surface of each layer vs fast fluence
     (matches the TECDOC benchmark comparison metric).

If a reference SVG exists in reference_figures/ for a given case/layer,
a comparison plot against the digitised TECDOC curves is also produced.
Axis scales are read from reference_figures/axes.py.

Usage
-----
    python post.py                    # all cases
    python post.py case_5 case_6      # subset
"""

import sys
import os
import glob
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from _config import CASES as _CASES, PYC_FAST_FLUX

HERE      = os.path.dirname(os.path.abspath(__file__))
REFDIR    = os.path.join(HERE, "..", "reference_figures")
EXTRACTOR = os.path.join(REFDIR, "extractSVGPlotData.py")

sys.path.insert(0, REFDIR)
from axes import AXES, get_series_style, OFFBEAT_STYLE

um   = 1e-6
MPa  = 1e6
fast_flux = PYC_FAST_FLUX["fastFlux"][0]

CRP6_RANGES = {
    ("case_5", "IPyC"): ( 40,  58),
    ("case_5", "SiC"):  (-56, -28),
    ("case_6", "IPyC"): ( 27,  38),
    ("case_6", "SiC"):  ( 28,  48),
    ("case_7", "IPyC"): ( 37,  50),
    ("case_7", "SiC"):  ( 10,  25),
}

# Maps (case_id, layer_name) → SVG stem in reference_figures/
REFS = {
    ("case_5", "IPyC"): "case_5_IPyC_tangential_stress",
    ("case_5", "SiC"):  "case_5_SiC_tangential_stress",
    ("case_6", "IPyC"): "case_6_IPyC_tangential_stress",
    ("case_6", "SiC"):  "case_6_SiC_tangential_stress",
    ("case_7", "IPyC"): "case_7_IPyC_tangential_stress",
    ("case_7", "SiC"):  "case_7_SiC_tangential_stress",
}


HERE = os.path.dirname(os.path.abspath(__file__))


def _all_csvs(case_id: str):
    """Return all radialProfile CSVs sorted by time."""
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


def post(case_id: str):
    cfg        = _CASES[case_id]
    layers_cfg = cfg["layers"]
    csvs       = _all_csvs(case_id)

    # ------------------------------------------------------------------ #
    # Build time series: max σ_t at inner surface of each layer
    # ------------------------------------------------------------------ #
    layer_names = [name for name, _, _ in layers_cfg]
    r = cfg["inner_radius"]
    layer_r_in = []
    for _, thickness, _ in layers_cfg:
        layer_r_in.append(r * um)
        r += thickness

    times      = []
    st_series  = {name: [] for name in layer_names}

    for csv_path in csvs:
        t = float(os.path.basename(os.path.dirname(csv_path)))
        df = pd.read_csv(csv_path)
        r_sim  = df["x"].to_numpy()
        st_sim = df["sigma_3"].to_numpy()   # σ_θθ [Pa]

        times.append(t)
        for name, r_in in zip(layer_names, layer_r_in):
            idx = np.argmin(np.abs(r_sim - r_in))
            st_series[name].append(st_sim[idx] / MPa)

    times = np.array(times)
    fluence = fast_flux * times * 1e4   # n/cm²/s × s × cm²/m² → n/m²

    # ------------------------------------------------------------------ #
    # Plot 1: radial stress profile at final time step
    # ------------------------------------------------------------------ #
    df_last = pd.read_csv(csvs[-1])
    r_sim   = df_last["x"].to_numpy()
    sr_sim  = df_last["sigma_0"].to_numpy()
    st_sim  = df_last["sigma_3"].to_numpy()

    fig1, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=False)
    t_last_days = times[-1] / 86400
    fig1.suptitle(f"TRISO {case_id} — radial profile at t = {t_last_days:.0f} days", fontsize=12)

    for ax, (sig, label) in zip(axes, [(sr_sim, "Radial σ_r"), (st_sim, "Tangential σ_θ")]):
        ax.step(r_sim / um, sig / MPa, where="mid", color="black", lw=1.2)
        r = cfg["inner_radius"]
        r_bounds = [r]
        for _, thickness, _ in layers_cfg:
            r += thickness
            r_bounds.append(r)
        for i, (rb, name) in enumerate(zip(r_bounds[:-1], layer_names)):
            rb2 = r_bounds[i + 1]
            ax.axvspan(rb, rb2, alpha=0.06, label=name)
        for rb in r_bounds:
            ax.axvline(rb, color="gray", lw=0.8, ls="--")
        ax.set_xlabel("r  [μm]")
        ax.set_ylabel("Stress  [MPa]")
        ax.set_title(label)
        ax.legend(fontsize=8)
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    figures_dir = os.path.join(HERE, "figures")
    os.makedirs(figures_dir, exist_ok=True)
    out1 = os.path.join(figures_dir, f"{case_id}_profile.png")
    fig1.savefig(out1, dpi=150, bbox_inches="tight")
    print(f"Saved {out1}")

    # ------------------------------------------------------------------ #
    # Plot 2: max σ_t at inner face of each layer vs fast fluence
    # ------------------------------------------------------------------ #
    fig2, ax2 = plt.subplots(figsize=(7, 4))
    fig2.suptitle(f"TRISO {case_id} — inner-face tangential stress vs fluence", fontsize=12)

    for name in layer_names:
        ax2.plot(fluence / 1e25, st_series[name], label=name, lw=1.5)

    ax2.set_xlabel("Fast fluence  [10²⁵ n/m²]")
    ax2.set_ylabel("σ_θ at inner face  [MPa]")
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    out2 = os.path.join(figures_dir, f"{case_id}_fluence.png")
    fig2.savefig(out2, dpi=150, bbox_inches="tight")
    print(f"Saved {out2}")

    # ------------------------------------------------------------------ #
    # Plot 3: comparison against digitised TECDOC reference curves
    # ------------------------------------------------------------------ #
    for (ref_case, layer_name), svg_stem in REFS.items():
        if ref_case != case_id or layer_name not in st_series:
            continue
        xlsx_path = os.path.join(REFDIR, svg_stem + ".xlsx")
        if not os.path.isfile(xlsx_path):
            svg_path = os.path.join(REFDIR, svg_stem + ".svg")
            if not os.path.isfile(svg_path):
                print(f"No Excel or SVG reference found, skipping {case_id}/{layer_name}")
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

        # ax_ref.plot(fluence / 1e25, st_series[layer_name], **OFFBEAT_STYLE, label="OFFBEAT")
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

    # ------------------------------------------------------------------ #
    # Terminal summary (final time step)
    # ------------------------------------------------------------------ #
    print(f"\n{case_id}  —  end of life  (t = {t_last_days:.0f} days,"
          f"  fluence = {fluence[-1]/1e25:.2f}×10²⁵ n/m²)")
    print(f"  {'layer':>6}  {'σ_θ [MPa]':>12}  {'CRP-6 [MPa]':>16}")
    print(f"  {'------':>6}  {'----------':>12}  {'----------------':>16}")
    for name in layer_names:
        val = st_series[name][-1]
        key = (case_id, name)
        if key in CRP6_RANGES:
            lo, hi = CRP6_RANGES[key]
            tol = 0.5
            marker = "✓" if lo - tol <= val <= hi + tol else "✗"
            print(f"  {name:>6}  {val:>12.1f}  [{lo:>6}, {hi:>6}]  {marker}")
        else:
            print(f"  {name:>6}  {val:>12.1f}")


if __name__ == "__main__":
    case_ids = sys.argv[1:] if len(sys.argv) > 1 else list(_CASES.keys())
    for case_id in case_ids:
        post(case_id)
    plt.show()
