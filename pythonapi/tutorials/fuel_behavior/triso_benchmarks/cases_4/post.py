"""
Post-processing for PyC TRISO benchmark cases 4a–4d.

Plots σ_θ at the inner face of each layer and σ_r at the IPyC/SiC interface
vs fast fluence. If a reference SVG exists in reference_figures/ a comparison
plot is produced against the digitised TECDOC curves.

Usage
-----
    python post.py                      # all cases
    python post.py case_4a case_4c      # subset
"""

import sys
import os
import glob
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from _config import CASES as _CASES, PYC_FAST_FLUX

HERE      = os.path.dirname(os.path.abspath(__file__))
REFDIR    = os.path.join(HERE, "..", "reference_figures")
EXTRACTOR = os.path.join(REFDIR, "extractSVGPlotData.py")

# Import axis scales from central metadata file
sys.path.insert(0, REFDIR)
from axes import AXES, get_series_style, OFFBEAT_STYLE

um        = 1e-6
MPa       = 1e6
fast_flux = PYC_FAST_FLUX["fastFlux"][0]

CRP6_RANGES = {
    ("case_4a", "IPyC"): ( 925,  970),
    ("case_4a", "SiC"):  (-850, -775),
    ("case_4b", "IPyC"): ( -25,  -25),
    ("case_4b", "SiC"):  ( 138,  142),
    ("case_4c", "IPyC"): (  25,   27),
    ("case_4c", "SiC"):  (  83,   92),
    ("case_4d", "IPyC"): (  25,   35),
    ("case_4d", "SiC"):  (  71,   88),
}

# Maps (case_id, series_key) → SVG stem in reference_figures/
# series_key is a layer name ("IPyC", "SiC") or "Interface_Radial"
REFS = {
    ("case_4a", "IPyC"):             "case_4a_IPyC_tangential_stress",
    ("case_4a", "SiC"):              "case_4a_SiC_tangential_stress",
    ("case_4a", "Interface_Radial"): "case_4a_interface_radial_stress",
    ("case_4b", "IPyC"):             "case_4b_IPyC_tangential_stress",
    ("case_4b", "SiC"):              "case_4b_SiC_tangential_stress",
    ("case_4b", "Interface_Radial"): "case_4b_interface_radial_stress",
    ("case_4c", "IPyC"):             "case_4c_IPyC_tangential_stress",
    ("case_4c", "SiC"):              "case_4c_SiC_tangential_stress",
    ("case_4c", "Interface_Radial"): "case_4c_interface_radial_stress",
    ("case_4d", "IPyC"):             "case_4d_IPyC_tangential_stress",
    ("case_4d", "SiC"):              "case_4d_SiC_tangential_stress",
    ("case_4d", "Interface_Radial"): "case_4d_interface_radial_stress",
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


def post(case_id: str):
    cfg        = _CASES[case_id]
    layers_cfg = cfg["layers"]
    csvs       = _all_csvs(case_id)

    layer_names = [name for name, _, _ in layers_cfg]
    r = cfg["inner_radius"]
    layer_r_in = []
    for _, thickness, _ in layers_cfg:
        layer_r_in.append(r * um)
        r += thickness

    r_interface = (cfg["inner_radius"] + layers_cfg[0][1]) * um

    times     = []
    st_series = {name: [] for name in layer_names}
    sr_iface  = []

    for csv_path in csvs:
        t  = float(os.path.basename(os.path.dirname(csv_path)))
        df = pd.read_csv(csv_path)
        r_sim  = df["x"].to_numpy()
        st_sim = df["sigma_3"].to_numpy()
        sr_sim = df["sigma_0"].to_numpy()

        times.append(t)
        for name, r_in in zip(layer_names, layer_r_in):
            idx = np.argmin(np.abs(r_sim - r_in))
            st_series[name].append(st_sim[idx] / MPa)

        idx_iface = np.argmin(np.abs(r_sim - r_interface))
        sr_iface.append(sr_sim[idx_iface] / MPa)

    times   = np.array(times)
    fluence = fast_flux * times * 1e4

    all_series = {**st_series, "Interface_Radial": sr_iface}

    figures_dir = os.path.join(HERE, "figures")
    os.makedirs(figures_dir, exist_ok=True)

    # ── tangential hoop stress at inner face of each layer ──────────────────
    fig, ax = plt.subplots(figsize=(7, 4))
    fig.suptitle(f"TRISO {case_id} — inner-face σ_θ vs fluence", fontsize=12)
    for name in layer_names:
        ax.plot(fluence / 1e25, st_series[name], label=name, lw=1.5)
    ax.set_xlabel("Fast fluence  [10²⁵ n/m²]")
    ax.set_ylabel("σ_θ at inner face  [MPa]")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out = os.path.join(figures_dir, f"{case_id}_tangential.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")

    # ── radial stress at IPyC/SiC interface ─────────────────────────────────
    fig2, ax2 = plt.subplots(figsize=(7, 4))
    fig2.suptitle(f"TRISO {case_id} — IPyC/SiC interface σ_r vs fluence", fontsize=12)
    ax2.plot(fluence / 1e25, sr_iface, lw=1.5)
    ax2.set_xlabel("Fast fluence  [10²⁵ n/m²]")
    ax2.set_ylabel("σ_r at IPyC/SiC interface  [MPa]")
    ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    out2 = os.path.join(figures_dir, f"{case_id}_interface_radial.png")
    fig2.savefig(out2, dpi=150, bbox_inches="tight")
    print(f"Saved {out2}")

    # ── comparison against digitised TECDOC reference curves ────────────────
    for (ref_case, series_key), svg_stem in REFS.items():
        if ref_case != case_id or series_key not in all_series:
            continue
        xlsx_path = os.path.join(REFDIR, svg_stem + ".xlsx")
        if not os.path.isfile(xlsx_path):
            svg_path = os.path.join(REFDIR, svg_stem + ".svg")
            if not os.path.isfile(svg_path):
                print(f"No Excel or SVG reference found, skipping {case_id}/{series_key}")
                continue
            subprocess.run(["python", EXTRACTOR, svg_path, REFDIR], check=True)

        axes = AXES[svg_stem]
        x0, x1 = axes["xlim"]
        y0, y1 = axes["ylim"]
        xl = pd.read_excel(xlsx_path, sheet_name=None)

        ylabel = ("σ_r at IPyC/SiC interface  [MPa]" if series_key == "Interface_Radial"
                  else f"σ_θ at {series_key} inner face  [MPa]")

        fig_ref, ax_ref = plt.subplots(figsize=(8, 5))
        for curve_name, df in xl.items():
            x_real = df["x_norm"] * (x1 - x0) + x0
            y_real = df["y_norm"] * (y1 - y0) + y0
            ax_ref.plot(x_real, y_real, **get_series_style(curve_name),
                        label=curve_name)

        # ax_ref.plot(fluence / 1e25, all_series[series_key], **OFFBEAT_STYLE, label="OFFBEAT")
        ax_ref.set_xlim(axes["xlim"])
        ax_ref.set_ylim(axes["ylim"])
        if "ytick_step" in axes:
            ax_ref.set_yticks(np.arange(y0, y1 + 1, axes["ytick_step"]))
        ax_ref.set_xlabel("Fast fluence  [10²⁵ n/m²]")
        ax_ref.set_ylabel(ylabel)
        ax_ref.set_title(f"TRISO {case_id} — {series_key}")
        ax_ref.legend(loc="lower left")
        ax_ref.grid(True, alpha=0.3)
        out_ref = os.path.join(figures_dir, f"{case_id}_{series_key}_comparison.png")
        fig_ref.savefig(out_ref, dpi=150, bbox_inches="tight")
        print(f"Saved {out_ref}")

    # ── terminal summary ─────────────────────────────────────────────────────
    t_last_days = times[-1] / 86400
    print(f"\n{case_id}  —  end of life  (t = {t_last_days:.0f} days,"
          f"  fluence = {fluence[-1]/1e25:.2f}×10²⁵ n/m²)")
    print(f"  {'layer':>6}  {'σ_θ [MPa]':>12}  {'CRP-6 [MPa]':>16}")
    print(f"  {'------':>6}  {'----------':>12}  {'----------------':>16}")
    for name in layer_names:
        val = st_series[name][-1]
        key = (case_id, name)
        if key in CRP6_RANGES:
            lo, hi = CRP6_RANGES[key]
            marker = "✓" if lo - 0.5 <= val <= hi + 0.5 else "✗"
            print(f"  {name:>6}  {val:>12.1f}  [{lo:>6}, {hi:>6}]  {marker}")
        else:
            print(f"  {name:>6}  {val:>12.1f}")
    print(f"  {'σ_r iface':>9}  {sr_iface[-1]:>12.1f}")


if __name__ == "__main__":
    case_ids = sys.argv[1:] if len(sys.argv) > 1 else list(_CASES.keys())
    for case_id in case_ids:
        post(case_id)
    plt.show()
