"""
Post-processing for IAEA HTGR FP diffusion benchmark cases 3a–3e.

Usage
-----
    python post.py                        # all cases
    python post.py case_3a case_3b        # subset
"""

import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from _config import CASES as _CASES

HERE   = os.path.dirname(os.path.abspath(__file__))
REFDIR = os.path.join(HERE, "..", "reference_figures")

s_per_h = 3600.0

REFERENCE_FIGS = {
    frozenset({"case_3a", "case_3b", "case_3c", "case_3d", "case_3e"}): dict(
        file=os.path.join(REFDIR, "case_3_fractional_Cs_0h_400h_1e-8_1_log.png"),
        xlim=(0, 400),
        ylim=(1e-8, 1),
        xtick_step=100,
        xlabel="Time  [h]",
        ylabel="Cs fractional release  [−]",
        yscale="log",
    ),
}

_CASE_COLORS = {
    "case_3a": "tab:purple",
    "case_3b": "tab:brown",
    "case_3c": "tab:pink",
    "case_3d": "tab:gray",
    "case_3e": "tab:olive",
}

_OVERLAY_STYLE = {
    "case_3a": dict(color="black", ls="--",              lw=2.5),
    "case_3b": dict(color="black", ls="-.",              lw=2.5),
    "case_3c": dict(color="black", ls=":",               lw=2.5),
    "case_3d": dict(color="black", ls=(0, (5, 2)),       lw=2.5),
    "case_3e": dict(color="black", ls=(0, (3, 1, 1, 1)), lw=2.5),
}


def _find_dat(case_id):
    files = glob.glob(os.path.join(
        case_id, "postProcessing", "fpRelease", "*", "writeFpRelease.dat"
    ))
    if not files:
        raise FileNotFoundError(f"No fpRelease dat under {case_id}/postProcessing/")
    files.sort(key=lambda p: float(os.path.basename(os.path.dirname(p))))
    return files[0]


def _read_dat(fpath):
    species, rows = [], []
    with open(fpath) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                parts = s.lstrip("#").split()
                if parts and parts[0] == "Time":
                    species = parts[1:]
                continue
            parts = s.split()
            if len(parts) >= 2:
                rows.append([float(v) for v in parts])
    arr = np.array(rows)
    return arr[:, 0], {sp: arr[:, i + 1] for i, sp in enumerate(species)}


def post(case_id):
    cfg = _CASES[case_id]
    expected = cfg.get("expected", {})
    try:
        fpath = _find_dat(case_id)
    except FileNotFoundError as e:
        print(f"  [{case_id}] skipped: {e}")
        return None

    times, data = _read_dat(fpath)
    times_h = times / s_per_h

    T_cfg = cfg["T"]
    T_str = (f"{T_cfg[0]:.0f}→{T_cfg[1]:.0f} K" if isinstance(T_cfg, tuple)
             else f"{T_cfg:.0f} K")

    n = len(data)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4), squeeze=False)
    fig.suptitle(f"FP release — {case_id}  (T = {T_str})", fontsize=11)

    print(f"\n{'─'*60}\n  {case_id}   T = {T_str}\n{'─'*60}")
    print(f"  {'Species':>8}  {'Final sim':>12}  {'Expected':>12}  {'Error %':>9}")

    for ax, sp in zip(axes[0], data):
        fracs = data[sp]
        ax.plot(times_h, fracs, color=_CASE_COLORS.get(case_id, "tab:blue"), lw=1.5, label="Offbeat")
        exp_val = expected.get(sp)
        final = fracs[-1]
        if exp_val is not None:
            ax.axhline(exp_val, color="tab:red", lw=1.2, ls="--", label=f"IAEA ref = {exp_val:.4g}")
            err = (final - exp_val) / exp_val * 100 if exp_val != 0 else float("nan")
            print(f"  {sp:>8}  {final:>12.6f}  {exp_val:>12.6f}  {err:>+9.2f}")
        else:
            print(f"  {sp:>8}  {final:>12.6f}  {'N/A':>12}")
        ax.set_xlim(0, cfg["endTime"] / s_per_h)
        ax.set_yscale("log")
        pos = fracs[fracs > 0]
        ax.set_ylim(pos.min() * 0.5 if len(pos) else 1e-8, 1.0)
        ax.set_xlabel("Time [h]")
        ax.set_ylabel(f"{sp} release fraction [−]")
        ax.set_title(sp)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    os.makedirs("figures", exist_ok=True)
    out = os.path.join("figures", f"fpRelease_{case_id}.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  Figure → {out}")
    return times, data


def post_overlay(case_data, ref):
    if not os.path.isfile(ref["file"]):
        print(f"  [overlay] not found: {ref['file']}")
        return
    is_log = ref.get("yscale") == "log"
    if is_log:
        log_ylim = (np.log10(ref["ylim"][0]), np.log10(ref["ylim"][1]))
        img_extent = [*ref["xlim"], *log_ylim]
        plot_ylim = log_ylim
    else:
        img_extent = [*ref["xlim"], *ref["ylim"]]
        plot_ylim = ref["ylim"]

    img = mpimg.imread(ref["file"])
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.imshow(img, extent=img_extent, aspect="auto", origin="upper", zorder=0)
    for case_id, (times, data) in sorted(case_data.items()):
        style = _OVERLAY_STYLE.get(case_id, dict(color="black", ls="--", lw=2.5))
        for fracs in data.values():
            y = np.log10(np.where(fracs > 0, fracs, np.nan)) if is_log else fracs
            ax.plot(times / s_per_h, y, **style, label=case_id, zorder=2)
    ax.set_xlim(ref["xlim"])
    ax.set_ylim(plot_ylim)
    ax.set_xticks(np.arange(ref["xlim"][0], ref["xlim"][1] + 1, ref["xtick_step"]))
    if is_log:
        yticks = np.arange(np.ceil(plot_ylim[0]), np.floor(plot_ylim[1]) + 1)
        ax.set_yticks(yticks)
        ax.set_yticklabels([f"$10^{{{int(t)}}}$" for t in yticks])
    ax.set_xlabel(ref["xlabel"])
    ax.set_ylabel(ref["ylabel"])
    ax.set_title("IAEA FP diffusion — cases 3a–3e (overlay)")
    ax.legend(fontsize=8, loc="lower right")
    os.makedirs("figures", exist_ok=True)
    stem = os.path.splitext(os.path.basename(ref["file"]))[0]
    out = os.path.join("figures", f"{stem}_overlay.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  Overlay → {out}")


if __name__ == "__main__":
    case_ids = sys.argv[1:] if len(sys.argv) > 1 else list(_CASES.keys())
    results = {}
    for cid in case_ids:
        r = post(cid)
        if r is not None:
            results[cid] = r
    for ref_cases, ref in REFERENCE_FIGS.items():
        available = ref_cases & set(results)
        if available:
            post_overlay({cid: results[cid] for cid in available}, ref)
    plt.show()
