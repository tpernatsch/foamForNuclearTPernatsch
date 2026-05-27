"""
Post-processing: compare Offbeat FP-diffusion simulation to IAEA benchmark values.

Usage
-----
    python post.py                    # all cases
    python post.py case_4a            # single case
    python post.py case_4a case_4b    # subset

Individual release-fraction-vs-time plots are written to figures/.
When all cases that share a reference figure have been processed, an
overlay plot is also produced.
"""

import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from _config import CASES as _CASES, T_IRRAD_END

HERE   = os.path.dirname(os.path.abspath(__file__))
REFDIR = os.path.join(HERE, "..", "reference_figures")

s_per_h = 3600.0

# x-axis shows time after start of heating (reference figure convention)
REFERENCE_FIGS = {
    "heating": dict(
        file=os.path.join(REFDIR, "case_4_fractional_Cs_0h_400h_1e-8_1_log.png"),
        cases=frozenset({"case_4a", "case_4b", "case_4c", "case_4d"}),
        xlim=(0, 400),
        ylim=(1e-8, 1),
        xtick_step=100,
        xlabel="Time after start of heating  [h]",
        ylabel="Cs / Ag fractional release  [−]",
        yscale="log",
        time_shift=T_IRRAD_END,  # subtract irradiation duration (s)
    ),
}

_CASE_COLORS = {
    "case_4a": "tab:blue",
    "case_4b": "tab:orange",
    "case_4c": "tab:green",
    "case_4d": "tab:red",
}

_OVERLAY_STYLE = {
    "case_4a": dict(color="black", ls="--",              lw=2.5),
    "case_4b": dict(color="black", ls="-.",              lw=2.5),
    "case_4c": dict(color="black", ls=":",               lw=2.5),
    "case_4d": dict(color="black", ls=(0, (5, 2)),       lw=2.5),
}


def _find_dat(case_id: str) -> str:
    pattern = os.path.join(
        case_id, "postProcessing", "fpRelease", "*", "writeFpRelease.dat"
    )
    files = glob.glob(pattern)
    if not files:
        raise FileNotFoundError(
            f"No fpRelease dat file found under {case_id}/postProcessing/fpRelease/"
        )
    files.sort(key=lambda p: float(os.path.basename(os.path.dirname(p))))
    return files[0]


def _read_dat(fpath: str) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    species: list[str] = []
    rows: list[list[float]] = []

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

    if not rows:
        raise RuntimeError(f"No data rows found in {fpath}")

    arr = np.array(rows)
    times = arr[:, 0]
    data = {sp: arr[:, i + 1] for i, sp in enumerate(species)}
    return times, data


def post(case_id: str) -> tuple[np.ndarray, dict[str, np.ndarray]] | None:
    """Plot one case (irradiation + heating as separate panels) and return (times_s, data)."""
    cfg = _CASES[case_id]
    expected = cfg.get("expected", {})

    try:
        fpath = _find_dat(case_id)
    except FileNotFoundError as e:
        print(f"  [{case_id}] skipped: {e}")
        return None

    times, data = _read_dat(fpath)
    times_h = times / s_per_h
    t_irrad_h = T_IRRAD_END / s_per_h          # 12 000 h
    heat_end_h = (cfg["endTime"] - T_IRRAD_END) / s_per_h   # 200 or 400 h

    species_list = list(data.keys())
    n = len(species_list)

    T_heat = cfg["T_heat"]
    T_str = (f"{T_heat[0]:.0f}→{T_heat[1]:.0f} K" if len(T_heat) == 2
             else f"{T_heat[0]:.0f} K")
    color = _CASE_COLORS.get(case_id, "tab:blue")

    print(f"\n{'─'*60}")
    print(f"  {case_id}   T_heat = {T_str}")
    print(f"{'─'*60}")
    print(f"  {'Species':>8}  {'Final sim':>12}  {'Expected':>12}  {'Error %':>9}")

    # Two-column layout: irradiation | heating
    fig, axes = plt.subplots(n, 2, figsize=(10, 4 * n), squeeze=False)
    fig.suptitle(f"FP release — {case_id}  (irradiation + heating @ {T_str})", fontsize=11)

    for i, sp in enumerate(species_list):
        fracs = data[sp]
        exp_val = expected.get(sp)
        final = fracs[-1] if len(fracs) else float("nan")

        if exp_val is not None:
            err = (final - exp_val) / exp_val * 100 if exp_val != 0 else float("nan")
            print(f"  {sp:>8}  {final:>12.6g}  {exp_val:>12.6g}  {err:>+9.2f}")
        else:
            print(f"  {sp:>8}  {final:>12.6g}  {'N/A':>12}")

        pos = fracs[fracs > 0]
        ymin = pos.min() * 0.5 if len(pos) else 1e-8

        for col, (mask, xlabel, xlim, ylim) in enumerate([
            (times_h <= t_irrad_h,  "Time [h]",                        (0, t_irrad_h),  (ymin, 1.0)),
            (times_h >= t_irrad_h,  "Time after start of heating [h]", (0, heat_end_h), (1e-8, 1.0)),
        ]):
            ax = axes[i][col]
            x = times_h[mask] if col == 0 else times_h[mask] - t_irrad_h
            ax.plot(x, fracs[mask], color=color, lw=1.5, label="Offbeat")
            if exp_val is not None and col == 1:
                ax.axhline(exp_val, color="tab:red", lw=1.2, ls="--",
                           label=f"IAEA ref = {exp_val:.4g}")
            ax.set_xlim(*xlim)
            ax.set_yscale("log")
            ax.set_ylim(*ylim)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(f"{sp} release fraction [−]")
            phase = "irradiation" if col == 0 else "heating"
            ax.set_title(f"{sp} — {phase}")
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    os.makedirs("figures", exist_ok=True)
    out = os.path.join("figures", f"fpRelease_{case_id}.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  Figure → {out}")

    return times, data


def post_overlay(
    case_data: dict[str, tuple[np.ndarray, dict[str, np.ndarray]]],
    ref: dict,
):
    """Overlay all cases in case_data on the reference figure."""
    if not os.path.isfile(ref["file"]):
        print(f"  [overlay] reference figure not found: {ref['file']}")
        return

    is_log = ref.get("yscale") == "log"
    time_shift_h = ref.get("time_shift", 0.0) / s_per_h

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
        times_h = times / s_per_h - time_shift_h
        mask = times_h >= 0 if time_shift_h > 0 else np.ones(len(times_h), dtype=bool)
        style = _OVERLAY_STYLE.get(case_id, dict(color="black", ls="--", lw=2.5))
        for sp, fracs in data.items():
            y = np.log10(np.where(fracs[mask] > 0, fracs[mask], np.nan)) if is_log else fracs[mask]
            ax.plot(times_h[mask], y, **style, label=f"{case_id} {sp}", zorder=2)

    ax.set_xlim(ref["xlim"])
    ax.set_ylim(plot_ylim)
    if "xtick_step" in ref:
        ax.set_xticks(np.arange(ref["xlim"][0], ref["xlim"][1] + 1, ref["xtick_step"]))
    if is_log:
        yticks = np.arange(np.ceil(plot_ylim[0]), np.floor(plot_ylim[1]) + 1)
        ax.set_yticks(yticks)
        ax.set_yticklabels([f"$10^{{{int(t)}}}$" for t in yticks])
    elif "ytick_step" in ref:
        ax.set_yticks(np.arange(ref["ylim"][0], ref["ylim"][1] + ref["ytick_step"] * 0.1,
                                ref["ytick_step"]))
    ax.set_xlabel(ref.get("xlabel", "Time [h]"))
    ax.set_ylabel(ref.get("ylabel", "Release fraction [−]"))
    case_ids_str = " & ".join(sorted(case_data.keys()))
    ax.set_title(f"IAEA FP diffusion benchmark — {case_ids_str} (overlay)")
    ax.legend(fontsize=8, loc="lower right")

    os.makedirs("figures", exist_ok=True)
    stem = os.path.splitext(os.path.basename(ref["file"]))[0]
    out = os.path.join("figures", f"{stem}_overlay.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"  Overlay → {out}")


if __name__ == "__main__":
    case_ids = sys.argv[1:] if len(sys.argv) > 1 else list(_CASES.keys())

    results: dict[str, tuple[np.ndarray, dict]] = {}
    for cid in case_ids:
        result = post(cid)
        if result is not None:
            results[cid] = result

    for ref_key, ref in REFERENCE_FIGS.items():
        available = ref["cases"] & set(results.keys())
        if not available:
            continue
        subset = {cid: results[cid] for cid in available}
        post_overlay(subset, ref)

    plt.show()
