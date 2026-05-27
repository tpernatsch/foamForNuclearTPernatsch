"""
Post-processing for IAEA HTGR FP diffusion benchmark case 5 (Cs + Ag,
cyclic irradiation then optional post-irradiation heating).

Usage
-----
    python post.py                    # run (case_5b only)
    python post.py case_5b            # explicit

Two reference-figure overlays are produced when both cases are available:
  1. Full timeline 0–25 000 h (irradiation + heating).
  2. Heating phase only 0–200 h (time-shifted, case_5b only).
"""

import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from _config import CASES as _CASES, T_IRRAD_END, T_HEATING

HERE   = os.path.dirname(os.path.abspath(__file__))
REFDIR = os.path.join(HERE, "..", "reference_figures")

s_per_h = 3600.0

# ── Reference figure definitions ─────────────────────────────────────────────
REFERENCE_FIGS = {
    # Full timeline (irradiation 0–24 000 h + heating)
    "full": dict(
        file=os.path.join(REFDIR, "case_5_fractional_Cs_0h_24000h_1e-8_1_log.png"),
        cases=frozenset({"case_5b"}),
        xlim=(0, 24000),
        ylim=(1e-8, 1),
        xtick_step=5000,
        xlabel="Time  [h]",
        ylabel="Cs / Ag fractional release  [−]",
        yscale="log",
        time_shift=0.0,
    ),
    # Heating phase only: x-axis shifted so 0 = start of heating
    "heating": dict(
        file=os.path.join(REFDIR, "case_5_fractional_Cs_0h_200h_1e-8_1_log.png"),
        cases=frozenset({"case_5b"}),
        xlim=(0, 200),
        ylim=(1e-8, 1),
        xtick_step=50,
        xlabel="Time after start of heating  [h]",
        ylabel="Cs / Ag fractional release  [−]",
        yscale="log",
        time_shift=T_IRRAD_END,
    ),
}

_CASE_STYLES = {
    "case_5b": dict(color="tab:blue", ls="-", lw=2.0),
}

_SP_COLORS = {
    "Cs": "tab:blue",
    "Ag": "tab:red",
}


def _find_dat(case_id: str) -> str:
    files = glob.glob(os.path.join(
        case_id, "postProcessing", "fpRelease", "*", "writeFpRelease.dat"
    ))
    if not files:
        raise FileNotFoundError(f"No fpRelease dat under {case_id}/postProcessing/")
    files.sort(key=lambda p: float(os.path.basename(os.path.dirname(p))))
    return files[0]


def _read_dat(fpath: str) -> tuple[np.ndarray, dict[str, np.ndarray]]:
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


def post(case_id: str = "case_5b") -> tuple[np.ndarray, dict[str, np.ndarray]] | None:
    """Plot case_5b: irradiation panel (0–24 000 h) + heating panel (0–200 h).

    The case_5a snapshot (end-of-irradiation release) is read from the same
    simulation output and reported alongside the end-of-heating result.
    """
    cfg = _CASES[case_id]
    expected_5b = cfg.get("expected", {})
    expected_5a = cfg.get("expected_5a", {})

    try:
        fpath = _find_dat(case_id)
    except FileNotFoundError as e:
        print(f"  [{case_id}] skipped: {e}")
        return None

    times, data = _read_dat(fpath)
    times_h = times / s_per_h
    t_irrad_h = T_IRRAD_END / s_per_h                           # 24 000 h
    heat_end_h = (cfg["endTime"] - T_IRRAD_END) / s_per_h       # 200 h

    species_list = list(data.keys())
    n = len(species_list)

    print(f"\n{'─'*60}")
    print(f"  {case_id}  (irradiation → heating @ {T_HEATING:.0f} K)")
    print(f"{'─'*60}")

    # Report end-of-irradiation values (case_5a snapshot)
    print(f"  At end of irradiation (t = {T_IRRAD_END/s_per_h:.0f} h):")
    print(f"  {'Species':>8}  {'Sim value':>12}  {'Expected':>12}  {'Error %':>9}")
    for sp in species_list:
        fracs = data[sp]
        # find closest time index to T_IRRAD_END
        idx_5a = int(np.argmin(np.abs(times - T_IRRAD_END)))
        val_5a = fracs[idx_5a]
        exp_5a = expected_5a.get(sp)
        if exp_5a is not None:
            err = (val_5a - exp_5a) / exp_5a * 100 if exp_5a != 0 else float("nan")
            print(f"  {sp:>8}  {val_5a:>12.6g}  {exp_5a:>12.6g}  {err:>+9.2f}")
        else:
            print(f"  {sp:>8}  {val_5a:>12.6g}  {'N/A':>12}")

    # Report end-of-heating values (case_5b result)
    print(f"\n  At end of heating (t = {cfg['endTime']/s_per_h:.0f} h):")
    print(f"  {'Species':>8}  {'Sim value':>12}  {'Expected':>12}  {'Error %':>9}")
    for sp in species_list:
        fracs = data[sp]
        final = fracs[-1]
        exp_5b = expected_5b.get(sp)
        if exp_5b is not None:
            err = (final - exp_5b) / exp_5b * 100 if exp_5b != 0 else float("nan")
            print(f"  {sp:>8}  {final:>12.6g}  {exp_5b:>12.6g}  {err:>+9.2f}")
        else:
            print(f"  {sp:>8}  {final:>12.6g}  {'N/A':>12}")

    fig, axes = plt.subplots(n, 2, figsize=(10, 4 * n), squeeze=False)
    fig.suptitle(
        f"FP release — case 5  (10 × 100 efpd cyclic + {heat_end_h:.0f} h @ {T_HEATING:.0f} K)",
        fontsize=11,
    )

    for i, sp in enumerate(species_list):
        fracs = data[sp]
        color = _SP_COLORS.get(sp, "tab:blue")
        exp_5b = expected_5b.get(sp)
        exp_5a = expected_5a.get(sp)

        # Irradiation panel (0 → 24 000 h)
        ax_irr = axes[i][0]
        mask_irr = times_h <= t_irrad_h
        ax_irr.plot(times_h[mask_irr], fracs[mask_irr], color=color, lw=1.5, label="Offbeat")
        if exp_5a is not None:
            ax_irr.axhline(exp_5a, color="black", lw=1.2, ls="--",
                           label=f"IAEA ref (5a) = {exp_5a:.4g}")
        ax_irr.set_xlim(0, t_irrad_h)
        ax_irr.set_yscale("log")
        ax_irr.set_ylim(1e-8, 1.0)
        ax_irr.set_xlabel("Time [h]")
        ax_irr.set_ylabel(f"{sp} release fraction [−]")
        ax_irr.set_title(f"{sp} — irradiation (case 5a snapshot)")
        ax_irr.legend(fontsize=8)
        ax_irr.grid(True, alpha=0.3, which="both")

        # Heating panel (0 → 200 h, shifted)
        ax_heat = axes[i][1]
        mask_heat = times_h >= t_irrad_h
        ax_heat.plot(times_h[mask_heat] - t_irrad_h, fracs[mask_heat],
                     color=color, lw=1.5, label="Offbeat")
        if exp_5b is not None:
            ax_heat.axhline(exp_5b, color="black", lw=1.2, ls="--",
                            label=f"IAEA ref (5b) = {exp_5b:.4g}")
        ax_heat.set_xlim(0, heat_end_h)
        ax_heat.set_yscale("log")
        ax_heat.set_ylim(1e-8, 1.0)
        ax_heat.set_xlabel("Time after start of heating [h]")
        ax_heat.set_ylabel(f"{sp} release fraction [−]")
        ax_heat.set_title(f"{sp} — heating (case 5b)")
        ax_heat.legend(fontsize=8)
        ax_heat.grid(True, alpha=0.3, which="both")

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
        # When time_shift is non-zero, only plot data after the shift point
        mask = times_h >= 0 if time_shift_h > 0 else np.ones(len(times_h), dtype=bool)
        style = _CASE_STYLES.get(case_id, dict(color="black", ls="--", lw=2.0))
        for sp, fracs in data.items():
            sp_color = _SP_COLORS.get(sp, style["color"])
            y = np.log10(np.where(fracs[mask] > 0, fracs[mask], np.nan)) if is_log else fracs[mask]
            ax.plot(times_h[mask], y, color=sp_color, ls=style["ls"], lw=style["lw"],
                    label=f"{case_id} {sp}", zorder=2)

    ax.set_xlim(ref["xlim"])
    ax.set_ylim(plot_ylim)
    if "xtick_step" in ref:
        ax.set_xticks(np.arange(ref["xlim"][0], ref["xlim"][1] + 1, ref["xtick_step"]))
    if is_log:
        yticks = np.arange(np.ceil(plot_ylim[0]), np.floor(plot_ylim[1]) + 1)
        ax.set_yticks(yticks)
        ax.set_yticklabels([f"$10^{{{int(t)}}}$" for t in yticks])
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
    # Only case_5b is run; it covers the full history.
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
