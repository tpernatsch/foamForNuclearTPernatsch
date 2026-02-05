# study.py
import matplotlib.pyplot as plt
import numpy as np
import os

from case import build_case
from run import r_grains, powers  # vectors defined in run.py

# make output folder
FIGDIR = "figures"
os.makedirs(FIGDIR, exist_ok=True)

if __name__ == "__main__":

    # store data: (r_grain, T0) -> (time, fgr)
    fgr_series = {}

    for P in powers:
        for r_g in r_grains:
            case_name = f"cases/P{int(P/1e3)}kWpm_grain{r_g}m"

            case, rod_mesh = build_case(
                r_grain=r_g,
                lhgr=P,
                case_name=case_name,
            )

            # get fgr function
            func_by_name = {f.name: f for f in case.functions}
            fgr = func_by_name["fgr"]   # adapt if your FO name differs

            data = fgr.read_from_case(startTime=0)

            time = np.array(list(data.keys()), dtype=float)
            fgr_vals = np.array(list(data.values()), dtype=float)

            fgr_series[(r_g, P)] = (time, fgr_vals)

    # -------------------------------------------------
    # 1) All fgr(t) together
    # -------------------------------------------------
    fig_all, ax_all = plt.subplots(figsize=(6, 4), dpi=200)

    for (r_g, P), (time, fgr_vals) in fgr_series.items():
        label = f"r_g={r_g}, P={P}"
        ax_all.plot(time, fgr_vals, label=label)

    ax_all.set_xlabel("Time")
    ax_all.set_ylabel("fgr")
    ax_all.set_title("FGR – all cases")
    ax_all.legend(fontsize=7)
    fig_all.tight_layout()
    fig_all.savefig(f"{FIGDIR}/fgr_all_cases.png")

    # -------------------------------------------------
    # 2) For a given r_grain: different T0
    # -------------------------------------------------
    for r_g in r_grains:
        fig_t, ax_t = plt.subplots(figsize=(6, 4), dpi=200)

        for P in powers:
            time, fgr_vals = fgr_series[(r_g, P)]
            label = f"P={P}"
            ax_t.plot(time, fgr_vals, label=label)

        ax_t.set_xlabel("Time")
        ax_t.set_ylabel("fgr")
        ax_t.set_title(f"FGR – r_grain={r_g}")
        ax_t.legend(fontsize=8)
        fig_t.tight_layout()
        fig_t.savefig(f"{FIGDIR}/fgr_rg{r_g}.png")
