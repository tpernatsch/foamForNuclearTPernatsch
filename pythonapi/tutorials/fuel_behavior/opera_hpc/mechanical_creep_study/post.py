# postprocess.py
import os
import numpy as np
import matplotlib.pyplot as plt

from case import build_case  # reuse geometry + names
from run import BASE, CASES   # reuse the same inputs as run.py

FIGDIR = "figures"
os.makedirs(FIGDIR, exist_ok=True)


if __name__ == "__main__":

    # Build the list of cases to postprocess (base + derived)
    cases_to_plot = [dict(tag="base"), *[c.copy() for c in CASES]]

    # ======================
    # Figure 1: Displacement (all cases)
    # ======================
    fig, ax_disp = plt.subplots(figsize=(5, 4), dpi=200)

    for c in cases_to_plot:
        tag = c["tag"]

        params = BASE.copy()
        params.update({k: v for k, v in c.items() if k != "tag"})
        case_name = params["case_name"] if tag == "base" else f"cases/{tag}"
        params["case_name"] = case_name

        case, rod_mesh = build_case(**params)
        func_by_name = {f.name: f for f in case.functions}
        clad_outer = func_by_name["cladOuter"]

        disp, points = clad_outer.read_from_case(startTime=0, fieldName="D")

        times = list(disp.keys())
        disp_list = [np.asarray(v, dtype=float)[0] * 1e3 for v in disp.values()]  # D_x in mm

        ax_disp.plot(times, disp_list, label=tag)

    ax_disp.set_ylabel("Displacement [mm]")
    ax_disp.set_ylim(0,0.04)
    ax_disp.set_xlabel("Time [hours]")
    ax_disp.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=3,
    )
    fig.subplots_adjust(bottom=0.25)
    fig.savefig(f"{FIGDIR}/fig_clad_outer_all_cases.png", bbox_inches="tight")
    plt.close(fig)

    # ======================
    # Figure 2: Stress (all cases)
    # ======================
    fig, ax_sigma = plt.subplots(figsize=(5, 4), dpi=200)

    for c in cases_to_plot:
        tag = c["tag"]

        params = BASE.copy()
        params.update({k: v for k, v in c.items() if k != "tag"})
        case_name = params["case_name"] if tag == "base" else f"cases/{tag}"
        params["case_name"] = case_name

        case, rod_mesh = build_case(**params)
        func_by_name = {f.name: f for f in case.functions}
        vol_average = func_by_name["volAverage"]

        data, times = vol_average.read_from_case(startTime=0)

        times = np.asarray(times, dtype=float)
        sigma = np.asarray(data["sigma"], dtype=float)  # expected (nTimes, nComp)

        avg_sigma_hoop = sigma[:, 3] / 1e6

        ax_sigma.plot(times, avg_sigma_hoop, label=tag)

    ax_sigma.set_xlabel("Time [hours]")
    ax_sigma.set_ylabel("Stress [MPa]")
    ax_sigma.set_ylim(60,90)
    ax_sigma.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=2,
    )
    fig.subplots_adjust(bottom=0.28)
    fig.savefig(f"{FIGDIR}/fig_avg_stress_all_cases.png", bbox_inches="tight")
    plt.close(fig)
