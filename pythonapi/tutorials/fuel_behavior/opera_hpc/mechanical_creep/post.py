# postprocess.py
import matplotlib.pyplot as plt
import numpy as np
import os

from case import build_case  # reuse geometry + names

# Make output folder
FIGDIR = "figures"
os.makedirs(FIGDIR, exist_ok=True)


if __name__ == "__main__":
    # Rebuild a "light" case just to get paths/names
    case, rod_mesh = build_case()

    # Build a simple lookup dict: name -> function object
    func_by_name = {f.name: f for f in case.functions}

    clad_outer = func_by_name["cladOuter"]
    vol_average = func_by_name["volAverage"]

    # ======================
    # Extract data from probe
    # ======================
    disp, points = clad_outer.read_from_case(
        startTime=0,
        fieldName="D"
    )

    times = list(disp.keys())
    disp_list = [np.asarray(v, dtype=float)[0] * 1e3 for v in disp.values()]  # D_x in mm

    fig, ax_disp = plt.subplots(figsize=(5, 4), dpi=200)

    ax_disp.plot(times, disp_list, label="D_x", color="tab:red")
    ax_disp.tick_params(axis="y",)
    ax_disp.set_ylabel("Displacement [mm]",)
    ax_disp.set_ylim(0,0.04)
    ax_disp.set_xlabel("Time [hours]")

    ax_disp.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=3,
    )
    fig.subplots_adjust(bottom=0.25)
    fig.savefig(f"{FIGDIR}/fig_clad_outer.png", bbox_inches="tight")
    plt.close(fig)

    # ======================
    # Extract volume-averaged stress
    # ======================
    data, times = vol_average.read_from_case(startTime=0)

    times = np.asarray(times, dtype=float)
    sigma = np.asarray(data["sigma"], dtype=float)  # expected (nTimes, nComp)

    avg_sigma_hoop = sigma[:, 3] / 1e6

    fig, ax_sigma = plt.subplots(figsize=(5, 4), dpi=200)

    ax_sigma.plot(times, avg_sigma_hoop, label="Hoop Stress", color="tab:red")
    ax_sigma.set_xlabel("Time [hours]")
    ax_sigma.set_ylabel("Stress [MPa]")
    ax_sigma.set_ylim(60,90)
    ax_sigma.legend()

    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig_avg_stress.png", bbox_inches="tight")
    plt.close(fig)
