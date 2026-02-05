# postprocess.py
import matplotlib.pyplot as plt
import numpy as np
import os

from case import build_case  # reuse geometry + names

# Make output folder
FIGDIR = "figures"
os.makedirs(FIGDIR, exist_ok=True)


if __name__ == "__main__":

    case, rod_mesh = build_case()

    # ======================
    # Time evolution: U235, Pu239, Pu241 (vs time)
    # ======================
    fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
    
    func_by_name = {f.name: f for f in case.functions}
    sumFields = func_by_name["sum"]

    sum_data, sum_times = sumFields.read_from_case(startTime="0")

    sum_times = np.asarray(sum_times, dtype=float)

    n_u_over_time = sum_data["N_U"]
    n_pu238_over_time = sum_data["N_Pu238"]
    n_pu239_over_time = sum_data["N_Pu239"]
    n_pu240_over_time = sum_data["N_Pu240"]
    n_pu241_over_time = sum_data["N_Pu241"]
    n_pu242_over_time = sum_data["N_Pu242"]

    n_HM_over_time = [
        sum(x)
        for x in zip(
            n_pu238_over_time, n_pu239_over_time, n_pu240_over_time,
            n_pu241_over_time, n_pu242_over_time, n_u_over_time
        )
    ]

    n_u235_normalized = [n / n_HM_over_time[0] for n in sum_data["N_U235"]]
    n_pu239_normalized = [n / n_HM_over_time[0] for n in n_pu239_over_time]
    n_pu241_normalized = [n / n_HM_over_time[0] for n in n_pu241_over_time]

    ax.plot(sum_times, n_u235_normalized, label=f"U_235", color="tab:blue")
    ax.plot(sum_times, n_pu239_normalized, label=f"Pu_239", color="tab:red")
    ax.plot(sum_times, n_pu241_normalized, label=f"Pu_241", color="tab:green")

    ax.set_xlabel("Time")
    ax.set_ylabel("At./At.HM (normalized)")
    ax.legend(loc='upper right')
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig_time_U235_Pu239_Pu241.png", bbox_inches="tight")
    plt.close(fig)

    # ======================
    # Time evolution: PuTot, Pu239, Pu240 (vs time)
    # ======================
    fig, ax = plt.subplots(figsize=(5, 4), dpi=200)

    func_by_name = {f.name: f for f in case.functions}
    sumFields = func_by_name["sum"]

    sum_data, sum_times = sumFields.read_from_case(startTime="0")

    sum_times = np.asarray(sum_times, dtype=float)

    n_u_over_time = sum_data["N_U"]
    n_pu238_over_time = sum_data["N_Pu238"]
    n_pu239_over_time = sum_data["N_Pu239"]
    n_pu240_over_time = sum_data["N_Pu240"]
    n_pu241_over_time = sum_data["N_Pu241"]
    n_pu242_over_time = sum_data["N_Pu242"]

    n_pu_over_time = [
        sum(x)
        for x in zip(
            n_pu238_over_time, n_pu239_over_time, n_pu240_over_time,
            n_pu241_over_time, n_pu242_over_time
        )
    ]
    n_HM_over_time = [
        sum(x)
        for x in zip(
            n_pu238_over_time, n_pu239_over_time, n_pu240_over_time,
            n_pu241_over_time, n_pu242_over_time, n_u_over_time
        )
    ]

    n_pu_normalized = [n / n_HM_over_time[0] for n in n_pu_over_time]
    n_pu239_normalized = [n / n_HM_over_time[0] for n in n_pu239_over_time]
    n_pu240_normalized = [n / n_HM_over_time[0] for n in n_pu240_over_time]

    ax.plot(sum_times, n_pu_normalized, label="Pu_tot", color="tab:orange")
    ax.plot(sum_times, n_pu239_normalized, label="Pu_239", color="tab:red")
    ax.plot(sum_times, n_pu240_normalized, label="Pu_240", color="tab:green")

    ax.set_xlabel("Time")
    ax.set_ylabel("At./At.HM (normalized)")
    ax.legend(loc='upper left')
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig_time_PuTot_Pu239_Pu240.png", bbox_inches="tight")
    plt.close(fig)

    # ======================
    # Radial profiles at lastTime: Pu239
    # ======================
    fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
    
    func_by_name = {f.name: f for f in case.functions}
    radialProfile = func_by_name["radialProfile"]
    volAverage = func_by_name["volAverage"]

    lastTime = case.get_time_steps()[-1]

    volAverageData, volAverageTimes = volAverage.read_from_case(startTime="0")
    avg_Pu_239 = volAverageData["N_Pu239"]

    data, locations = radialProfile.read_from_case(startTime=0)

    radial_Pu_239 = [x / avg_Pu_239[-1] for x in data["N_Pu239"][lastTime]]
    radial_locations = [r * 1000 for r in locations]

    ax.plot(
        radial_locations,
        radial_Pu_239,
        color="tab:red",
    )

    ax.set_xlabel("Radial locations, mm")
    ax.set_ylabel("Normalized profile")
    ax.set_title("Pu-239 radial profile (normalized)")
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig_radial_Pu239.png", bbox_inches="tight")
    plt.close(fig)

    # ======================
    # Radial profiles at lastTime: U235
    # ======================
    fig, ax = plt.subplots(figsize=(5, 4), dpi=200)

    func_by_name = {f.name: f for f in case.functions}
    radialProfile = func_by_name["radialProfile"]
    volAverage = func_by_name["volAverage"]

    lastTime = case.get_time_steps()[-1]

    volAverageData, volAverageTimes = volAverage.read_from_case(startTime="0")
    avg_U_235 = volAverageData["N_U235"]

    data, locations = radialProfile.read_from_case(startTime=0)

    radial_U_235 = [x / avg_U_235[-1] for x in data["N_U235"][lastTime]]
    radial_locations = [r * 1000 for r in locations]

    ax.plot(
        radial_locations,
        radial_U_235,
        color="tab:blue",
    )

    ax.set_xlabel("Radial locations, mm")
    ax.set_ylabel("Normalized profile")
    ax.set_title("U-235 radial profile (normalized)")
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig_radial_U235.png", bbox_inches="tight")
    plt.close(fig)

    # ======================
    # Radial profiles at lastTime: Bu
    # ======================
    fig, ax = plt.subplots(figsize=(5, 4), dpi=200)

    func_by_name = {f.name: f for f in case.functions}
    radialProfile = func_by_name["radialProfile"]
    volAverage = func_by_name["volAverage"]

    lastTime = case.get_time_steps()[-1]

    volAverageData, volAverageTimes = volAverage.read_from_case(startTime="0")
    avg_bu = np.array(volAverageData["Bu"])

    data, locations = radialProfile.read_from_case(startTime=0)

    radial_Bu = [x / avg_bu[-1] for x in data["Bu"][lastTime]]
    radial_locations = [r * 1000 for r in locations]

    ax.plot(
        radial_locations,
        radial_Bu,
        color="tab:red",
    )

    ax.set_xlabel("Radial locations, mm")
    ax.set_ylabel("Normalized profile")
    ax.set_title("Burnup radial profile (normalized)")
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig_radial_Bu.png", bbox_inches="tight")
    plt.close(fig)
