# study.py
import matplotlib.pyplot as plt
import numpy as np
import os

from case import build_case
from run import powers, h_gaps

# make output folder
FIGDIR = "figures"
os.makedirs(FIGDIR, exist_ok=True)

if __name__ == "__main__":

    # store data: (P, h) -> (r_mm, T(r) at last time)
    radial_profiles = {}

    for P in powers:
        for h in h_gaps:
            case_name = f"cases/P{int(P/1e3)}kWpm_h{int(h/1e3)}e3"

            case, rod_mesh = build_case(
                lhgr=P,
                h_gap=h,
                case_name=case_name,
            )

            # get radialProfile function
            func_by_name = {f.name: f for f in case.functions}
            radialProfile = func_by_name["radialProfile"]

            lastTime = case.get_time_steps()[-1]
            data, locations = radialProfile.read_from_case(
                startTime=0)

            T_radial = np.array(data["T"][lastTime])
            r_mm     = np.array([r * 1000.0 for r in locations])

            radial_profiles[(P, h)] = (r_mm, T_radial)

    # -------------------------------------------------
    # 1) All radial profiles together
    # -------------------------------------------------
    fig_all, ax_all = plt.subplots(figsize=(6, 4), dpi=200)

    for (P, h), (r_mm, T_prof) in radial_profiles.items():
        label = f"P={P/1e3:.0f} kW/m, h={h/1e3:.0f}·10³"
        ax_all.plot(r_mm, T_prof, label=label)

    ax_all.set_xlabel("Radius [mm]")
    ax_all.set_ylabel("T [K]")
    ax_all.set_title("Radial temperature profiles – all cases")
    ax_all.legend(fontsize=7)
    fig_all.tight_layout()
    fig_all.savefig(f"{FIGDIR}/radial_all_cases.png")

    # -------------------------------------------------
    # 2) For a given power: different hGap
    # -------------------------------------------------
    for P in powers:
        fig_h, ax_h = plt.subplots(figsize=(6, 4), dpi=200)
        for h in h_gaps:
            r_mm, T_prof = radial_profiles[(P, h)]
            label = f"hGap={h/1e3:.0f}·10³ W/m²K"
            ax_h.plot(r_mm, T_prof, label=label)

        ax_h.set_xlabel("Radius [mm]")
        ax_h.set_ylabel("T [K]")
        ax_h.set_title(f"Radial T profiles – P={P/1e3:.0f} kW/m")
        ax_h.legend(fontsize=8)
        fig_h.tight_layout()
        fig_h.savefig(f"{FIGDIR}/radial_P{int(P/1e3)}kWpm.png")
