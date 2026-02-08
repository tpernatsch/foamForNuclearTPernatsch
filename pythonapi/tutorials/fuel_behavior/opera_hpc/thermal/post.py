# study.py
import matplotlib.pyplot as plt
import numpy as np
import os

from case import build_case

# make output folder
FIGDIR = "figures"
os.makedirs(FIGDIR, exist_ok=True)

if __name__ == "__main__":

    # store data: (P, h) -> (r_mm, T(r) at last time)
    radial_profiles = {}
            
    case, rod_mesh = build_case()

    # get radialProfile function
    func_by_name = {f.name: f for f in case.functions}
    radialProfile = func_by_name["radialProfile"]

    lastTime = case.get_time_steps()[-1]
    data, locations = radialProfile.read_from_case(
        startTime=0)

    T_radial = np.array(data["T"][lastTime])
    r_mm     = np.array([r * 1000.0 for r in locations])

    # -------------------------------------------------
    # 1) Radial profile
    # -------------------------------------------------
    fig_all, ax_all = plt.subplots(figsize=(6, 4), dpi=200)

    ax_all.plot(r_mm, T_radial)

    ax_all.set_xlabel("Radius [mm]")
    ax_all.set_ylabel("T [K]")
    ax_all.set_title("Radial temperature profile")
    fig_all.tight_layout()
    fig_all.savefig(f"{FIGDIR}/radial_profile.png")