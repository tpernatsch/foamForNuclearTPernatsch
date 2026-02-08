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

    fuel_centerline = func_by_name["fuelCenterline"]
    radial_profile  = func_by_name["radialProfile"]

    # Extract data from probe
    temp, points = fuel_centerline.read_from_case(
        startTime=0,
        fieldName="T"
    )
    
    bu, points = fuel_centerline.read_from_case(
        startTime=0,
        fieldName="Bu"
    )

    time = [t for t in list(temp.keys())]
    temp_list = [val[0] for val in temp.values()]
    bu_list = [val[0]/1000 for val in bu.values()]

    fig, ax_temp = plt.subplots(figsize=(5, 4), dpi=200)

    # ---- Temperature (left axis) ----
    ax_temp.plot(time, temp_list, label="T", color='tab:red')
    ax_temp.tick_params(axis='y', labelcolor='tab:red')
    ax_temp.set_ylabel('Temperature [K]', color='tab:red')

    # ---- Burnup (right axis) ----
    ax_bu = ax_temp.twinx()
    ax_bu.plot(time, bu_list, label="Bu", color='tab:blue')
    ax_bu.tick_params(axis='y', labelcolor='tab:blue')
    ax_bu.set_ylabel('Burnup [GWd/tU]', color='tab:blue')

    ax_temp.set_xlabel('Time, days')

    fig.legend()

    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig_center.png")

    #====== Plot temperature radial profiles

    # Extract data radial profile
    data, locations = radial_profile.read_from_case(
        startTime=0
    )

    radial_T_7 = data["T"][7.04]
    radial_T_168 = data["T"][168.04]
    radial_T_364 = data["T"][364]
    radial_locations = [r*1000 for r in locations]

    fig, ax_radial = plt.subplots(figsize=(5, 4), dpi=200)

    ax_radial.plot(radial_locations, radial_T_7, label= "T @7.04days", color='black', linestyle='-')
    ax_radial.plot(radial_locations, radial_T_168, label= "T @168.04days", color='black', linestyle='--')
    ax_radial.plot(radial_locations, radial_T_364, label= "T @364days", color='black', linestyle='-.')
    
    ax_radial.tick_params(axis='y', labelcolor='black')
    ax_radial.set_xlabel('Radial locations, mm')
    ax_radial.set_ylabel('Temprature, K', color='black')

    fig.legend()

    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig_radial.png")

    # #========

