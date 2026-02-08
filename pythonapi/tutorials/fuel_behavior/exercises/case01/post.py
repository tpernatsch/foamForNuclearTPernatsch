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

    time = [t for t in list(temp.keys())]
    temp_List = [val[0] for val in temp.values()]

    fig, ax_temp = plt.subplots(figsize=(5, 4), dpi=200)

    color = 'tab:red'
    ax_temp.plot(time, temp_List, label= "T", color=color)
    ax_temp.tick_params(axis='y', labelcolor=color)
    ax_temp.set_xlabel('Time, days')
    ax_temp.set_ylabel('Temperature, K', color=color)

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

