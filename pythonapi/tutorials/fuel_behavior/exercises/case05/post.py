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
    fuel_outer = func_by_name["fuelOuter"]
    radial_profile  = func_by_name["radialProfile"]

    # Extract data from probes
    temp, points = fuel_centerline.read_from_case(
        startTime=0,
        fieldName="T"
    )
    
    bu, points = fuel_centerline.read_from_case(
        startTime=0,
        fieldName="Bu"
    )
    
    gw, points = fuel_outer.read_from_case(
        startTime=0,
        fieldName="gapWidth"
    )

    time = [t for t in list(temp.keys())]
    temp_list = [val[0] for val in temp.values()]
    bu_list = [val[0]/1000 for val in bu.values()]
    gw_list = [val[0]*1000 for val in gw.values()]

    fig, ax_temp = plt.subplots(figsize=(5, 4), dpi=200)

    # ---- Temperature (left axis) ----
    l1, = ax_temp.plot(time, temp_list, label="T", color='tab:red')
    ax_temp.tick_params(axis='y', labelcolor='tab:red')
    ax_temp.set_ylabel('Temperature [K]', color='tab:red')

    # ---- Burnup (right axis #1) ----
    ax_bu = ax_temp.twinx()
    l2, = ax_bu.plot(time, bu_list, label="Bu", color='tab:blue')
    ax_bu.tick_params(axis='y', labelcolor='tab:blue')
    ax_bu.set_ylabel('Burnup [GWd/tU]', color='tab:blue')

    # ---- Gap width (right axis #2, offset) ----
    ax_gw = ax_temp.twinx()
    ax_gw.spines['right'].set_position(('axes', 1.2))   # shift outward
    ax_gw.spines['right'].set_visible(True)             # make sure it draws
    l3, = ax_gw.plot(time, gw_list, label="Gap width", color='tab:green')
    ax_gw.tick_params(axis='y', labelcolor='tab:green')
    ax_gw.set_ylabel('Gap width [mm]', color='tab:green')

    # ---- X axis ----
    ax_temp.set_xlabel('Time [days]')

    # Legend: collect handles explicitly from all axes
    lines = [l1, l2, l3]
    labels = [ln.get_label() for ln in lines]
    
    # Legend outside
    ax_temp.legend(
        lines, labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=3,
    )
    fig.subplots_adjust(bottom=0.25)

    # fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig_center.png", bbox_inches="tight")

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
    fig.savefig(f"{FIGDIR}/fig_radial_T.png")

    #====== Plot radial stress radial profiles

    # Extract data radial profile
    data, locations = radial_profile.read_from_case(
        startTime=0
    )

    radial_sigma_r_7 = [sigma[0]/1e6 for sigma in data["sigma"][7.04]]
    radial_sigma_r_168 = [sigma[0]/1e6 for sigma in data["sigma"][168.04]]
    radial_sigma_r_364 = [sigma[0]/1e6 for sigma in data["sigma"][364]]
    radial_locations = [r*1000 for r in locations]

    fig, ax_radial = plt.subplots(figsize=(5, 4), dpi=200)

    ax_radial.plot(radial_locations, radial_sigma_r_7, label= "Radial stress @7.04days", color='red', linestyle='-')
    ax_radial.plot(radial_locations, radial_sigma_r_168, label= "Radial stress @168.04days", color='red', linestyle='--')
    ax_radial.plot(radial_locations, radial_sigma_r_364, label= "Radial stress @364days", color='red', linestyle='-.')
    
    ax_radial.tick_params(axis='y',)
    ax_radial.set_xlabel('Radial locations, mm')
    ax_radial.set_ylabel('Radial stress, MPa',)

    fig.legend()

    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig_radial_stress_radial.png")

    #====== Plot hoop stress radial profiles

    # Extract data radial profile
    data, locations = radial_profile.read_from_case(
        startTime=0
    )

    radial_sigma_hoop_7 = [sigma[3]/1e6 for sigma in data["sigma"][7.04]]
    radial_sigma_hoop_168 = [sigma[3]/1e6 for sigma in data["sigma"][168.04]]
    radial_sigma_hoop_364 = [sigma[3]/1e6 for sigma in data["sigma"][364]]
    radial_locations = [r*1000 for r in locations]

    fig, ax_radial = plt.subplots(figsize=(5, 4), dpi=200)

    ax_radial.plot(radial_locations, radial_sigma_hoop_7, label= "Hoop stress @7.04days", color='blue', linestyle='-')
    ax_radial.plot(radial_locations, radial_sigma_hoop_168, label= "Hoop stress @168.04days", color='blue', linestyle='--')
    ax_radial.plot(radial_locations, radial_sigma_hoop_364, label= "Hoop stress @364days", color='blue', linestyle='-.')
    
    ax_radial.tick_params(axis='y',)
    ax_radial.set_xlabel('Radial locations, mm')
    ax_radial.set_ylabel('Hoop stress, MPa',)

    fig.legend()

    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig_radial_stress_hoop.png")

    # #========

