# post.py
from case import build_case
import os

import matplotlib.pyplot as plt

# Make output folder
FIGDIR = "figures"
os.makedirs(FIGDIR, exist_ok=True)


if __name__ == "__main__":
    case, fluid_mesh, solid_mesh = build_case()

    # ------------------------------------------------------------------
    # Function objects: build lookup dict and extract the one we need
    # ------------------------------------------------------------------
    func_by_name = {f.name: f for f in case.functions}

    probesFunctionObject = func_by_name["probes"]   # change name if needed

    # ------------------------------------------------------------------
    # Mesh plots
    # ------------------------------------------------------------------
    case.plot_mesh(region=[fluid_mesh, solid_mesh], show_edges=True, normal="z", fig_dir=FIGDIR)
    case.plot_mesh(region=fluid_mesh, show_edges=True, normal="z", fig_dir=FIGDIR)
    case.plot_mesh(region=solid_mesh, show_edges=True, normal="z", fig_dir=FIGDIR)

    # ------------------------------------------------------------------
    # Residuals
    # ------------------------------------------------------------------
    case.plot_residuals(
        parameters=["h", "p_fluid"],
        title="Thermal-hydraulics",
        fig_dir=FIGDIR
    )

    # ------------------------------------------------------------------
    # Fluid fields
    # ------------------------------------------------------------------
    for fieldName, unit in [("T", "K"), ("U", "m/s")]:
        case.plot_slice(
            region=fluid_mesh,
            time=case.settings.endTime,
            fieldName=fieldName,
            cmap="RdBu_r",
            unit=unit,
            normal="z",
            fig_dir=FIGDIR
        )

        case.plot_animation(
            region=fluid_mesh,
            fieldName=fieldName,
            cmap="RdBu_r",
            unit=unit,
            fps=2,
            normal="z",
            fig_dir=FIGDIR
        )

    # ------------------------------------------------------------------
    # Combined temperature slice
    # ------------------------------------------------------------------
    case.plot_slice(
        region=[solid_mesh, fluid_mesh],
        time=case.settings.endTime,
        fieldName="T",
        cmap="RdBu_r",
        unit="K",
        normal="z",
        fig_dir=FIGDIR
    )

    case.plot_animation(
        region=solid_mesh,
        fieldName="T",
        cmap="RdBu_r",
        unit="K",
        fps=2,
        normal="z",
        fig_dir=FIGDIR
    )

    # ------------------------------------------------------------------
    # Probe data
    # ------------------------------------------------------------------
    data, locations = probesFunctionObject.read_from_case(
        startTime=0,
        fieldName="T"
    )

    xRange = [loc.x for loc in locations]
    lastTime = case.settings.endTime

    # Normalize temperatures
    T_norm = (data[lastTime] - 300) / 10

    # Plot comparison
    plt.figure(figsize=(8, 4))
    plt.plot(xRange, T_norm, marker="o", label="foamForNuclear")
    # plt.plot(numerical["x"], numerical["y"], label="Expected Numerical", linestyle="--")
    # plt.plot(analytical["x"], analytical["y"], label="Expected Analytical", linestyle="-.")
    plt.xlabel("x [m]")
    plt.ylabel(r"Normalized Temperature $\frac{T - 300}{T_s - 300}$")
    plt.title(f"Normalized Temperature at t = {lastTime}s")
    plt.xlim((0, 1))
    plt.grid(True)
    plt.tight_layout()
    plt.legend()
    plt.savefig(os.path.join(FIGDIR, "fig_results_temperatureComparison.png"), bbox_inches="tight")
    plt.close()