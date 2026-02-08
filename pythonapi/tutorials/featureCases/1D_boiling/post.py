# post.py
from case import build_case
import os

# Make output folder
FIGDIR = "figures"
os.makedirs(FIGDIR, exist_ok=True)

if __name__ == "__main__":
    case, th_mesh = build_case(region="fluidRegion")
    settings = case.settings

    case.plot_slice(region=th_mesh.region, time=settings.endTime, fieldName="T.liquid", unit="K")
    case.plot_slice(region=th_mesh.region, time=settings.endTime, fieldName="T.vapour", unit="K")
    case.plot_slice(region=th_mesh.region, time=settings.endTime, fieldName="alpha.liquid", unit="-")
    case.plot_slice(region=th_mesh.region, time=settings.endTime, fieldName="alpha.vapour", unit="-")

    case.plot_animation(region=th_mesh, fieldName="alpha.liquid", normal="y", unit="-")
    case.plot_animation(region=th_mesh, fieldName="T.liquid", normal="y", unit="K")
    case.plot_animation(region=th_mesh, fieldName="T.vapour", normal="y", unit="K")

    case.plot_residuals(parameters=["h.liquid", "h.vapour", "p_rgh"], title="Thermal-hydraulics")
