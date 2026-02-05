# run.py
from case import build_case

# Parameters defined here
powers = [20e3, 25e3]         # W/m
r_grains = [5e-6, 50e-6]         # m

# Run function
def run_case(
    lhgr: float,
    r_grain: float,
    case_name: str,
):
    """Build case with given (lhgr, r_g, case_name), export & run it."""
    case, rod_mesh = build_case(lhgr=lhgr, r_grain=r_grain, case_name=case_name)

    print(f"Running case: P={lhgr/1e3:.0f} kW/m, grain radius={r_grain} m")
    case.clean()
    case.run()

    return case, rod_mesh


if __name__ == "__main__":
    
    for P in powers:
        for r_g in r_grains:
            case_name = f"cases/P{int(P/1e3)}kWpm_grain{r_g}m"

            case, rod_mesh = run_case(
                lhgr=P,
                r_grain=r_g,
                case_name=case_name,
            )