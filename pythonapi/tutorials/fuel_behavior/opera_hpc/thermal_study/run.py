# run.py
from case import build_case

# Parameters defined here
powers = [15e3, 40e3]         # W/m
h_gaps = [5e3, 25e3]          # W/m2K

# Run function
def run_case(
    lhgr: float,
    h_gap: float,
    case_name: str,
):
    """Build case with given (lhgr, h_gap, case_name), export & run it."""
    case, rod_mesh = build_case(lhgr=lhgr, h_gap=h_gap, case_name=case_name)

    print(f"Running case: P={lhgr/1e3:.0f} kW/m, hGap={h_gap:.2e}")
    case.clean()
    case.run()

    return case, rod_mesh


if __name__ == "__main__":
    
    for P in powers:
        for h in h_gaps:
            case_name = f"cases/P{int(P/1e3)}kWpm_h{int(h/1e3)}e3"

            case, rod_mesh = run_case(
                lhgr=P,
                h_gap=h,
                case_name=case_name,
            )