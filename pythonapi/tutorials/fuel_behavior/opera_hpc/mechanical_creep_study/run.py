# run.py
from case import build_case

BASE = dict(
    case_name="cases/base",
    # default input variables
    clad_ri=4.56e-3,
    clad_thickness=0.7e-3,
    fast_flux=1e13,
    T0=320 + 273,
)

CASES = [
    # case1: change geometry together
    dict(tag="geom", clad_ri=4.91e-3, clad_thickness=0.6e-3),

    # case2: change fast flux
    dict(tag="flux", fast_flux=2.0e13),

    # case3: change temperature
    dict(tag="temp", T0=350 + 273),
]


def run_case(case_name: str, **kwargs):
    params = BASE.copy()
    params.update(kwargs)
    params["case_name"] = case_name

    case, rod_mesh = build_case(**params)

    print(f"Running case: {case_name}")
    case.clean()
    case.run()

    return case, rod_mesh


if __name__ == "__main__":
    # ---- run base case ----
    run_case(case_name=BASE["case_name"])

    # ---- run derived cases ----
    for c in CASES:
        tag = c.pop("tag")
        case_name = f"cases/{tag}"
        run_case(case_name=case_name, **c)
