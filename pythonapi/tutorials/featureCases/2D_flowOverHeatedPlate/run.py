# run.py
from case import build_case

if __name__ == "__main__":
    case, fluid_mesh, solid_mesh = build_case()
    case.clean()
    case.run()
