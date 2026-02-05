# run.py
from case import build_case

if __name__ == "__main__":
    case, rod_mesh = build_case()
    case.clean()
    case.run()
