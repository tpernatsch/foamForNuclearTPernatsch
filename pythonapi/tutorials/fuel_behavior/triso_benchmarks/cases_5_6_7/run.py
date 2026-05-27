import sys
from case import build_case
from _config import CASES

if __name__ == "__main__":
    case_ids = sys.argv[1:] if len(sys.argv) > 1 else list(CASES.keys())
    for case_id in case_ids:
        case, _ = build_case(case_id)
        case.clean()
        case.run()
