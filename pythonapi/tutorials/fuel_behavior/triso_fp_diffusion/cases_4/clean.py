import sys
import os
import shutil
from _config import CASES

if __name__ == "__main__":
    case_ids = sys.argv[1:] if len(sys.argv) > 1 else list(CASES.keys())
    for case_id in case_ids:
        if os.path.isdir(case_id):
            shutil.rmtree(case_id)

    if os.path.isdir("figures"):
        shutil.rmtree("figures")
