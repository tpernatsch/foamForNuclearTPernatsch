import subprocess
import sys
import os
import glob

_GROUPS = ["cases_1_2_3", "cases_4", "cases_5_6_7", "cases_8"]
_HERE = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    for group in _GROUPS:
        subprocess.run(
            [sys.executable, "clean.py", *sys.argv[1:]],
            cwd=os.path.join(_HERE, group),
            check=True,
        )

    # Remove generated XLSXs from the central reference_figures/ directory
    for xlsx in glob.glob(os.path.join(_HERE, "reference_figures", "*.xlsx")):
        os.remove(xlsx)
        print(f"Removed {xlsx}")
