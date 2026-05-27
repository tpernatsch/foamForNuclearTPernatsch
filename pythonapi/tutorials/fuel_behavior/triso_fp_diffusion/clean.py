import subprocess
import sys
import os

_GROUPS = ["cases_1_2", "cases_3", "cases_4", "cases_5"]
_HERE = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    for group in _GROUPS:
        subprocess.run(
            [sys.executable, "clean.py", *sys.argv[1:]],
            cwd=os.path.join(_HERE, group),
            check=True,
        )
