import runpy
import pathlib

def main():
    # Location of this file: pythonAPI/foamForNuclear/tools/rodmaker_entry.py
    pkg_tools = pathlib.Path(__file__).resolve().parent

    # pythonAPI/
    python_api_root = pkg_tools.parents[1]

    # repository root = one level above pythonAPI
    repo_root = python_api_root.parent

    # real script location: repo_root/tools/rodMaker/rodMaker.py
    script_path = repo_root / "tools" / "rodMaker" / "rodMaker.py"

    if not script_path.exists():
        raise FileNotFoundError(f"Could not find rodMaker.py at {script_path}")

    # Execute the external script as __main__
    runpy.run_path(str(script_path), run_name="__main__")