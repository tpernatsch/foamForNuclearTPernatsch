# clean.py
import foamForNuclear as ffn
import os
import shutil

if __name__ == "__main__":
    ffn.allclean()

    # Clean cases
    if os.path.isdir("cases"):
        shutil.rmtree("cases")

    # Clean figures
    if os.path.isdir("figures"):
        shutil.rmtree("figures")