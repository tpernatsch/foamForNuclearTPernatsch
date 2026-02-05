# clean.py
import foamForNuclear as ffn
import os
import shutil

if __name__ == "__main__":
    ffn.allclean()

    # Clean case
    if os.path.isdir("case"):
        shutil.rmtree("case")

    # Clean figures
    if os.path.isdir("figures"):
        shutil.rmtree("figures")

    # Clean extracted MATLAB data
    if os.path.isfile("extracted_data.txt"):
        os.remove("extracted_data.txt")