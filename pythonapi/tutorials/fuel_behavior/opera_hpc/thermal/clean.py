import foamForNuclear as ffn
import os
import shutil

if __name__ == "__main__":    
    ffn.allclean()

    # Clean all parametric case
    if os.path.isdir("case"):
        shutil.rmtree("case")

    # Clean all figures
    if os.path.isdir("figures"):
        shutil.rmtree("figures")