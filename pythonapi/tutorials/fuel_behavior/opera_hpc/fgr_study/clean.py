import foamForNuclear as ffn
import os
import shutil

if __name__ == "__main__":    
    ffn.allclean()

    # Clean all parametric cases
    if os.path.isdir("cases"):
        shutil.rmtree("cases")

    # Clean all figures
    if os.path.isdir("figures"):
        shutil.rmtree("figures")