import foamForNuclear as ffn
import os
import shutil

ffn.allclean()

# Clean all parametric cases
if os.path.isdir("cases"):
    shutil.rmtree("cases")