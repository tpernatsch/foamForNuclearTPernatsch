import foamForNuclear as ffn
import shutil
import os

ffn.allclean()

# Clean folders
for folder in ['steadyState', 'transient']:
    if os.path.exists(folder):
        shutil.rmtree(folder)
