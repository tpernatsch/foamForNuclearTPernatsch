import os
import shutil

if __name__ == "__main__":
    if os.path.isdir("case_8"):
        shutil.rmtree("case_8")
    if os.path.isdir("figures"):
        shutil.rmtree("figures")
