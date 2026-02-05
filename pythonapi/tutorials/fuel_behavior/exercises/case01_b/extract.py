#!/usr/bin/env python3
from pathlib import Path

out = Path("extracted_data.txt")


from pathlib import Path

def read_xy(path):
    """
    Read a two-column file (CSV or whitespace).
    Skips:
      - empty lines
      - lines starting with '#'
      - the first non-comment line ONLY for .csv files (CSV header)
    """
    x, y = [], []

    is_csv = Path(path).suffix.lower() == ".csv"
    first_data_line = True

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            # Skip header only for CSV
            if is_csv and first_data_line:
                first_data_line = False
                continue

            first_data_line = False

            # allow comma or whitespace
            if "," in line:
                a, b = line.split(",")[:2]
            else:
                a, b = line.split()[:2]

            x.append(float(a))
            y.append(float(b))

    return x, y



def write_array(f, name, arr):
    f.write(f"{name} = [")
    f.write(" ".join(f"{v:g}" for v in arr))
    f.write("];\n\n")


# -----------------------------------------------------------------
# Input files
# -----------------------------------------------------------------

centerline_file = Path("case/postProcessing/fuelCenterline/0/T")

radial_7_04   = Path("case/postProcessing/radialProfile/7.04/line_T.csv")
radial_168_04 = Path("case/postProcessing/radialProfile/168.04/line_T.csv")
radial_364    = Path("case/postProcessing/radialProfile/364/line_T.csv")

# -----------------------------------------------------------------
# Read data
# -----------------------------------------------------------------

time, temp = read_xy(centerline_file)
rPos, temp_7_04   = read_xy(radial_7_04)
_,    temp_168_04 = read_xy(radial_168_04)
_,    temp_364    = read_xy(radial_364)

# -----------------------------------------------------------------
# Write MATLAB-style file
# -----------------------------------------------------------------

with open(out, "w") as f:
    write_array(f, "time", time)
    write_array(f, "temp", temp)

    write_array(f, "rPos", rPos)
    write_array(f, "radialTemp_7_04",   temp_7_04)
    write_array(f, "radialTemp_168_04", temp_168_04)
    write_array(f, "radialTemp_364",    temp_364)

print(f"Written {out}")
