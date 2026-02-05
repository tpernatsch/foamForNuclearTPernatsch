#!/usr/bin/env python3
from pathlib import Path

out = Path("extracted_data.txt")


def read_xy(path, col=1):
    """
    Read a multi-column file (CSV or whitespace).

    Parameters
    ----------
    path : str or Path
        File to read
    col : int
        Column index to extract as y (0-based). Default = 1.

    Skips:
      - empty lines
      - lines starting with '#'
      - the first non-comment line (CSV header)
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

            # Split line (CSV or whitespace)
            if "," in line:
                parts = [p.strip() for p in line.split(",")]
            else:
                parts = line.split()

            # Safety check
            if len(parts) <= col:
                raise ValueError(
                    f"{path}: requested column {col} but only {len(parts)} columns found"
                )

            x.append(float(parts[0]))
            y.append(float(parts[col]))

    return x, y


def write_array(f, name, arr):
    f.write(f"{name} = [")
    f.write(" ".join(f"{v:g}" for v in arr))
    f.write("];\n\n")


# -----------------------------------------------------------------
# Input files
# -----------------------------------------------------------------

centerline_file_T = Path("case/postProcessing/fuelCenterline/0/T")
centerline_file_gapW = Path("case/postProcessing/fuelOuter/0/gapWidth")
centerline_file_ip = Path("case/postProcessing/fuelOuter/0/interfaceP")
centerline_file_hGap = Path("case/postProcessing/fuelOuter/0/hGap")

radial_7_04   = Path("case/postProcessing/radialProfile/7.04/line_T_sigma.csv")
radial_168_04 = Path("case/postProcessing/radialProfile/168.04/line_T_sigma.csv")
radial_364    = Path("case/postProcessing/radialProfile/364/line_T_sigma.csv")

# -----------------------------------------------------------------
# Read data
# -----------------------------------------------------------------

time, temp = read_xy(centerline_file_T)
_, gapW       = read_xy(centerline_file_gapW)
_, interfaceP = read_xy(centerline_file_ip)
_, hGap       = read_xy(centerline_file_hGap)
rPos, temp_7_04     = read_xy(radial_7_04)
_,    temp_168_04   = read_xy(radial_168_04)
_,    temp_364      = read_xy(radial_364)
_,    sigmaR_7_04   = read_xy(radial_7_04, 2)
_,    sigmaR_168_04 = read_xy(radial_168_04, 2)
_,    sigmaR_364    = read_xy(radial_364, 2)
_,    sigmaH_7_04   = read_xy(radial_7_04, 5)
_,    sigmaH_168_04 = read_xy(radial_168_04, 5)
_,    sigmaH_364    = read_xy(radial_364, 5)

# -----------------------------------------------------------------
# Write MATLAB-style file
# -----------------------------------------------------------------

with open(out, "w") as f:
    write_array(f, "time", time)
    write_array(f, "temp", temp)
    write_array(f, "gapW", gapW)
    write_array(f, "hGap", hGap)
    write_array(f, "interfaceP", interfaceP)

    write_array(f, "rPos", rPos)
    write_array(f, "radialTemp_7_04",   temp_7_04)
    write_array(f, "radialTemp_168_04", temp_168_04)
    write_array(f, "radialTemp_364",    temp_364)
    write_array(f, "radialSigmaR_7_04",   sigmaR_7_04)
    write_array(f, "radialSigmaR_168_04", sigmaR_168_04)
    write_array(f, "radialSigmaR_364",    sigmaR_364)
    write_array(f, "radialSigmaH_7_04",   sigmaH_7_04)
    write_array(f, "radialSigmaH_168_04", sigmaH_168_04)
    write_array(f, "radialSigmaH_364",    sigmaH_364)

print(f"Written {out}")