"""
From ThermoPower

The perfomance characteristics are specified by two characteristic equations:
the first relates the flow number phic, the pressure ratio PR and the referred
speed N_T; the second relates the efficiency eta, the flow number phic, and the
referred speed N_T [1]. To avoid singularities, the two characteristic equations
are expressed in parametric form by adding a further variable beta (method of
beta lines [2]).

The performance maps are thus tabulated into three differents tables, tablePhic,
tablePR and tableEta, which express phic, PR and eta as a function of N_T and
beta, respectively, where N_T is the first row while beta is the first column.
The referred speed N_T is defined as a percentage of the design referred speed
and beta are arbitrary lines, usually drawn parallel to the surge-line on the
performance maps.

The tables have been extended using ChatGPT and are therefore not to be used in
production. These are only representative results.
"""
#==============================================================================*
# Imports

import matplotlib.pyplot as plt
import numpy as np


#==============================================================================*
# Functions

def convertModelicaTable(rawData):
    data = []
    for line in rawData.split('\n'):
        if (line == ''):
            continue

        line = line.replace(';', '')
        data.append([float(e) for e in line.split(',')])

    return(np.array(data))


def plot(
        charac: np.array,
        dataName: str,
        xlabel: str,
        ylabel: str,
        zlabel: str,
        zscale: float=1
    ):
    fig, ax = plt.subplots(figsize=(5, 4), dpi=200)

    print(f"{dataName} table shape = {charac.shape}")

    xvalues = charac[0][1:]
    yvalues = charac[1:][:,0]
    zvalues = zscale*charac[1:][:,1:]
    CS = ax.contourf(xvalues, yvalues, zvalues, 20)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    cbar = fig.colorbar(CS)
    cbar.ax.set_ylabel(zlabel)

    fig.tight_layout()
    fig.savefig(f"fig_input_{dataName}.png")


#==============================================================================*
# Data

#                                 beta, N_T1, N_T2, N_T3
# parameter Real tableEtaC[6, 4] = [
#   0, 95,      100,      105;
#   1, 82.5e-2, 81e-2,    80.5e-2;
#   2, 84e-2,   82.9e-2,  82e-2;
#   3, 83.2e-2, 82.2e-2,  81.5e-2;
#   4, 82.5e-2, 81.2e-2,  79e-2;
#   5, 79.5e-2, 78e-2,    76.5e-2
# ];
tableEtaC = """
0,    0.1    50,    95,    100,   105,   125,   150;
0.5,  0.60,  0.68,  0.810, 0.79,  0.79,  0.75,  0.70;
1.0,  0.62,  0.70,  0.825, 0.81,  0.805, 0.76,  0.71;
2.0,  0.65,  0.73,  0.84,  0.829, 0.82,  0.78,  0.73;
3.0,  0.64,  0.72,  0.832, 0.822, 0.815, 0.775, 0.72;
4.0,  0.62,  0.70,  0.825, 0.812, 0.79,  0.74,  0.70;
5.0,  0.59,  0.66,  0.795, 0.78,  0.765, 0.72,  0.68;
6.0,  0.55,  0.63,  0.765, 0.75,  0.735, 0.70,  0.65
"""

# parameter Real tablePhicC[6, 4] = [
# 0, 95,          100,         105;
# 1, 38.3e-3*0.3, 43e-3*0.3,   46.8e-3*0.3;
# 2, 39.3e-3*0.3, 43.8e-3*0.3, 47.9e-3*0.3;
# 3, 40.6e-3*0.3, 45.2e-3*0.3, 48.4e-3*0.3;
# 4, 41.6e-3*0.3, 46.1e-3*0.3, 48.9e-3*0.3;
# 5, 42.3e-3*0.3, 46.6e-3*0.3, 49.3e-3*0.3
# ];
tablePhicC = """
0,    0.1,   50,    95,     100,    105,    125,   150;
0.5,  0.009, 0.019, 0.0373, 0.0422, 0.0457, 0.050, 0.053;
1.0,  0.010, 0.020, 0.0383, 0.0430, 0.0468, 0.051, 0.055;
2.0,  0.011, 0.021, 0.0393, 0.0438, 0.0479, 0.052, 0.057;
3.0,  0.012, 0.022, 0.0406, 0.0452, 0.0484, 0.053, 0.058;
4.0,  0.013, 0.023, 0.0416, 0.0461, 0.0489, 0.054, 0.059;
5.0,  0.014, 0.024, 0.0423, 0.0466, 0.0493, 0.055, 0.060;
6.0,  0.015, 0.025, 0.0430, 0.0470, 0.0498, 0.056, 0.061;
7.0,  0.016, 0.026, 0.0435, 0.0474, 0.0500, 0.0565, 0.062
"""

# parameter Real tablePR[6, 4] = [
# 0, 95,    100,  105;
# 1, 22.6,  27,   32;
# 2, 22,    26.6, 30.8;
# 3, 20.8,  25.5, 29;
# 4, 19,    24.3, 27.1;
# 5, 17,    21.5, 24.2
# ];
tablePR = """
0,     0.1,   50,    95,     100,    105,    125,    150;
0.5,   1.1,   4.0,   23.2,   27.4,   33.2,   40.0,   43.0;
1.0,   1.2,   5.0,   22.6,   27.0,   32.0,   38.0,   42.0;
2.0,   1.3,   6.0,   22.0,   26.6,   30.8,   36.0,   41.0;
3.0,   1.5,   7.5,   20.8,   25.5,   29.0,   34.5,   39.5;
4.0,   1.6,   8.5,   19.0,   24.3,   27.1,   32.0,   37.0;
5.0,   1.6,   9.0,   17.0,   21.5,   24.2,   29.0,   34.0;
6.0,   1.5,   8.8,   15.0,   19.5,   22.0,   26.5,   31.0;
7.0,   1.3,   8.0,   13.0,   17.0,   20.0,   24.0,   29.0
"""


#   PR,     N_T1,       N_T2,       N_T3
# parameter Real tableEtaT[5, 4] = [
#   1,      90,         100,        110;
#   2.36,   89e-2,      89.5e-2,    89.3e-2;
#   2.88,   90e-2,      90.6e-2,    90.5e-2;
#   3.56,   90.5e-2,    90.6e-2,    90.5e-2;
#   4.46,   90.2e-2,    90.3e-2,    90e-2
# ];
tableEtaT = """
1,     0.1,   10,    20,    30,    40,    50,    60,    70,    80,    90,     100,    110,    120,    130,    140,    150;
2.36,  0.30,  0.45,  0.60,  0.72,  0.80,  0.85,  0.87,  0.88,  0.885, 0.89,   0.895,  0.893,  0.885,  0.87,   0.85,   0.82;
2.88,  0.35,  0.50,  0.65,  0.77,  0.84,  0.88,  0.89,  0.895, 0.898, 0.90,   0.906,  0.905,  0.895,  0.88,   0.86,   0.83;
3.56,  0.38,  0.52,  0.67,  0.79,  0.86,  0.90,  0.905, 0.906, 0.905, 0.905,  0.906,  0.905,  0.895,  0.88,   0.86,   0.83;
4.46,  0.36,  0.50,  0.65,  0.76,  0.83,  0.88,  0.895, 0.898, 0.90,  0.902,  0.903,  0.900,  0.89,   0.87,   0.84,   0.81
"""

#==============================================================================*
# Plot

plot(
    charac=convertModelicaTable(tableEtaT),
    dataName="turbine_tableEtaT",
    xlabel="Speed ratio [%]",
    ylabel="Pressure ratio",
    zlabel="Turbine efficiency [%]",
    zscale=100
)
plot(
    charac=convertModelicaTable(tableEtaC),
    dataName="compressor_tableEtaC",
    xlabel="Speed ratio [%]",
    ylabel="Beta",
    zlabel="Compressor efficiency [%]",
    zscale=100
)
plot(
    charac=convertModelicaTable(tablePhicC),
    dataName="compressor_tablePhicC",
    xlabel="Speed ratio [%]",
    ylabel="Beta",
    zlabel="Flow Coefficient"
)
plot(
    charac=convertModelicaTable(tablePR),
    dataName="compressor_tablePR",
    xlabel="Speed ratio [%]",
    ylabel="Beta",
    zlabel="Pressure ratio"
)

#==============================================================================*
