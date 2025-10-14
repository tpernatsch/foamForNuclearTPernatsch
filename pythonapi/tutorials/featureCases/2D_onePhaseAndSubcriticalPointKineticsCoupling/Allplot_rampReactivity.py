"""
Script to compare the output from GeN-Foam to analytic and finite difference
methods for a ramp increase of reactivity.

Author: Thomas Guilbaud, EPFL/Transmutex SA, 01/06/2023
"""

#=============================================================================*
# Imports
#=============================================================================*

import pandas as pd
import matplotlib.pyplot as plt
import sys
import re
import numpy as np
import numericalSolver

# Test give an output file to analyse
if (len(sys.argv) <= 1):
    print(f'Usage: python {sys.argv[0]} path/to/log.GeN-Foam')
    sys.exit(1)

#=============================================================================*

s = 1
pcm = 1e-5
W = 1
MW = 1e6 * W

MEDIUM_FONT_SIZE = 16

#=============================================================================*

def extractFromLogPK(filename: str, pattern: str, idx: int=2, tmin: float=0, tmax: float=1e9):
    time = 0
    times, values = [], []
    with open(filename, 'r') as file:
        for line in file.readlines():
            if ("Time = " in line and "ExecutionTime" not in line):
                time = float(line.split()[2])
            elif (pattern in line):
                if (tmin <= time and time <= tmax):
                    times.append(time)
                    values.append(float(line.split()[idx]))
    return(times, values)

def errorFormula(target: float, ref: float) -> float:
    return(abs(target-ref)/ref)

def error(target: list, ref: list) -> list:
    """ Relative ratio in % """
    return(list(map(lambda x, y: errorFormula(x, y)*100, target, ref)))

#=============================================================================*
# Reactivity

t0 = 0.01*s
rho_0, r_rho = 0, 1020.408163265307/t0 * pcm/s

def rho(t):
    if (t < t0):
        return(rho_0 + r_rho * t)
    else:
        return(rho_0 + r_rho * t0)
drho_dt = lambda t: r_rho if (t <= t0) else 0

# Prompt Generation Time
Lambda = 0.0015 * s

# One-Precursor Group
betaList   = [0.000075]
lambdaList = [1e-5/s]
betaTot = sum(betaList)  # Delayed neutron fraction
lambda_ = lambdaList[0]  # Decay constant
nPrec = len(lambdaList)

# Subcriticality index
rho_src = lambda ksrc: (ksrc-1)/ksrc
zeta = -rho_src(0.98)
print("zeta = {:.2f} pcm".format(zeta/pcm))

#=============================================================================*
# Initial Neutron Density
N0 = 10 * MW
# --- Precursor equilibrium
initialPopulation = [N0]+[betaList[i]*N0/(lambdaList[i]*Lambda) for i in range(len(lambdaList))]

# --- Population Matrix
def populationMatrix(t):
    matrix = np.zeros((nPrec+1, nPrec+1))
    matrix[0][0] = (rho(t)-zeta-betaTot)/Lambda
    for i, lambdaTemp in enumerate(lambdaList):
        matrix[0][i+1] = lambdaTemp
        matrix[i+1][0] = betaList[i]/Lambda
        matrix[i+1][i+1] = -lambdaTemp
    return(matrix)

# --- External Source Vector
q_0 = N0/Lambda * (abs(rho_0) + zeta) # neutrons/s
r_q, r_sq = 0, 0
print("Source definition: S_0 = {:.2f}; r_S = {:.2f}".format(q_0*Lambda, r_q*Lambda))
print("    Modulation: m_S(t=t0) = {:.6f}".format(t0*r_q/q_0))
omega = 1
def q(t):
    if (t <= t0):
        return(q_0 + r_q * t + r_sq*np.sin(omega*t))
    else:
        return(q_0 + r_q * t0 + r_sq*np.sin(omega*t0))

def dq_dt(t):
    if (t <= t0):
        return(r_q + r_sq*omega*np.cos(omega*t))
    else:
        return(0)

def extSrcVec(t):
    return(np.array([q(t)] + [0]*nPrec))


#=============================================================================*
# Substitute Point-Kinetics equation

# From [1], works only for 1 precursor group
def dN_dt_substitute(t: float, N: float) -> float:
    A = lambda_ * (rho(t) - zeta) + drho_dt(t)
    B = -rho(t) + betaTot + lambda_*Lambda + zeta
    C = lambda_*Lambda*q(t) + Lambda*dq_dt(t)
    return(A/B * N + C/B)

def solverNsubstitute(t: float, dt: float, N: float):
    return(N + dt * dN_dt_substitute(t, N))


#=============================================================================*
# Analytical solution
# From [3]
def NconstSrcAnalytic(t):
    f = lambda tp: (betaTot*N0 + Lambda*q_0)/(betaTot - r_rho*tp + abs(rho_0-zeta))
    if (t >= t0):
        a = Lambda*q_0/(abs(rho_0-zeta) - r_rho*t0)
        k1 = (abs(rho_0-zeta)*lambda_ - r_rho*t0*lambda_)/(betaTot+abs(rho_0-zeta)-r_rho*t0)
        return(a*(1 - np.exp(-k1 * (t-t0))) + f(t0)*np.exp(-k1 * (t-t0)))
    else:
        return(f(t))

#=============================================================================*

def plot(ax, axZoom, axError, tRange: list, yRange: list, label: str,
    marker: str='',
    markevery: int=None
):
    ax.plot(tRange, yRange, label=label, marker=marker, markevery=markevery)

    yRef = [NconstSrcAnalytic(t) for t in tRange]
    errorRange = error(yRange, yRef)
    axError.plot(tRange, errorRange, label=label, marker=marker, markevery=markevery)

    axZoom.plot(
        tRange, yRange, label=label, marker=marker,
        markevery=(markevery[0], int(markevery[1]/10)) if markevery != None else None
    )


#=============================================================================*

logFileName = sys.argv[1]

# Extract from log
timeLog, powerPK = extractFromLogPK(logFileName, "totalPower", 2)

# Offset the time to start at 0
timeLog = [t-timeLog[0] for t in timeLog]

fig, (ax, axError) = plt.subplots(2, sharex=True, figsize=(14, 9))
axZoom = fig.add_axes([0.4, 0.6, 0.2, 0.15]) # [x0, y0, width, height]

plot(ax, axZoom, axError, timeLog, powerPK, label="GeN-Foam Subcritical Point-Kinetics")

#=============================================================================*

# Set time range for the numerical solvers
dt = 1e-2 * s
tmin, tmax = 0, timeLog[-1]
tRange = np.arange(tmin, tmax+dt, dt)

# --- Euler
# tRangeDirect, NdirectRecord = numericalSolver.solverNeutronPopulation(
#     tmin, tmax, dt,
#     lambda t, t0, x: numericalSolver.solverEuler(t, t0, x, populationMatrix(t), extSrcVec(t)),
#     initialPopulation=initialPopulation
# )
# plt.plot(tRangeDirect, [e[0] for e in NdirectRecord], label="Finite difference Euler")

# --- Runge-Kutta 4
tRangeRK4, NRK4Record = numericalSolver.solverNeutronPopulation(
    tmin, tmax, dt,
    lambda t, t0, x: numericalSolver.solverRungeKutta4(t, t0, x, populationMatrix(t), extSrcVec(t)),
    initialPopulation=initialPopulation
)
NRK4Record = [e[0] for e in NRK4Record]
plot(
    ax, axZoom, axError, tRangeRK4, NRK4Record,
    label="Finite difference Runge-Kutta 4", marker="o", markevery=(250, 500)
)

# --- Finite Difference Substitution
tRangeNsubstitute, NsubstituteRecord = numericalSolver.solverNeutronPopulation(
    tmin, tmax, dt, solverNsubstitute, initialPopulation=N0
)
plot(ax, axZoom, axError, tRangeNsubstitute, NsubstituteRecord, label="Finite difference substitution (E. Henrice et al.)", marker="x", markevery=(0, 500))

# --- Matrix Exponantial
tRangeNconstSrcExpMat, NconstSrcExpMatRecord = numericalSolver.solverNeutronPopulation(
    tmin, tmax, dt,
    lambda t, t0, x: numericalSolver.NconstSrcExpMat(t, t0, x, populationMatrix(t), extSrcVec(t)),
    initialPopulation=initialPopulation
)
NconstSrcExpMatRecord = [e[0] for e in NconstSrcExpMatRecord]
plot(ax, axZoom, axError, tRangeNconstSrcExpMat, NconstSrcExpMatRecord, label="Finite difference exponential (M. Kinard)", marker="*", markevery=(25, 500))


# --- Analytic (reference)
powerRef = [NconstSrcAnalytic(t) for t in tRange]
ax.plot(tRange, powerRef, '--', label="Analytic (E. Henrice et al.)")
axZoom.plot(tRange, powerRef, '--', label="Analytic (E. Henrice et al.)")



#=============================================================================*
# --- Final step comparison to analytic formula

NlastTh = NconstSrcAnalytic(tmax)
print("| Method              | Value [W]    | Error [%]    |")
print("|:--------------------|:------------:|:------------:|")
print("| Infinite analytic   | {:.6e} | -            |".format(NconstSrcAnalytic(1e16*s)))
print("| Last value analytic | {:.6e} | -            |".format(NlastTh))
print("|            GeN-Foam | {:.6e} | {:.6e} |".format(powerPK[-1], errorFormula(powerPK[-1], NlastTh)*100))
print("|       Runge-Kutta 4 | {:.6e} | {:.6e} |".format(NRK4Record[-1], errorFormula(NRK4Record[-1], NlastTh)*100))
print("|  Euler Substitution | {:.6e} | {:.6e} |".format(NsubstituteRecord[-1], errorFormula(NsubstituteRecord[-1], NlastTh)*100))
print("|  Matrix Exponential | {:.6e} | {:.6e} |".format(NconstSrcExpMatRecord[-1], errorFormula(NconstSrcExpMatRecord[-1], NlastTh)*100))

#=============================================================================*

ax.set_xlabel("Time [s]", fontsize=MEDIUM_FONT_SIZE)
ax.set_ylabel("Fission power [W]", fontsize=MEDIUM_FONT_SIZE)
ax.tick_params(axis='x', labelsize=MEDIUM_FONT_SIZE)
ax.tick_params(axis='y', labelsize=MEDIUM_FONT_SIZE)
ax.grid(True)
ax.legend()

axZoom.set_xlim(-0.1*s, 1.7*s)
axZoom.set_ylim(0.994*NconstSrcAnalytic(t0), 1.003*NconstSrcAnalytic(t0))

axError.set_xlabel("Time [s]", fontsize=MEDIUM_FONT_SIZE)
axError.set_ylabel("Error [%]", fontsize=MEDIUM_FONT_SIZE)
axError.set_yscale("log")
axError.tick_params(axis='x', labelsize=MEDIUM_FONT_SIZE)
axError.tick_params(axis='y', labelsize=MEDIUM_FONT_SIZE)
axError.grid(True)
axError.legend()

plt.rcParams['font.size'] = 12
fig.savefig("fig_results_rampReactivity.png")
plt.show()

#=============================================================================*
