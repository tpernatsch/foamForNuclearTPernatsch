"""
Analytical solution for the subcritical point-kinetics solver development in
GeN-Foam
Author: Thomas Guilbaud, EPFL/Transmutex SA
Last Update: 2022/04/05

Reference:
[2] Matthew Kinard, Efficient numerical solution of the point kinetics equation,
Thesis, December 2003
"""

# Imports ----------------------------------------------------------------------*
import numpy as np
import time
from scipy.linalg import expm

# Classical expression ---------------------------------------------------------*

# Euler
def solverEuler(t: float, dt: float, x: list, mat: list, src: list) -> list:
    dx_dt = np.dot(mat, x) + src
    return(x + dt * dx_dt)

# Runge-Kutta 4
def solverRungeKutta4(t: float, dt: float, x: list, mat: list, src: list) -> list:
    dx_dt = lambda t_, x_: np.dot(mat, x_) + src
    k1 = dx_dt(t, x)
    k2 = dx_dt(t+dt/2, x+k1*dt/2)
    k3 = dx_dt(t+dt/2, x+k2*dt/2)
    k4 = dx_dt(t+dt, x+k3*dt)
    return(x + dt/6 * (k1 + 2*k2 + 2*k3 + k4))

#-------------------------------------------------------------------------------*
# Analytical solution
# From [2]
def NconstSrcExpMat(t: float, dt: float, x, mat: list, src: list) -> list:
    popMat = mat
    expMat = expm(popMat*dt)
    return(
        np.dot(expMat, x)
        + np.dot(
            expMat - np.eye(2),
            np.dot(np.linalg.inv(popMat), src)
        )
    )


# General Solver ---------------------------------------------------------------*

def solverNeutronPopulation(tmin, tmax, dt, solver, initialPopulation):
    t1 = time.time()
    tRange = np.arange(tmin, tmax+dt, dt)
    popRecord = [initialPopulation]
    for t in tRange[1:]:
        popRecord.append(solver(t, dt, popRecord[-1]))
    print("time = {:.6f} s  ".format(time.time()-t1))
    return(tRange, popRecord)
