#!/usr/bin/env python3

import numpy as np
import os

def bending(x, y):
    L = 4
    d = 0.3
    h = 0.6 - d*x/L
    b = 0.3
    I = b*h**3 / 12
    M = 2000*(L-x)
    sigma_max = M*h/2/I
    return sigma_max*(1 - 2*y/h)

def shear(x, y):
    L = 4
    d = 0.3
    h = 0.6 - d*x/L
    b = 0.3
    I = b*h**3 / 12
    V = 2000
    M = V*(L-x)
    dHdx = d/L
    return -M/I/2*(3*y**2/h - 2*y)*dHdx - V/I/2*(-y**2 + h*y)

# Read solution
WM_PROJECT_VERSION = os.environ.get("WM_PROJECT_VERSION", "")
success = True
for x, i in [(0.8, 1), (2.4, 2), (3.8, 3)]:
    if "v" in WM_PROJECT_VERSION:     # equivalent of [[ $WM_PROJECT_VERSION == *v* ]]
        filename = f"postProcessing/xyGraph_ESI/1/line{i}_sigma.xy"
    else:
        filename = f"postProcessing/xyGraph/1/line{i}_sigma.xy"
    data = np.loadtxt(filename)[1:-1].T
    y, sigxx, tauxy, tauxz, sigyy, tauyz, sigzz = data
    sig_vm = np.sqrt(sigxx**2 + sigyy**2 + sigzz**2
                     - sigxx*sigyy - sigxx*sigzz - sigyy*sigzz
                     + 3*(tauxy**2 + tauxz**2 + tauyz**2))
    
    sigxx_analytical = bending(x, y)
    Exx = (sigxx - sigxx_analytical)/5e5
    exx = np.linalg.norm(Exx)/len(Exx)
    
    tauxy_analytical = shear(x, y)
    Exy = (tauxy - tauxy_analytical)/5e4
    exy = np.linalg.norm(Exy)/len(Exy)

    # from matplotlib import pyplot
    # pyplot.plot(y, sigxx/1e3)
    # pyplot.plot(y, sigxx_analytical/1e3)
    # pyplot.plot(y, tauxy/1e3)
    # pyplot.plot(y, tauxy_analytical/1e3)
    
    success &= ((exx < 0.005) & (exy < 0.005))
    print(f"x = {x}, bending stress error norm = {100*exx:.2g}%, shear stress error norm = {100*exy:.2g}%")

if success:
    print("Test passed...")
else:
    print("Test failed...")
    exit(1)
