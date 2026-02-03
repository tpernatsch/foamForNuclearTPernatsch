#!/usr/bin/env python3
#
from matplotlib import pyplot
import numpy as np
import os

# Kinetics parameters
e = 1.6022e-19
NA = 6.0221e23
R = 8.31446
Tssp0 = 30853.0
Tssd0 = 1.02e5
Qp = 0.26*e*NA
Qs = Qp
Qd = 0.37*e*NA
L = 0.0254
Tssp_ = lambda T_: Tssp0 * np.exp(-Qp / (R*T_))
Tssd_ = lambda T_: Tssd0 * np.exp(-Qd / (R*T_))

# Diffusion parameters
D0 = 1.08e-6
Ed = 4.4e4    # Q


def is_float(s):
    try:
        float(s)
    except:
        return False
    
    return True

root = 'postProcessing/graphs'
subdirs = [f.name for f in os.scandir(root) if f.is_dir() and is_float(f.name)]
timeNames = sorted(subdirs, key=lambda x: float(x))

for timeName in timeNames:
    print(f'Creating graph for t={timeName} s')
    x, T, Css, Cpp, Ctot, KG, KN, KD, TSSp, TSSd = np.loadtxt(f'postProcessing/graphs/{timeName}/line_T_Css_Cpp_Ctot_KG_KN_KD_TSSp_TSSd.xy').T
    x = 100*(L-x)
    t = float(timeName)
    Tmin = min(T)
    Tmax = max(T)
    K = (Tmax - Tmin)/L # dT/dx
    
    pyplot.figure(figsize=(5,4))
    pyplot.plot(x, Cpp, label='$C_{pp}$')
    pyplot.plot(x, Css, label='$C_{ss}$')
    pyplot.plot(x, Ctot, label='$C_{tot}$')
    pyplot.plot(x, Tssp_(T), 'k--', label='$T_{ssp}$')
    pyplot.plot(x, Tssd_(T), 'k', label='$T_{ssd}$')
    pyplot.ylabel('C (wt.ppm)')
    pyplot.xlabel('x (cm)')

    # Overlay Sawatsky's analytical solution profiles
    if t >= 3542400:
        # Initial concentration
        NI = -np.sum(0.5*(Ctot[:-1] + Ctot[1:])*np.diff(x)) / 2.54
        
        # Find optimal N0
        N0 = 1.03e5
        def fn(N0):
            N = NI + K**2*D0*N0*(Qd+Qp)/(R*T**4)*((Qd + Ed)/R - 2*T)*t*np.exp(-(Qd+Ed)/(R*T))
            return np.sum(np.abs((Cpp-N)[Cpp>0]))
        
        from scipy.optimize import minimize_scalar
        N0 = minimize_scalar(fn, bounds=[0.1*N0, 10*N0]).x
        
        # 2 phase region
        N = NI + K**2*D0*N0*(Qd+Qp)/(R*T**4)*((Qd + Ed)/R - 2*T)*t*np.exp(-(Qd+Ed)/(R*T))
        N[Css<TSSd] = np.nan
        pyplot.plot(x, N, 'r', label='Analyical')
        
        # Single phase region
        C0 = Ctot[0]/np.exp(Qs/R/T[0])
        N = C0*np.exp(Qs/R/T)
        N[Css>TSSd] = np.nan
        pyplot.plot(x, N, 'r')

    # Finalise plots    
    pyplot.legend(loc=1, framealpha=1)
    pyplot.tight_layout()
    pyplot.savefig(f'C_{t/3600/24:.1f}d.png', dpi=300)
    pyplot.close()

