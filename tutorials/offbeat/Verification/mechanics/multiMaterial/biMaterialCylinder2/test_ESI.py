#!/bin/env python3
import numpy as np

# Parameters
tau_o = 1.59174e+06
E1 = 200e9
E2 = 50e9
nu1 = 0.3
nu2 = 0.35
K1 = E1/3/(1-2*nu1)
K2 = E2/3/(1-2*nu2)
lamda1 = E1*nu1/(1+nu1)/(1-2*nu1)
lamda2 = E2*nu2/(1+nu2)/(1-2*nu2)
mu1 = 1.5*(K1-lamda1)
mu2 = 1.5*(K2-lamda2)

r1 = 0.01
r2 = 0.015
r3 = 0.02

# Analytical solution for bimaterial pipe under torque
def tau(r):
    return tau_o*r3**2 / r**2

def dy1(r):
    return r3**2*tau_o / (2*mu1) * (r/r1**2 - 1/r)

def dy2(r):
    return r3**2*tau_o/2 * ((1/mu1/r1**2 - 1/mu1/r2**2 + 1/mu2/r2**2)*r - 1/mu2/r)

def dy(r):
    return (r < r2)*dy1(r) + (r >= r2)*dy2(r)

# Shear stress
filename = 'postProcessing/xyGraph_ESI/1/radialProfile_sigmaEq_D_epsilon_sigmaCyl_gradD.xy'
data = np.loadtxt(filename)[1:-1].T

# Verify
r = data[0] + r1
tau_offbeat = data[12]
tau_analytical = tau(r)

#from matplotlib import pyplot
#pyplot.figure()
#pyplot.plot(r, tau_analytical)
#pyplot.plot(r, tau_offbeat, 'x')

e = (tau_analytical - tau_offbeat)
maxRelativeErrorTau = np.max(np.abs(e))/tau_o

# Displacement
filename = 'postProcessing/xyGraph_ESI/1/radialProfile_sigmaEq_D_epsilon_sigmaCyl_gradD.xy'
data = np.loadtxt(filename)[1:-1].T

r = data[0] + r1
D_offbeat = data[3]
D_analytical = dy(r)

#pyplot.figure()
#pyplot.plot(r, D_analytical)
#pyplot.plot(r, D_offbeat, 'x')

e = (D_analytical - D_offbeat)
maxRelativeErrorD = np.max(np.abs(e))/np.average(D_analytical)

print('Shear stress, maximum relative error: %g %%' % (100*maxRelativeErrorTau))
print('Displacement, maximum relative error: %g %%' % (100*maxRelativeErrorD))
if maxRelativeErrorTau < 0.01 and maxRelativeErrorD < 0.001:
    print('Test passed')
else:
    print('Test failed')
    exit(1)
