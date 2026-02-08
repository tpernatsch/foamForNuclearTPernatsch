#!/bin/env python3
import numpy as np

# Parameters
Pin = 1e6
E1 = 200e9
E2 = 50e9
nu1 = 0.3
nu2 = 0.35
lam1 = E1/2/(1+nu1)
lam2 = E2/2/(1+nu2)

def sigma_xx(x):
    return -1e6*np.ones_like(x)

def sigma_yy(x):
    return (x < 0)*(-4.28571e+05) + (x >= 0)*(-5.38462e+05)

def sigma_xy(x):
    return np.zeros_like(x)


# Read solution
filename = 'postProcessing/xyGraph_ESI/1/radialProfile_sigmaEq_epsilon_sigma_gradD.xy'
data = np.loadtxt(filename)[1:-1].T

# Verify
x = data[0]
sigma_xx_offbeat = data[8]
sigma_yy_offbeat = data[11]
sigma_xx_analytical = sigma_xx(x)
sigma_yy_analytical = sigma_yy(x)

#from matplotlib import pyplot
#pyplot.figure()
#pyplot.plot(x, sigma_xx_analytical)
#pyplot.plot(x, sigma_xx_offbeat, 'x')
#
#pyplot.figure()
#pyplot.plot(x, sigma_yy_analytical)
#pyplot.plot(x, sigma_yy_offbeat, 'x')

e = (sigma_xx_analytical - sigma_xx_offbeat)
maxRelativeErrorXX = np.max(np.abs(e))/Pin
e = (sigma_yy_analytical - sigma_yy_offbeat)
maxRelativeErrorYY = np.max(np.abs(e))/Pin

print('sigma_xx, maximum relative error: %g %%' % (100*maxRelativeErrorXX))
print('sigma_yy, maximum relative error: %g %%' % (100*maxRelativeErrorYY))
if maxRelativeErrorXX < 1e-4 and maxRelativeErrorYY < 1e-4:
    print('Test passed')
else:
    print('Test failed')
    exit(1)
