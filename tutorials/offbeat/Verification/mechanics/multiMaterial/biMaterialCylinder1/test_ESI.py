#!/bin/env python3
import numpy as np

# Parameters
Pin = 1e7
E1 = 200e9
E2 = 50e9
nu1 = 0.3
nu2 = 0.35

r1 = 0.01
r2 = 0.015
r3 = 0.02

# Analytical solution
P12 = (2*r1**2*Pin/(E1*(r2**2-r1**2)) / (
        1/E2*((r3**2+r2**2)/(r3**2-r2**2)+nu2)
      + 1/E1*((r2**2+r1**2)/(r2**2-r1**2)-nu1)
))

def sigma_r1(r):
    return ((r1**2*Pin - r2**2*P12 + (P12-Pin)*(r1*r2/r)**2)\
            / (r2**2 - r1**2))

def sigma_r2(r):
    return ((r2**2*P12 - P12*(r2*r3/r)**2) / (r3**2 - r2**2))

def sigma_r(r):
    return (r < r2)*sigma_r1(r) + (r >= r2)*sigma_r2(r)

def sigma_h1(r):
    return ((r1**2*Pin - r2**2*P12 - (P12-Pin)*(r1*r2/r)**2)\
            / (r2**2 - r1**2))

def sigma_h2(r):
    return ((r2**2*P12 + P12*(r2*r3/r)**2) / (r3**2 - r2**2))

def sigma_h(r):
    return (r < r2)*sigma_h1(r) + (r >= r2)*sigma_h2(r)


# Read solution
filename = 'postProcessing/xyGraph_ESI/1/radialProfile_sigmaEq_D_epsilon_sigmaCyl_gradD.xy'
data = np.loadtxt(filename)[1:-1].T

# Verify
r = data[0] + r1
sigma_r_offbeat = data[11]
sigma_h_offbeat = data[14]
sigma_r_analytical = sigma_r(r)
sigma_h_analytical = sigma_h(r)

#from matplotlib import pyplot
#pyplot.figure()
#pyplot.plot(r, sigma_r_analytical)
#pyplot.plot(r, sigma_r_offbeat, 'x')
#
#pyplot.figure()
#pyplot.plot(r, sigma_h_analytical)
#pyplot.plot(r, sigma_h_offbeat, 'x')

e = (sigma_r_analytical - sigma_r_offbeat)
maxRelativeErrorR = np.max(np.abs(e))/Pin
e = (sigma_h_analytical - sigma_h_offbeat)
maxRelativeErrorH = np.max(np.abs(e))/Pin

print('Radial stress, maximum relative error: %g %%' % (100*maxRelativeErrorR))
print('Hoop stress, maximum relative error: %g %%' % (100*maxRelativeErrorH))
if maxRelativeErrorR < 0.01 and maxRelativeErrorH < 0.01:
    print('Test passed')
else:
    print('Test failed')
    exit(1)
