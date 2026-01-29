#!/bin/env python3
import numpy as np

# Parameters
sigma0 = 1e6
R = 0.5
H = 4

# Analytical stress concentration factor
Kc = 3 - 3.14*R/H + 3.667*(R/H)**2 - 1.527*(R/H)**3
sigma_avg = H/(H-R)*sigma0
sigma_max = Kc*sigma_avg

# Analytical solution for infinite flat plate
def sigma_xx(y):
    return sigma0*(1 + R**2/2/y**2 + 3*R**4/2/y**4)

# Read solution
filename = 'postProcessing/xyGraph_ESI/1/symmetry_sigmaEq_epsilon_sigma_gradD.xy'
data = np.loadtxt(filename)[1:-1].T

# Verify
y = data[0]
sigma_xx_offbeat = data[8]
sigma_xx_analytical = sigma_xx(y)

#from matplotlib import pyplot
#pyplot.figure()
#pyplot.plot(y, sigma_xx_analytical)
#pyplot.plot(y, sigma_xx_offbeat, 'x')

Kc_offbeat = np.max(sigma_xx_offbeat)/sigma_avg
maxRelativeError = abs(Kc - Kc_offbeat)/Kc

print('Analytical Kc, %.3f' % Kc)
print('Kc predicted by OFFBEAT, %.3f' % Kc_offbeat)
print('Maximum relative error: %g %%' % (100*maxRelativeError))
if maxRelativeError < 0.01:
    print('Test passed')
else:
    print('Test failed')
    exit(1)
