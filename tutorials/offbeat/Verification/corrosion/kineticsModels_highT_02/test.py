# Verification for high temperature oxidation models
import numpy as np
from matplotlib import pyplot

R = 1.987
t = np.linspace(0,10000,500)
T = 1073.15

#Baker/Just
dOxide = 1e6*np.sqrt(1.883e-4*t*np.exp(-45500/R/T))
pyplot.plot(t, dOxide, 'b', label='Baker-Just')

t1, dOxide1 = np.loadtxt('Baker-Just.out').T
dOxide1 *= 1e6
pyplot.plot(t1, dOxide1, 'bx', mfc="none", label='OFFBEAT', markevery=(1,3))

# Cathcart Pawel
dOxide = 1e6*np.sqrt(2*1.12569e-6*np.exp(-18062/T)*t)
pyplot.plot(t, dOxide, 'r', label='Cathcart-Pawel')

t2, dOxide2 = np.loadtxt('Cathcart-Pawel-Prater-Courtright.out').T
dOxide2 *= 1e6
pyplot.plot(t2, dOxide2, 'ro', mfc="none", label='OFFBEAT', markevery=(2,3))

# Leistikow
dOxide = 1e6*np.sqrt(7.82e-6*t*np.exp(-20214/T))
pyplot.plot(t, dOxide, 'g--', label='Leistikow')

t3, dOxide3 = np.loadtxt('Leistikov-Prater-Courtright.out').T
dOxide3 *= 1e6
pyplot.plot(t3, dOxide3, 'gd', mfc="none", label='OFFBEAT', markevery=3)

pyplot.gca().set_xscale('function', functions=(lambda x: np.power(x,0.5), lambda x: np.power(x, 2)))
pyplot.xlabel("Time (s)")
pyplot.ylabel('Oxide Thickness (um)')
pyplot.legend()
pyplot.xlim(xmin=0)
pyplot.ylim(ymin=0)
pyplot.savefig('dOxide.png', dpi=300)

dOxide = (1e6*np.sqrt(1.883e-4*t1*np.exp(-45500/R/T)))
e1 = (dOxide1-dOxide)/np.maximum(dOxide, 1e-12)
assert np.max(np.abs(e1)) < 1e-3

dOxide = (1e6*np.sqrt(2*1.12569e-6*np.exp(-18062/T)*t2))
e2 = (dOxide2-dOxide)/np.maximum(dOxide, 1e-12)
assert np.max(np.abs(e2[:35])) < 1e-3

dOxide = (1e6*np.sqrt(7.82e-6*t3*np.exp(-20214/T)))
e3 = (dOxide3-dOxide)/np.maximum(dOxide, 1e-12)
assert np.max(np.abs(e3[:35])) < 1e-3

print("Success")
