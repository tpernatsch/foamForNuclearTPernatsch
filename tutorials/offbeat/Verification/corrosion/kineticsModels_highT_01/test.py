# Verification for high temperature oxidation models
import numpy as np
from matplotlib import pyplot

t = 1000
R = 1.987
rho = 0.72

#Baker/Just
T = np.linspace(700, 2500, 100)
dOxide = 1e6*np.sqrt(1.883e-4*t*np.exp(-45500/R/T))
pyplot.plot(1e4/T, dOxide, 'b', label='Baker-Just')

T1, dOxide1 = np.loadtxt('Baker-Just.out').T
dOxide1 *= 1e6
pyplot.plot(1e4/T1, dOxide1, 'bx', mfc="none", label='OFFBEAT', markevery=3)

# Cathcart Pawel
T = np.linspace(700, 1900, 10)
dOxide = 1e6*np.sqrt(2*1.12569e-6*np.exp(-18062/T)*t)
pyplot.plot(1e4/T, dOxide, 'r', label='Cathcart-Pawel')

# Prater-Courtright
T = np.linspace(1800, 2500, 10)
dOxide = 1e6*np.sqrt(2.98116e-3*t*np.exp(-28420/T))
pyplot.plot(1e4/T, dOxide, 'r--', label='Prater-Courtright')

T2, dOxide2 = np.loadtxt('Cathcart-Pawel-Prater-Courtright.out').T
dOxide2 *= 1e6
pyplot.plot(1e4/T2, dOxide2, 'ro', mfc="none", label='OFFBEAT', markevery=(1,3))

# Leistikow
T = np.linspace(700, 1900, 10)
dOxide = 1e6*np.sqrt(7.82e-6*t*np.exp(-20214/T))
pyplot.plot(1e4/T, dOxide, 'g--', label='Leistikow')

T3, dOxide3 = np.loadtxt('Leistikov-Prater-Courtright.out').T
dOxide3 *= 1e6
pyplot.plot(1e4/T3, dOxide3, 'gd', mfc="none", label='OFFBEAT', markevery=(2,3))

pyplot.semilogy()
pyplot.xlabel(r"10'000/T $\left(\frac{10^4}{K}\right)$")
pyplot.ylabel('Oxide Thickness (um)')
pyplot.legend()
pyplot.savefig('dOxide.png', dpi=300)

e1 = 1 - dOxide1/(1e6*np.sqrt(1.883e-4*t*np.exp(-45500/R/T1)))
assert np.max(np.abs(e1)) < 1e-3

e2a = 1 - dOxide2/(1e6*np.sqrt(2*1.12569e-6*np.exp(-18062/T2)*t))
assert np.max(np.abs(e2a[:35])) < 1e-3

e2b = 1 - dOxide2/(1e6*np.sqrt(2.98116e-3*t*np.exp(-28420/T2)))
assert np.max(np.abs(e2b[38:])) < 1e-3

e3 = 1 - dOxide3/(1e6*np.sqrt(7.82e-6*t*np.exp(-20214/T3)))
assert np.max(np.abs(e3[:35])) < 1e-3

print("Success")
