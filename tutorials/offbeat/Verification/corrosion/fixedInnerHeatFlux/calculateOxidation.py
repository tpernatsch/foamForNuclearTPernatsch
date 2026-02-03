# Importing libraries:
import numpy as np
import subprocess
import math
import matplotlib.pyplot as plt
# plt.style.use(['science','ieee'])

# --------------------------------------------------------------------------- #
#                                MAIN                                         #
# --------------------------------------------------------------------------- #

t = 0
deltaT = 1

S = 0
Sold = 0
Stransition = 2e-6

Q1R = 16266
Q2R = 13775

C1 = 6.30E-09
C2 = 215.758168983928

heatflux = 298900
To = 600
phi = 1e13
kox = 0
Ti = 550

while t < 1e3:

	Sold = S

	err = 1
	counter = 0

	while err > 1e-16:

		Sprev = S

		kox = 0.835 + 1.81e-4 * Ti; 

		Ti = To + heatflux*S/kox

		if S <= Stransition:
			
			S = pow(C1*math.exp(-Q1R/Ti)*deltaT + pow(Sold, 3.0), 1.0/3.0)

		else:
			
			S = pow(C2*math.exp(-Q2R/Ti)*deltaT + pow(Sold, 1.0), 1.0/1.0)

		err = (S-Sprev)/S
		counter += 1

	t += deltaT

print 'Final oxide thickness: ' + str(S) + ' m'
print 'Final interface temperature: ' + str(Ti) + ' K'