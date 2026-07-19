import numpy as np
import matplotlib.pyplot as plt

B = np.genfromtxt('betaFraction', skip_header=2)
T = np.genfromtxt('T'           , skip_header=2)

Be = np.genfromtxt('betaFractionEq', skip_header=2)
Te = np.genfromtxt('TEq'           , skip_header=2)

Meq = np.genfromtxt('modelEq'   	 , skip_header=0)
Mpm = np.genfromtxt('modelPlusMinus' , skip_header=0)

# Print out value for verification
idx = np.where((Te[:,1] == 1150))[0][0]
relErr = (float(Be[idx,1]) - 0.401681)/0.401681
with open("results", "w") as f:
	f.write("Eq. beta fraction at 1150 K (computed): {:f}".format(float(Be[idx,1])))
	f.write("\n")
	f.write("Eq. beta fraction at 1150 K (expected): 0.401681")
	f.write("\n")
	f.write("Relative error: {:f}\n".format(relErr))

with open("results", "a") as f:
	if ( relErr < 0.005 ):
		f.write("Test passed")
		f.write("\n")
	else:
		f.write("Test failed")
		f.write("\n")

plt.figure()

plt.plot(Meq[:,0], Meq[:,1], label='Massih model', marker='o', ls='', mfc='none', color ='k')
plt.plot(Mpm[:,0], Mpm[:,1], marker='o', ls='', mfc='none', color ='k')

plt.plot(Te[:,1], Be[:,1], label='OFFBEAT, eq', lw = 2, ls='-')
plt.plot(T[:,1], B[:,1], label='OFFBEAT, +-10 K/s', lw = 2, ls='--')

plt.xlim(1000,1300)
plt.ylim(0,1)

plt.grid(ls=':')

plt.xlabel("Temperature (K)")
plt.ylabel("Beta phase fraction (-)")

plt.legend()

plt.title("Zircaloy phase transition model, w=0 ")

plt.savefig("verificationBetaModel.png")
#plt.show()