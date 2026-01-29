import matplotlib.pyplot as plt
import numpy as np
import sys

# === Functions ===

# Temperature profile
def T(r, q3, k, T0):

	Ri=r[0]
	Ro=r[-1]
	
	return T0 + q3/(4*k)*(Ro**2 - r**2) + q3/(2*k)*Ri**2*np.log(r/Ro)

# Gradient T 
def gradT(r, q3, k):

	Ri=r[0]

	return q3/2/k*(Ri**2/r-r)

# Velocity Shape imposed (const in time)
def vel(r, A, l):

	return A*np.exp(-l*r)

# Clement solution for porosity (with given v shape)
def ClementPorosity(r, A, l, t, p0):

	return p0*(1+A*l*t*np.exp(-l*r))**(-1)*(1+1/(l*r)*np.log(1+A*l*t*np.exp(-l*r)))


# === Main ===

porosityTerm=sys.argv[1]

Ri = 0.75e-3
Ro = 3.00e-3
t = Ro - Ri

r = np.linspace(Ri + t/200,Ro - t/200, 200)
r = np.insert(r, 0, Ri)
r = np.append(r, Ro)

# Data for temperature profile
q1 = 35e3
k = 1.5665839180086671
T0 = 1100

# Compute q'''
q3 = q1/((r[-1]**2-r[0]**2)*3.14)

# Compute temperature profile
Temp = T(r, q3, k, T0)

# Parameters for porosity analytical profile
p0 = 0.15
A  = 1e-6
l = 1000

# Plots
times = [ 180, 360]

# Initialize percentage root-mean-square deviation
# PRD is calculated as the RMSE divided by the mean of the observed values and 
# multiplied by 100 to get a percentage value.
PRD = []

fig, axs = plt.subplots(2, sharex=True)
for time in times:
	
	data = np.genfromtxt('postProcessing/sampleRadialProfile/{}/sample_porosity_poreVelocityx.xy'.format(time))
	radius = data[:,0]
	p = data[:,1]
	v = data[:,2]

	# compute PRD
	y_pred = p
	y_true = ClementPorosity(radius, A, l, time, p0)
	prd = 100 * np.sqrt(1/2.0 * np.sum((y_pred - y_true)**2.0) / np.mean(y_true)**2.0)
	PRD.append(prd)

	axs[0].plot(radius, p, ls='', marker='o', mfc='None', markersize='3', label='OFFBEAT - {}'.format(porosityTerm), color = 'r')
	axs[0].plot(r, ClementPorosity(r, A, l, time, p0), label='Clement', color = 'k')
	axs[0].text(0.7*radius[0], max(p), '{:.0f} s'.format(time))
	axs[0].set_ylabel("Porosity (-)")
	axs[0].set_xlim(0.6*radius[0],1.02*radius[-1] )

	if(time==180): 
		axs[0].legend(frameon=False)
		axs[1].plot(radius, v, color ='k', label='OFFBEAT')
		axs[1].plot(r, -vel(r, A, l), color ='r', ls='--', label='Clement')
		axs[1].legend(frameon=False)
		axs[1].set_xlabel('Radius (m)')
		axs[1].set_ylabel('Pore Velocity (m/s)')

axs[0].set_title('Verification Porosity Solver OFFBEAT')
fig.tight_layout()
plt.savefig("VerificationProfiles_{}.pdf".format(porosityTerm))
# plt.show()


# Print PRD results in a file
with open("./results_{}".format(porosityTerm), "w") as f:
	f.write("\tPRD wrt analytical - {} Case (t=180s): {:.4} %\n".format(porosityTerm, PRD[0]))
	f.write("\tPRD wrt analytical - {} Case (t=360s): {:.4} %".format(porosityTerm, PRD[1]))

