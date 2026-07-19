import matplotlib.pyplot as plt
import numpy as np

# === Functions ===

# Temperature profile
def T(r, q3, k, T0):

	Ri=r[0]
	Ro=r[-1]
	
	return T0 + q3/(4*k)*(Ro**2 - r**2) + q3/(2*k)*Ri**2*np.log(r/Ro)

# Gradient 
def gradT(q3, k, r):

	Ri=r[0]

	return q3/2/k*(Ri**2/r-r)

# Laplacian
def laplT(q3, k):

	return -q3/k

# Clement's analytical solution
def ClementSolution(T, c0, q3, k, r, t, D):

	U = -17623/T*((55891/T-2)*(gradT(q3, k,r)/T)**2+laplT(q3, k)/T)

	return (c0 + c0*(1-c0)*D*U*t)


# === Main ===

Ri = 0.75e-3
Ro = 3.00e-3
t  = Ro-Ri

r = np.linspace(Ri + t/200,Ro - t/200, 100)
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

# Print out file with temperature for OFFBEAT
# np.savetxt('radialTempProfileInput.csv', Temp, fmt='%.2f')

#plt.step(r, Temp, where='mid')
#plt.show()
#quit()

# Initial concentrations
c0Pu = 0.2
c0Am = 0.02

# Diffusion coeffs
DPu = 3.4e-5*np.exp(-55891/Temp)
DAm = 2.0e-6*np.exp(-46324/Temp)

# Plots
# 			1h       5h
times = [ 360000, 1800000]

fig, axs = plt.subplots(3, sharex=True)

for time in times:

	data = np.genfromtxt('postProcessing/sampleRadialProfile/{}/sample_AmFormFactor_PuFormFactor_T.xy'.format(time))

	radius = data[:,0]
	Pu = data[:,2]
	Am = data[:,1]
	T = data[:,3]

	axs[0].plot(radius, Pu, ls='', marker='o', mfc='None', markersize='3', label='OFFBEAT', color = 'r')
	axs[0].plot(r, ClementSolution(Temp, c0Pu, q3, k, r, time, DPu), label='Analitycal', color = 'k')
	axs[0].text(0.7*radius[0], 1.002*max(Pu), '{:.0f} h'.format(time/3600))
	axs[0].set_ylabel("Pu (-)")
	axs[0].set_xlim(0.6*radius[0],1.02*radius[-1] )

	axs[1].plot(radius, Am, ls='', marker='o', mfc='None', markersize='3', label='OFFBEAT', color = 'r')
	axs[1].plot(r, ClementSolution(Temp, c0Am, q3, k, r, time, DAm), label='Analitycal', color = 'k')
	axs[1].text(0.7*radius[0], 1.002*max(Am), '{:.0f} h'.format(time/3600))
	axs[1].set_ylabel("Am(-)")
	axs[1].set_xlim(0.6*radius[0],1.02*radius[-1] )

	if(time==360000): 
		axs[0].legend(frameon=False)

		# Interpolate analytical solution to match OFFBEAT radii
		Pu_interpolated = np.interp(radius, r, ClementSolution(Temp, c0Pu, q3, k, r, time, DPu))
		PuRMSE = np.sqrt(np.mean((Pu - Pu_interpolated) ** 2))

		# Interpolate analytical solution to match OFFBEAT radii
		Am_interpolated = np.interp(radius, r, ClementSolution(Temp, c0Am, q3, k, r, time, DAm))
		AmRMSE = np.sqrt(np.mean((Am - Am_interpolated) ** 2))


axs[2].plot(radius, T,  ls='', marker='o', mfc='None', markersize='3', label='OFFBEAT', color = 'r')
axs[2].plot(r, Temp, color ='k')
axs[2].set_xlabel('Radius (m)')
axs[2].set_ylabel('Temperature (K)')

axs[0].set_title('Verification Redistribution Model OFFBEAT')
fig.tight_layout()
plt.savefig("VerificationProfiles.png")

with open("results","w") as f:
	f.write(f"RMSE Pu Distribution @ 100 h : {PuRMSE*100:.2e} %")
	f.write("\n")
	f.write(f"RMSE Am Distribution @ 100 h : {AmRMSE*100:.2e} %")
	f.write("\n")


