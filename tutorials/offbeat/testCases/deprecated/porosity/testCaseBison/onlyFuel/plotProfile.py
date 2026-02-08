import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import os
import sys

time = 10000
if(len(sys.argv) > 1):
	time = sys.argv[1]

cases = ['UO2Sens', 'MOXLackey', 'vaporPressureMOX']
ls = ['-', '--', '-.']
col = ['b', 'r', 'g']

resultsFile = "sample_porosity_poreVelocityx_T_Q_k_gradTx.xy"


plots = ["p", "v", "T", "Q", "k", "gradT"]
labels = ["Porosity (-)", "Pore v (m/s)", "T (K)", "Q norm (-)", "k (W/m/K)", "gradT"]
cm=1/2.54
fig, axs = plt.subplots(len(plots), 1, sharex=True, figsize = (10, 20))

font = {'family' : 'sans-serif	',
        'size'   : 16}

matplotlib.rc('font', **font)
matplotlib.rcParams['lines.linewidth'] = 3

for i in range(0,len(plots)):
	path = "../resultsBison/"+plots[i]+".csv"
	BIS  = np.genfromtxt(path, delimiter=",")

	if (plots[i] == "Q"):
		BIS[:,1] /= max(BIS[:,1])

	axs[i].plot(BIS[:,0], BIS[:,1], label='Bison', color = 'k')


for j, case in enumerate(cases):
		
	data = np.genfromtxt('postProcessing_{}/sampleRadialProfile/'.format(case)+str(time)+'/'+resultsFile)

	radius = data[:,0]/max(data[:,0])

	# Normalize power
	data[:,4] = data[:,4]/max(data[:,4])

	for i in range(0,len(plots)):

		axs[i].plot(radius, data[:,i+1], label='OFFBEAT - {}'.format(case), ls=ls[j], color=col[j])
		axs[i].set_ylabel(labels[i])

		
		if (i == len(plots)):
			axs[i].set_xlabel("Radius (-)")


	axs[0].set_title('time = {} s'.format(time))
	axs[0].legend(frameon=False)

fig.tight_layout()
plt.savefig("results.pdf")
plt.savefig("results.png")
