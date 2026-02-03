import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import os
import sys

time = 25200
if(len(sys.argv) > 1):
	time = sys.argv[1]

resultsFile = "sample_porosity_poreVelocityx_T_Q_k.xy"

plots = ["p", "v", "T", "Q", "k"]
labels = ["Porosity (-)", "Pore v (m/s)", "T (K)", "Q norm (-)", "k (W/m/K)"]
cm=1/2.54
fig, axs = plt.subplots(len(plots), 1, sharex=True, figsize = (10, 20))

font = {'family' : 'sans-serif	',
        'size'   : 16}

matplotlib.rc('font', **font)
matplotlib.rcParams['lines.linewidth'] = 3


data = np.genfromtxt(f'postProcessing/sampleRadialProfile/{str(time)}/{resultsFile}')

radius = data[:,0]/max(data[:,0])

# Normalize power
data[:,4] = data[:,4]/max(data[:,4])

for i in range(0,len(plots)):

	axs[i].plot(radius, data[:,i+1])
	axs[i].set_ylabel(labels[i])
	

axs[-1].set_xlabel("Radius (-)")
axs[0].set_title('time = {} s'.format(time))

fig.tight_layout()
plt.savefig("results.pdf")
