import numpy as np

timesteps = []

for i in range(0,10000):
	if(i > 7e3 and i%10 ==0):
		timesteps.append(i)

np.savetxt("timesteps", timesteps, fmt="%d")
