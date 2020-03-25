from matplotlib import pyplot as plt
import os

startTime = 200
fig, (ax1, ax2) = plt.subplots(2,1)

with open("GeN-Foam.dat") as f:
	time = []
	power = []
	avFuelTemp = []
	lines = f.readlines()
	for line in lines :
		splitLine = line.split(";")
		try :
			time.append(float(splitLine[0])-startTime)
			power.append(float(splitLine[2]))
			avFuelTemp.append(float(splitLine[6]))
		except :
			pass

	ax1.semilogx(time, power)
	ax1.set_xlabel("Time (s)")
	ax1.set_ylabel("Power (W)")
	ax1.grid(True)
	ax2.semilogx(time, avFuelTemp)
	ax2.set_xlabel("Time (s)")
	ax2.set_ylabel("Average Fuel Temperature (K)")
	ax2.grid(True)
	plt.subplots_adjust(hspace=0.5)
	plt.show()
