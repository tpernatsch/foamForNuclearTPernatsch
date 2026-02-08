import numpy as np
import matplotlib.pyplot as plt

Vss = "./Vss"
Vls = "./Vls"

# Import volumes and normalize them
ssData = np.genfromtxt("Vss", delimiter='\n')
lsData = np.genfromtxt("Vls", delimiter='\n')
ssData /= ssData[0]
lsData /= lsData[0]

# Import displacements @ inner radius
Dss = []
file = "./Dss"
with open(file) as f:
    lines = f.readlines()
    rangeLines = len(lines)
    for line in lines[2:rangeLines]:
        Dss.append(float((line.split()[1]).lstrip("(")))

Dls = []
file = "./Dss"
with open(file) as f:
    lines = f.readlines()
    rangeLines = len(lines)
    for line in lines[2:rangeLines]:
        Dls.append(float((line.split()[1]).lstrip("(")))

# Initial radius
r0 = 0.01

# Deformed radius / initial radius
rss = 1 + np.divide(Dss,r0)
rls = 1 + np.divide(Dls,r0)

plt.figure()
plt.hlines(1,0.8*rss[0],rss[-1]*1.2, label='initial volume', ls=':', color='k')
plt.plot(rss,ssData, label="smallStrain")
plt.plot(rls,lsData, label="largeStrainUpd")
plt.xlim(0.95*rss[0],1.05*rss[-1])
plt.ylim(0.7,1.3)
plt.xlabel(r'$R_{in}^{def}/R_{in}^{0}$ (-)')
plt.ylabel('Total Mesh Volume Normalized (-)')
plt.title('Expansion of a cylinder - volume conservation')
plt.legend()
plt.savefig("volumes.png")
#plt.show()

