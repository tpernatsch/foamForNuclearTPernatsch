# Importing libraries:
import numpy as np
import matplotlib.pyplot as plt
import sys

lsDataTime = []
lsDataStress = []

plasticity = str(sys.argv[1])

file = "postProcessing/innerStress/0/sigmaCyl"
with open(file) as f:
    lines = f.readlines()

    rangeLines = len(lines)
    for line in lines[4:rangeLines]:
        # if(float(line.split()[1]) > 1000):
        lsDataTime.append(float(line.split()[0])+10)
        lsDataStress.append(-1*float((line.split()[1]).lstrip("(")))

# Computation analytical profile
sigmaY = 0.5e6
a0 = 10e-3
b0 = 20e-3
r0 = 10e-3
a  = np.linspace(a0,a0+75e-3,100);


def analyticalProfile(x):
    return -sigmaY/np.sqrt(3)*np.log(( (r0/a0)**2 + (x/a0)**2 - 1) /( (b0/a0)**2 + (x/a0)**2 - 1) )

# Evaluate error between analytical and numerical solution
relErr = 0
for i in range(0,len(lsDataStress)):
    relErr = max(
                    abs(lsDataStress[i]-analyticalProfile(lsDataTime[i]*1e-3))
                    /analyticalProfile(lsDataTime[i]*1e-3)
                    ,relErr
                )

f=open('relativeErrors','a')
f.write('max relErr '+plasticity+' : {:.2f} %\n'.format(relErr*100))
f.close()

plt.figure()
# plt.plot(ssData[:,0]*1e3,ssData[:,1],label='smallStrain',color='b')
plt.plot(lsDataTime,lsDataStress,label=plasticity,color='r', linewidth= 10 )
plt.plot(a*1000,analyticalProfile(a),label='analytical',color='k')
plt.xlabel('radial coordinate (mm)')
plt.ylabel('Stress (Pa)')
plt.legend()
plt.savefig('Stresses_'+plasticity+'.png')
