# - Internal Gas Pressure (MPa) :
import subprocess as sp
import numpy as np
import os

rip_string  = sp.check_output("cat log.offbeat | grep -E 'Gas pressure: [+-]?[0-9]+([.][0-9]+)?'| cut -d' ' -f3 > rip.tmp", shell=True).decode("utf-8")
iter_string = sp.check_output("cat log.offbeat | grep -E 'OuterIteration n.'| cut -d' ' -f3 > iter.tmp", shell=True).decode("utf-8")

ripData  = np.genfromtxt("rip.tmp", skip_header = 1)
iterData = np.genfromtxt("iter.tmp")
os.remove('rip.tmp')
os.remove('iter.tmp')

pGap = []
for i in range(len(iterData)):
    if( i == len(iterData)-1):
        pGap.append(ripData[-1])
        break
    if( iterData[i] > iterData[i+1]):
        pGap.append(ripData[i])

with open("gasPressures", "a") as f:
    f.write("Internal gas pressure (MPa) \n")
    f.write("time(days)/internal gas pressure \n")
    for i in range(0,len(pGap)):
        f.write(str(pGap[i]) + " \n")
    f.write("\n")
