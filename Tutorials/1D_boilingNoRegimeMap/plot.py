### IMPORTS

import sys
import matplotlib.pyplot as plt

### FUNCTIONS

def limit(x, y) :
    
    tmpX = []
    i = 0
    while len(tmpX) < len(y) :
        tmpX.append(x[i])
        i += 1
    return tmpX

def filter(x, s) :

    tmpX = []
    n = len(x)
    i = 0
    while i < n :
        j = i
        avg = 0
        ss = 0
        while j < i + s :
            try :
                avg += x[j]
                ss += 1
                j += 1
            except :
                break
        avg /= ss
        j = i
        while j < i + ss :
            tmpX.append(avg)
            j += 1
        i = j
    return tmpX

### MAIN

fig, ax1 = plt.subplots(1)

ax1.set_ylabel("initialResidual (-)")
ax1.set_xlabel("time (s)")

for q in range(len(sys.argv)-1) :

    filename = sys.argv[q+1]
    file = open(filename, "r")
    lines = file.readlines()
    n = 0;
    for line in lines :
        if "Outer iteration " in line :
            nn = int(line.split()[2])
            if nn > n :
                n = nn
            else :
                break

    times = []
    runTimes = []
    
    pEqnIters = 0
    regimeMapTime = 0
    regimeModelsTime = 0
    alphaEqnsTime = 0
    pEqnTime = 0
    EEqnsTime = 0
    nSubCycles = 0
    maxCo = 0
    maxUrCo = 0

    pEqnIterss = []
    regimeMapTimes = []
    regimeModelsTimes = []
    alphaEqnsTimes = []
    pEqnTimes = []
    EEqnsTimes = []
    nSubCycless = []
    maxCos = []
    
    initialResiduals = []
    
    readResidual = False
    firstTime = True
    
    for line in lines :
        if "Courant Number (" in line:
            maxCo = float(line.split()[6])
        elif "Max Ur" in line :
            maxUrCo = float(line.split()[5])
        elif "Time = " in line :
            splitLine = line.split()
            if len(splitLine) == 3 :
                time = float(splitLine[2])
                if not firstTime :
                    times.append(time)
                    runTime = regimeMapTime + regimeModelsTime + alphaEqnsTime + pEqnTime + EEqnsTime
                    runTimes.append(runTime)
                    
                    maxCos.append(max(maxCo, maxUrCo))
                    maxCo = 0
                    maxUrCo = 0
                    pEqnIterss.append(pEqnIters)
                    pEqnIters = 0
                    regimeMapTimes.append(regimeMapTime) 
                    regimeMapTime = 0
                    regimeModelsTimes.append(regimeModelsTime)
                    regimeModelsTime = 0
                    alphaEqnsTimes.append(alphaEqnsTime)
                    alphaEqnsTime = 0
                    pEqnTimes.append(pEqnTime)
                    pEqnTime = 0
                    EEqnsTimes.append(EEqnsTime)
                    EEqnsTime = 0
                    nSubCycless.append(nSubCycles)
                    nSubCycles = 0
            else :
                firstTime = False
        elif "Outer iteration "+str(n) in line :
            readResidual = True
        elif "GAMG:  Solving for p_rgh" in line :
            splitLine = line.split()
            pEqnIters += int(splitLine[14])
            if readResidual :
                residual = float(splitLine[7].split(",")[0])
                initialResiduals.append(residual)
                readResidual = False
        elif "regime map" in line :
            regimeMapTime += float(line.split()[5])
        elif "regime models" in line :
            regimeModelsTime += float(line.split()[5])
        elif "alphaEqns took" in line :
            alphaEqnsTime += float(line.split()[3])
        elif "pEqn-UEqns" in line :
            pEqnTime += float(line.split()[3])
        elif "EEqnsTime" in line :
            EEqnsTime += float(line.split()[3])
        elif "No SubCycles" in line :
            nSubCycles = float(line.split()[7])
    
    ax1.semilogy(limit(times, initialResiduals), initialResiduals, label=filename)
    fRunTimes = filter(runTimes, 10)
    #ax2.plot(limit(times, fRunTimes), fRunTimes, label=filename)
    
ax1.legend()
ax1.grid(True)
plt.show()
