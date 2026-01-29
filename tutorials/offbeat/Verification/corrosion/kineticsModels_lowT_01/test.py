# Verification for low temperature oxidation models
# Compare against model of Almarshad, Oregon State University, 1990
import numpy as np
from matplotlib import pyplot
import re

def logFileIter(filename):
    buffer = open(filename, 'r').read()
    for m in re.finditer(r"Time = (\S+) (?:seconds|days)$(.*?)Execution Time",
                         buffer, re.DOTALL | re.MULTILINE):
        t = float(m.group(1))
        buffer2 = m.group(2)
        matches = re.findall(r'Interface temperature: (\S+)$', buffer2, re.MULTILINE)
        Ti = float(matches[-1])
        matches = re.findall(r'Average oxide layer thickness: (\S+)$', buffer2, re.MULTILINE)
        dOxide = float(matches[-1])
        yield (t, Ti, dOxide)

refValues = {
    'MATPRO-CORROS-PWR': {
        0: (620.0, 1.65733e-05),
        50: (624.174, 1.78768e-05),
        80: (627.022, 1.8796e-05),
        100: (629.097, 1.94809e-05),
        150: (635.073, 2.15224e-05)
    },
    'MATPRO-CORROS-BWR': {
        0: (620.0, 3.18747e-05),
        50: (628.678, 3.71712e-05),
        80: (635.565, 4.16757e-05),
        100: (641.292, 4.56144e-05),
        150: (663.841, 6.26609e-05),
    },
    'EPRI-KWU-CE': {
        0: (620.0, 2.49062e-05),
        50: (626.615, 2.83581e-05),
        80: (631.594, 3.10695e-05),
        100: (635.52, 3.32759e-05),
        150: (648.778, 4.11597e-05)
    },
    'EPRI-SLI': {
        0: (620.0, 2.0886e-05),
        50: (625.412, 2.31784e-05),
        80: (629.29, 2.48698e-05),
        100: (632.219, 2.61686e-05),
        150: (641.217, 3.02994e-05)
    }
}

success = True

for model in ['MATPRO-CORROS-PWR', 'MATPRO-CORROS-BWR', 'EPRI-KWU-CE', 
              'EPRI-SLI']:

    # Plot of interface temperature versus oxide thickness
    pyplot.figure(layout='constrained', figsize=(5,4))
    pyplot.suptitle(f'{model}, T=620 K, Time=1000 days')
    
    for Q in [0, 50, 80, 100, 150]:
        t,Ti,dOxide = np.array(list(logFileIter(f'{model}_Q{Q}.out'))).T
        pyplot.plot(1e6*dOxide, Ti, label=f"$q''={Q}\\ W/cm^2$")
        if not ((np.abs(Ti[-1] - refValues[model][Q][0]) < 1e-2) and 
               (np.abs(dOxide[-1]/refValues[model][Q][1] - 1) < 1e-4)):
            success = False
            print(f'Code regression error for Q={Q} W/cm2')
    
    pyplot.legend()
    pyplot.xlabel('Oxide Thickness $(\\mu m)$')
    pyplot.ylabel('Oxide-Metal Interface Temperature (K)')
    pyplot.xlim(xmin=0)
    pyplot.savefig(f'heatFluxInfluence_{model}.png', dpi=300)
    
    
    # Time-dependent plot of interface temperature
    pyplot.figure(layout='constrained', figsize=(5,4))
    pyplot.suptitle(f'{model}, T=620 K')
    
    for Q in [0, 50, 80, 100, 150]:
        t,Ti,dOxide = np.array(list(logFileIter(f'{model}_Q{Q}.out'))).T
        pyplot.plot(t, Ti, label=f"$q''={Q}\\ W/cm^2$")
    
    pyplot.ylabel('Oxide-Metal Interface Temperature (K)')
    pyplot.xlabel('Time (days)')
    pyplot.xlim(xmin=0)
    pyplot.legend()
    pyplot.savefig(f"interfaceT_{model}.png", dpi=300)
    
    
    # Time-dependent plot of oxide thickness
    pyplot.figure(layout='constrained', figsize=(5,4))
    pyplot.suptitle(f'{model}, T=620 K')
    for Q in [0, 50, 80, 100, 150]:
        t,Ti,dOxide = np.array(list(logFileIter(f'{model}_Q{Q}.out'))).T
        r = 1e6*np.diff(dOxide)/np.diff(t)
        T = 0.5*(Ti[:-1]+Ti[1:])
        pyplot.plot(t[1:][r<50], 1e6*dOxide[1:][r<50], label=f"$q''={Q}\\ W/cm^2$")
    
    pyplot.xlim(xmin=0)
    pyplot.ylim(ymin=0)
    pyplot.ylabel(r'Oxide Thickness ($\mu$m)')
    pyplot.xlabel('Time (days)')
    pyplot.legend()
    pyplot.savefig(f"dOxide_{model}.png", dpi=300)

if success:
    print("Test successful")
else:
    data = {}
    for model in ['MATPRO-CORROS-PWR', 'MATPRO-CORROS-BWR', 'EPRI-KWU-CE', 
                  'EPRI-SLI']:
        data[model] = {}
        for Q in [0, 50, 80, 100, 150]:
            t,Ti,dOxide = np.array(list(logFileIter(f'{model}_Q{Q}.out'))).T
            data[model][Q] = (Ti[-1], dOxide[-1])
    
    print(data)
    
    raise RuntimeError("Test failed")
