#!/usr/bin/env python3

from matplotlib import pyplot
import numpy as np
import os
import pandas

def is_float(s):
    try:
        v = float(s)
    except:
        return False
    
    return True

# Plots of concentration profiles
root = 'postProcessing/graphs'
subdirs = [f.name for f in os.scandir(root) if f.is_dir() and is_float(f.name)]
timeNames = sorted(subdirs, key=lambda x: float(x))

for timeName in timeNames:
    print(f'Creating plot for t={timeName}s')
    x, T, Css, Cpp, Ctot, KG, KN, KD, TSSp, TSSd = np.loadtxt(f'postProcessing/graphs/{timeName}/line_T_Css_Cpp_Ctot_KG_KN_KD_TSSp_TSSd.xy').T
    x *= 1000

    pyplot.plot(x, Cpp, label='$C_{pp}$')
    pyplot.plot(x, Css, label='$C_{ss}$')
    pyplot.plot(x, Ctot, label='$C_{tot}$')
    pyplot.plot(x, TSSp, 'k--', label='$T_{ssp}$')
    pyplot.plot(x, TSSd, 'k', label='$T_{ssd}$')
    pyplot.ylabel('C (wt.ppm)')
    pyplot.xlabel('x (mm)')
    pyplot.legend()
    pyplot.tight_layout()
    pyplot.savefig(f'postProcessing/graphs/C_{timeName}.png', dpi=300)
    pyplot.close()


# Time-dependent plot of volume-average contrations
data1 = pandas.read_csv('postProcessing/average/0/volFieldValue.dat', skiprows=3, sep='\t')
t = data1['# Time        '].to_numpy()
Css = data1['volAverage(Css)'].to_numpy()
Cpp = data1['volAverage(Cpp)'].to_numpy()

pyplot.plot(t, Css, label='Css')
pyplot.plot(t, Cpp, label='Cpp')
pyplot.plot(t, Css+Cpp, label='Ctot')
pyplot.ylabel('C (wt.ppm)')
pyplot.xlabel('Time (day)')
pyplot.legend()
pyplot.tight_layout()
pyplot.savefig(f'timePlot.png', dpi=300)
pyplot.close()

# Check final concentration value
C_ref = 200.0
Css_ref = 183.8244
Cpp_ref = 200 - Css_ref
Css_final = Css[-1]*1000/t[-1]
Cpp_final = Cpp[-1]*1000/t[-1]
C_final = Css_final + Cpp_final

print(f'Final Css: {Css_final:.5f} wt.ppm')
print(f'Final Cpp: {Cpp_final:.5f} wt.ppm')
print(f'Final total hydrogen concentration: {C_final:.5f} wt.ppm')
print(f'Expected average value: {C_ref} wt.ppm')
print(f'Average Error: {100*abs(C_final/C_ref - 1):g}%')

if abs(C_final/C_ref - 1) > 1e-3:
    raise RuntimeError(f'Final total concentration should equal {C_ref} wt.ppm')

if abs(Css_final/Css_ref - 1) > 1e-3:
    raise RuntimeError(f'Final Css should equal {Css_ref} wt.ppm')

if abs(Cpp_final/Cpp_ref - 1) > 1e-3:
    raise RuntimeError(f'Final Cpp should equal {Cpp_ref} wt.ppm')

print("Test successful")

