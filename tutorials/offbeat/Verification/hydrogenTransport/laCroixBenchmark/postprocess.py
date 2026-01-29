#!/usr/bin/env python3

from matplotlib import pyplot
import pandas
import numpy as np

data1 = pandas.read_fwf('postProcessing/samples/0/T', skiprows=1)
data2 = pandas.read_fwf('postProcessing/samples/0/Css', skiprows=1)
data3 = pandas.read_fwf('postProcessing/samples/0/Cpp', skiprows=1)

expDataC1 = pandas.read_excel('experimentalData.xlsx', sheet_name='cycle1')
expDataC2 = pandas.read_excel('experimentalData.xlsx', sheet_name='cycle2')
expDataC3 = pandas.read_excel('experimentalData.xlsx', sheet_name='cycle3')

keyPoints = {
    "0":    (100, 25, 0),
    "1":    (1600, 425, 240),
    "2":    (3400, 170, 5),
    "3":    (6200, 384, 181),
    "3'":   (5900, 285, 130),
    "4":    (7900, 285, 100),
    "5":   (9100, 425, 260),
    "6":    (10300, 170, 30)
}


t = data1['# Time']
T = data1['0'] - 273.15
Css = data2['0']
Cpp = data3['0']

steps = [0, 4444, 7710, 12000]
indices = [round(i) for i in np.interp(steps, t, np.arange(len(t)))]


pyplot.figure(figsize=(4,3))
for (i0, i1, label, s, w) in zip(indices[:-1], indices[1:], ['A', 'B', 'C'], ['g-', 'r--', 'b-.'], [1.25,1,0.75]):
    pyplot.plot(T[i0:i1], Css[i0:i1], s, lw=w, label=label)

pyplot.plot(expDataC1['T (C)'], expDataC1['Css (ppm)'], 'g^', mfc='none', alpha=0.7, label='Exp. A', markersize=3, markeredgewidth=0.5)
pyplot.plot(expDataC2['T (C)'], expDataC2['Css (ppm)'], 'r>', mfc='none', alpha=0.7, label='Exp. B', markersize=3, markeredgewidth=0.5)
pyplot.plot(expDataC3['T (C)'], expDataC3['Css (ppm)'], 'bo', mfc='none', alpha=0.7, label='Exp. C', markersize=3, markeredgewidth=0.5)

for label, (t_,T_,Css_) in keyPoints.items():
    pyplot.annotate(label, (T_,Css_))

pyplot.legend(loc=2)
pyplot.xlabel('Temperature ($^o C$)')
pyplot.ylabel('Css (wt.ppm)')
pyplot.ylim(ymin=0, ymax=300)
pyplot.xlim(xmin=150, xmax=450)
pyplot.tight_layout()
pyplot.savefig('Css.png', dpi=300)

pyplot.figure(figsize=(4,3))
pyplot.plot(t, Css, label='Css')
pyplot.plot(t, Cpp, label='Cpp')
pyplot.plot(t, Cpp+Css, label='Ctot')
pyplot.xlim(xmin=0, xmax=12000)
pyplot.xlabel('Time (s)')
pyplot.ylabel('Hydrogen Concentration (wt.ppm)')
pyplot.legend(loc=4)
pyplot.tight_layout()
pyplot.savefig('timePlot.png', dpi=300)


pyplot.figure(figsize=(4,3))
for (i0, i1, label, s, w) in zip(indices[:-1], indices[1:], ['A', 'B', 'C'], ['g-', 'r--', 'b-.'], [1, 1, 1]):
    pyplot.plot(t[i0:i1], T[i0:i1], s, lw=w, label=label)

for label, (t_,T_,Css_) in keyPoints.items():
    pyplot.annotate(label, (t_, T_))

pyplot.xlim(xmin=0, xmax=12000)
pyplot.ylim(ymax=500)
pyplot.xlabel('Time (s)')
pyplot.ylabel('Temperature ($^o C$)')
pyplot.legend(loc=8)
pyplot.tight_layout()
pyplot.savefig('T.png', dpi=300)
