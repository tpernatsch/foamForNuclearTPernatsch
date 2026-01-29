#!/usr/bin/env python3

from matplotlib import pyplot
import pandas
import numpy as np

data1 = pandas.read_fwf('postProcessing/samples/0/T', skiprows=1)
data2 = pandas.read_fwf('postProcessing/samples/0/Css', skiprows=1)
data3 = pandas.read_fwf('postProcessing/samples/0/Cpp', skiprows=1)

t = data1['# Time']/3600
T = data1['0']
Css = data2['0']
Cpp = data3['0']

Tmax = round(T.max())
Tmin = round(T.min())
steps = np.r_[210, 23310, 65310, 296310]/3600

indices = [round(i) for i in np.interp(steps, t, np.arange(len(t)))]

pyplot.figure(figsize=(5,4))
for (i0, i1, label) in zip(indices[:-1], indices[1:], ['100 K/min', '10 K/min', '1 K/min', '0.1 K/min']):
    pyplot.plot(T[i0:i1], Css[i0:i1], label=label)

# Overlay plots of Tssd and Tssp
T_ = np.linspace(Tmin-50, Tmax+50, 100)
e = 1.6022e-19
NA = 6.0221e23
R = 8.31446
Tssp0 = 30853.0
Tssd0 = 1.02e5
Qp = 0.26*e*NA
Qd = 0.37*e*NA
Tssp_ = lambda T_: Tssp0 * np.exp(-Qp / (R*T_))
Tssd_ = lambda T_: Tssd0 * np.exp(-Qd / (R*T_))
Tssp = Tssp_(T_)
Tssd = Tssd_(T_)

pyplot.plot(T_, Tssd, 'k', label=None, lw=0.5, ls=(0, (8, 8)))
pyplot.plot(T_, Tssp, 'k', label=None, lw=0.5, ls=(0, (8, 8)))

ymin = 0
ymax = 1.2*Css.max()
xmin = Tmin-20
xmax = Tmax+20
pyplot.plot([Tmin, Tmin], [ymin, ymax], 'k', label=None, lw=0.5, ls=(0, (8, 8)))
pyplot.plot([Tmax, Tmax], [ymin, ymax], 'k', label=None, lw=0.5, ls=(0, (8, 8)))

pyplot.legend(framealpha=1)
pyplot.xlim(xmin=xmin, xmax=xmax)
pyplot.ylim(ymin=ymin, ymax=ymax)
pyplot.xlabel('Temperature (K)')
pyplot.ylabel('Css (wt.ppm)')
ax = pyplot.gca()
ax.annotate('$T_{min}$', xy=(Tmin, 40), xytext=(Tmin+35, 50), ha='left', va='bottom', 
            arrowprops=dict(arrowstyle="->", lw=0.5),
            bbox=dict(pad=1, facecolor="none", edgecolor="none"))
ax.annotate('$T_{max}$', xy=(Tmax, 10), xytext=(Tmax-35, 5), ha='right', va='bottom', 
            arrowprops=dict(arrowstyle="->", lw=0.5),
            bbox=dict(pad=1, facecolor="none", edgecolor="none"))
ax.annotate('$T_{SSP}$', xy=(380, Tssp_(380)), xytext=(380-35, Tssp_(380)+5), 
            ha='center', va='bottom', 
            arrowprops=dict(arrowstyle="->", lw=0.5),
            bbox=dict(pad=1, facecolor="none", edgecolor="none"))
ax.annotate('$T_{SSD}$', xy=(550, Tssd_(550)), xytext=(550+45, Tssd_(550)-10), 
            ha='center', va='top', 
            arrowprops=dict(arrowstyle="->", lw=0.5), 
            bbox=dict(pad=1, facecolor="none", edgecolor="none"))

pyplot.tight_layout()
pyplot.savefig('Css.png', dpi=300)

#pyplot.gca().set_xscale('function', functions=(lambda x: 1/x, lambda x: 1/x))
ax.set_yscale('function', functions=(lambda x: np.power(x,0.5), lambda x: np.power(x,2)))
ax.yaxis.get_ticklocs(minor=True)
ax.minorticks_on()
pyplot.tight_layout()
pyplot.savefig('Css_log.png', dpi=300)

pyplot.figure(figsize=(5,4))
pyplot.plot(t, Css, label='Css')
pyplot.plot(t, Cpp, label='Cpp')
pyplot.plot(t, Cpp+Css, label='Ctot')
pyplot.xlabel('Time (h)')
pyplot.ylabel('Hydrogen Concentration (wt.ppm)')
pyplot.legend()
pyplot.tight_layout()
pyplot.savefig('timePlot_Css.png', dpi=300)

pyplot.figure(figsize=(5,4))
pyplot.plot(t, T, label='T')
pyplot.xlabel('Time (h)')
pyplot.ylabel('Temperature (K)')
pyplot.legend()
pyplot.tight_layout()
pyplot.savefig('timePlot_T.png', dpi=300)
