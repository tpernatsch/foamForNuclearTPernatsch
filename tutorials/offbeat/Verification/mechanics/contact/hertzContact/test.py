# Test the contact forces
import numpy as np
from matplotlib import pyplot

d1 = 1e-1
d2 = 1e12
L = 1.11803e-3
P = 1e7
F0 = P*d1*L
    
E = 50e9
nu = 0.3

b = np.sqrt(2*F0/np.pi/L*((1-nu**2)/E + (1-nu**2)/E)/(1/d1 + 1/d2))
pMax = 2*F0/np.pi/b/L

y = 5e-2*np.linspace(0, 1, 100)**2
sigY = pMax/np.sqrt(1+y**2/b**2)

filename = "postProcessing/xyGraph/1/base_sigma.xy"
data = np.loadtxt(filename).T
y1 = data[0]
sigma1 = data[1:]
sigY1 = sigma1[3]
e1 = (np.interp(y1, y, sigY)+sigY1) / sigY.max()

filename = "postProcessing/xyGraph/1/cylinder_sigma.xy"
data = np.loadtxt(filename).T
y2 = data[0]
sigma2 = data[1:]
sigY2 = sigma2[3]
e2 = (np.interp(y2, y, sigY)+sigY2) / sigY.max()


pyplot.figure(layout="constrained")
pyplot.plot(y, sigY/1e6, "k", label="analytical")
pyplot.plot(y1, -sigma1[3]/1e6, "x", label="base")
pyplot.plot(y2, -sigma2[3]/1e6, "+", label="cylinder")
pyplot.legend()
pyplot.ylabel(r"$\sigma_{yy}$ (MPa)")
pyplot.xlabel("y (m)")
pyplot.savefig("sigma.png")

assert np.max(np.abs(e1)) < 0.15
assert np.average(np.abs(e1)) < 0.05
assert np.max(np.abs(e2)) < 0.15
assert np.average(np.abs(e2)) < 0.05

print("Test Succeeded")
