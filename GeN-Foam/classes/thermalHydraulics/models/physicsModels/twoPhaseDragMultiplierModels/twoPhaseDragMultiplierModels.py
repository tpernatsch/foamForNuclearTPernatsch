import math
from matplotlib import pyplot as plt

### PHYSICAL PROPERTIES OF FLUIDS

rhol = 850
rhov = 0.5
mul = 1e-4
muv = 1e-7

### CALCULATIONS SETTINGS

# When plotting against vapour quality
xN = 500
xRange = [x/xN for x in range(xN+1)]

# When plotting against the sqrt of the Lockhart-Martinelli parameter
XLMMax = 30
XLMMin = 4.51999e-2
XLMN = 100
XLMRange = [XLMMin+(dXLM/XLMN) for dXLM in range(XLMMax*XLMN)]

### FUNCTIONS

def XLM(x) :

	x = max(x, 1e-6)
	XLM = (mul/muv)**0.1*((1-x)/x)**0.9*math.sqrt(rhov/rhol)
	return max(min(XLM, XLMMax), XLMMin)

# Kottowski-Savatteri

def KottowskiSavatteri_x(x) :

	log10X = math.log10(XLM(x))
	exp = 2*(0.1046*(log10X**2) - 0.5098*log10X + 0.6252)
	return 10**exp

def KottowskiSavatteri_XLM(XLM) :

	log10X = math.log10(XLM)
	exp = 2*(0.1046*(log10X**2) - 0.5098*log10X + 0.6252)
	return 10**exp

# Kaiser (1988)

def Kaiser88_x(x) :

	pass

def Kaiser88_XLM(XLM) :

	pass

def plot(xAxis, func) :

	yAxis = []
	for x in xAxis :
		yAxis.append(func(x))
	plt.semilogy(xAxis, yAxis, label=func.__name__)

### MAIN

plot(XLMRange, KottowskiSavatteri_XLM)

plt.legend()
plt.grid()
plt.show()
