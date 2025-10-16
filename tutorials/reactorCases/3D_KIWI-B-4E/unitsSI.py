# Author: Thomas Guilbaud

import numpy as np

# --- Constants
pi = np.pi

# --- Basis
kg, s, m, J, mol, K, rad = 1, 1, 1, 1, 1, 1, 1

# --- Mass
g = 1e-3 * kg
mg = 1e-3 * g
pounds = 0.453592 * kg

# --- Time
minutes = 60 * s
hour = 60 * minutes
day = 24 * hour

# --- Distance
cm = 1e-2 * m
mm = 1e-3 * m
inch = 2.54 * cm

# --- Area
cm2 = cm**2
m2 = m**2

# --- Volume
cm3 = cm**3
m3 = m**3

# --- Angle
deg = np.pi/180 * rad

# --- Energy
eV = 1.602e-19 * J
µeV = 1e-6 * eV
meV = 1e-3 * eV
keV = 1e3 * eV
MeV = 1e6 * eV

# --- Power
W =  J/s
kW = 1e3 * W
MW = 1e6 * W

# --- Force
N = kg * m * s**2

# --- Pressure
Pa = N/m2
MPa = 1e6 * Pa
psi = 0.00689457 * MPa
bar = 1e5 * Pa
psia = psi

# --- Percent
percent = 1e-2
pcm = 1e-5

# --- Temperature
degR = 0.555555 * K # Rankine
