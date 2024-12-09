"""
Script to extract the main parameters of an NTP rocket engine using the last 
time in a log.GeN-Foam file.

Author: Thomas Guilbaud, EPFL, 03/12/2022
"""

#==============================================================================*
# Imports
#==============================================================================*

from numpy import pi, sqrt
import os
import sys

#==============================================================================*
# Units and Constants
#==============================================================================*

m, s, K, kg, mol = 1, 1, 1, 1, 1
inch = 2.54e-2*m
g = 1e-3*kg
m2 = m**2
N = kg*m/s**2
kN = 1e3*N
J = kg*(m*s)**2
W = J/s
MW = 1e6*W
Pa = N/m2
MPa = 1e6*Pa

gConst = 9.8 * m/(s**2)
Rconst = 8.314 * J/mol/K
gamma_H2 = 1.4
M_H2 = 2.014 * g/mol


#==============================================================================*
# Reactor Parameters
#==============================================================================*

nFuelAssembly = 1500/6 # Number of fuel assembly in core
nominalPower = 937*MW/nFuelAssembly
fuelElementStructureFraction = 0.696911
SfuelElement = 0.00190587*m2 * (1-fuelElementStructureFraction) # Fluid Sp
SPlenum = 0.00222351*m2

outletNozzlePres = 0 * MPa
nozzleMinArea = pi * (8.720 * inch/m /2)**2
nozzleOutletArea = pi * (30.370 * inch/m /2)**2


#==============================================================================*
# Usefull Functions and Classes
#==============================================================================*

def computeNozzleExhaustVelocity(
    T_chamber: float, P_chamber: float, P_outletNozzle: float, gamma: float,
    M_gas: float, V_chamber: float=0
) -> float:
    """
    Compute the nozzle exhaust velocity in m/s.
        T_chamber: chamber temperature.
        P_chamber: chamber pressure.
        P_outletNozzle: nozzle outlet pressure.
        gamma: gamma constant of the gas.
        M_gas: molar mass of the propellant.
    """
    try:
        gm1og = (gamma-1)/gamma
        fracP = P_outletNozzle / P_chamber
        factor = 2/gm1og * Rconst*T_chamber/M_gas
        return(sqrt(factor * (1 - fracP**gm1og) + V_chamber))
    except Exception as e:
        return(0)


def computeIsp(nozzleExhaustVelocity: float) -> float:
    """ Compute the specific impulse in s. """
    return(nozzleExhaustVelocity / gConst)


def computeThrust(massFlow: float, nozzleExhaustVelocity: float) -> float:
    """ Compute the thrust in N. """
    return(massFlow * nozzleExhaustVelocity)


def relErr(th: float, val: float) -> float:
    """ Return the relative error in %. """
    return(100 * abs(th-val)/th)


class SuperDict:
    def __init__(self) -> None:
        self.dict = {}

    def __getitem__(self, key):
        try:
            return(self.dict[key])
        except:
            return(0)


    def __setitem__(self, key, value) -> None:
        self.dict[key] = value

    def any(self) -> bool:
        return(any(self.dict))

    def computeDensity(self) -> None:
        keys = self.dict.keys()
        if ('alphaRhoPhi' in keys and 'alphaPhi' in keys):
            self['rho'] = self['alphaRhoPhi'] / self['alphaPhi']

    def findInLine(self, line: str, pattern: str, key: str=None, pos: int=4) -> None:
        if (pattern in line):
            if (not key):
                key = pattern.split()[-1]
            self[key] = float(line.split()[pos])


#==============================================================================*
# User information on usage
#==============================================================================*

if (len(sys.argv) != 2):
    print(f"\n  Usage: python3 {sys.argv[0]} path/to/log.GeN-Foam\n")
    sys.exit(0)


#==============================================================================*
# Data extraction
#==============================================================================*

# Open the log file
filename = sys.argv[1]
logfile = open(filename, 'r')

# Data
time, executionTime, clockTime = 0, 0, 0
powers = SuperDict()
inlet,        outlet        = SuperDict(), SuperDict()
inletFuel,    outletFuel    = SuperDict(), SuperDict()
inletCentral, outletCentral = SuperDict(), SuperDict()

# Read the file in reverse
for line in logfile.readlines()[::-1]:
    powers.findInLine(line, "volIntegrate(fluidRegion) of powerDensityNeutronics", "fluid")
    powers.findInLine(line, "volIntegrate(neutroRegion) of powerDensity", "neutro")
    powers.findInLine(line, "volIntegrate(fluidRegion) of hdeltaT", "AhdeltaT")


    inlet.findInLine(line, "patch inlet massFlow")
    inlet.findInLine(line, "patch inlet massFlow", 'S', pos=7)
    inlet.findInLine(line, "areaAverage(inlet) of magU", 'U')
    inlet.findInLine(line, "areaAverage(inlet) of p")
    inlet.findInLine(line, "areaAverage(inlet) of T")
    inlet.findInLine(line, "areaAverage(inlet) of alphaRhoPhi")
    inlet.findInLine(line, "areaAverage(inlet) of alphaPhi")
    inlet.findInLine(line, "patch inlet TBulk")
    inlet.findInLine(line, "min(T) = ", 'Tmin', pos=2)

    outlet.findInLine(line, "patch outlet massFlow")
    outlet.findInLine(line, "patch outlet massFlow", 'S', pos=7)
    outlet.findInLine(line, "areaAverage(outlet) of magU", 'U')
    outlet.findInLine(line, "areaAverage(outlet) of p")
    outlet.findInLine(line, "areaAverage(outlet) of T")
    outlet.findInLine(line, "areaAverage(outlet) of alphaRhoPhi")
    outlet.findInLine(line, "areaAverage(outlet) of alphaPhi")
    outlet.findInLine(line, "patch outlet TBulk")
    outlet.findInLine(line, "max(T) = ", 'Tmax', pos=2)


    inletFuel.findInLine(line, "faceZone inletFuelElement massFlow")
    inletFuel.findInLine(line, "faceZone inletFuelElement massFlow", 'S', pos=7)
    inletFuel.findInLine(line, "areaAverage(inletFuelElement) of magU", 'U')
    inletFuel.findInLine(line, "areaAverage(inletFuelElement) of p")
    inletFuel.findInLine(line, "areaAverage(inletFuelElement) of T")
    inletFuel.findInLine(line, "areaAverage(inletFuelElement) of alphaRhoPhi")
    inletFuel.findInLine(line, "areaAverage(inletFuelElement) of alphaPhi")
    inletFuel.findInLine(line, "faceZone inletFuelElement TBulk")

    outletFuel.findInLine(line, "faceZone outletFuelElement massFlow")
    outletFuel.findInLine(line, "faceZone outletFuelElement massFlow", 'S', pos=7)
    outletFuel.findInLine(line, "areaAverage(outletFuelElement) of magU", 'U')
    outletFuel.findInLine(line, "areaAverage(outletFuelElement) of p")
    outletFuel.findInLine(line, "areaAverage(outletFuelElement) of T")
    outletFuel.findInLine(line, "areaAverage(outletFuelElement) of alphaRhoPhi")
    outletFuel.findInLine(line, "areaAverage(outletFuelElement) of alphaPhi")
    outletFuel.findInLine(line, "faceZone outletFuelElement TBulk")


    inletCentral.findInLine(line, "faceZone inletCentralUnloaded massFlow")
    inletCentral.findInLine(line, "faceZone inletCentralUnloaded massFlow", 'S', pos=7)
    inletCentral.findInLine(line, "areaAverage(inletCentralUnloaded) of magU", 'U')
    inletCentral.findInLine(line, "areaAverage(inletCentralUnloaded) of p")
    inletCentral.findInLine(line, "areaAverage(inletCentralUnloaded) of T")
    inletCentral.findInLine(line, "areaAverage(inletCentralUnloaded) of alphaRhoPhi")
    inletCentral.findInLine(line, "areaAverage(inletCentralUnloaded) of alphaPhi")
    inletCentral.findInLine(line, "faceZone inletCentralUnloaded TBulk")

    outletCentral.findInLine(line, "faceZone outletCentralUnloaded massFlow")
    outletCentral.findInLine(line, "faceZone outletCentralUnloaded massFlow", 'S', pos=7)
    outletCentral.findInLine(line, "areaAverage(outletCentralUnloaded) of magU", 'U')
    outletCentral.findInLine(line, "areaAverage(outletCentralUnloaded) of p")
    outletCentral.findInLine(line, "areaAverage(outletCentralUnloaded) of T")
    outletCentral.findInLine(line, "areaAverage(outletCentralUnloaded) of alphaRhoPhi")
    outletCentral.findInLine(line, "areaAverage(outletCentralUnloaded) of alphaPhi")
    outletCentral.findInLine(line, "faceZone outletCentralUnloaded TBulk")


    if (powers.any() and inlet.any() and outlet.any() and "Time =" == line[:6]):
        time = float(line.split()[2])
        break
    if (powers.any() and inlet.any() and outlet.any() and "ExecutionTime" == line[:13]):
        executionTime = float(line.split()[2])
        clockTime = float(line.split()[6])


#==============================================================================*
# Compute parameters
#==============================================================================*

print(f"Time = {time} s   Execution Time = {executionTime} s   Clock Time = {clockTime} s")

# Compute the densities
inlet.computeDensity()
outlet.computeDensity()
inletFuel.computeDensity()
outletFuel.computeDensity()
inletCentral.computeDensity()
outletCentral.computeDensity()


# Compute nozzle exhaust velocity
nozzleExhaustVelocity = computeNozzleExhaustVelocity(
    outlet['TBulk'], outlet['p'], outletNozzlePres, gamma_H2, M_H2
)

idealNozzleExhaustVelocity = computeNozzleExhaustVelocity(
    outlet['Tmax'], outlet['p'], outletNozzlePres, gamma_H2, M_H2
)

# Compute the Isp
IspNozzle = computeIsp(nozzleExhaustVelocity)
idealIspNozzle = computeIsp(idealNozzleExhaustVelocity)

# Compute the thrust
thrustNozzle = computeThrust(nFuelAssembly * inlet['massFlow'], nozzleExhaustVelocity)
idealThrustNozzle = computeThrust(nFuelAssembly * inlet['massFlow'], idealNozzleExhaustVelocity)

effectiveExhaustVelocity, characteristicVelocity = 0, 0
try:
    effectiveExhaustVelocity = nozzleExhaustVelocity + outletNozzlePres * nozzleOutletArea/(nFuelAssembly*inlet['massFlow'])
    characteristicVelocity = outlet['p'] * nozzleMinArea / (nFuelAssembly*inlet['massFlow'])
except:
    print("inlet['massFlow'] is 0 ! Maybe the simulation didn't print a result in the log.")
    sys.exit(1)


#==============================================================================*
# Print results
#==============================================================================*
# The expected values are taken from Finseth 1991 Final Report, Table C-1: Reactor Test Summary

print(f"""
|                             | Unit | Ref.  | GeN-Foam | Error [%] |
|:----------------------------|:----:|:-----:|:--------:|:---------:|
| Core power                  | MW   |   937 | {nFuelAssembly*powers['neutro']/MW:8.0f} | {relErr(937, nFuelAssembly*powers['neutro']/MW):9.2f} |
| Core inlet mass flowrate    | kg/s |    31 | {nFuelAssembly*inlet['massFlow']:8.1f} | {relErr(31, nFuelAssembly*inlet['massFlow']):9.2f} |
| Core plenum  inlet pressure | MPa  | 3.861 | {inlet['p']/MPa:8.3f} | {relErr(3.861, inlet['p']/MPa):9.2f} |
| Core plenum outlet pressure | MPa  | 3.427 | {outlet['p']/MPa:8.3f} | {relErr(3.427, outlet['p']/MPa):9.2f} |
| Core plenum  inlet temperat | K    |  83.3 | {inlet['TBulk']:8.1f} | {relErr(83.3, inlet['TBulk']):9.2f} |
| Core plenum outlet temperat | K    |  1972 | {outlet['TBulk']:8.1f} | {relErr(1972, outlet['TBulk']):9.2f} |
| Exhaust velocity            | m/s  |  7549 | {nozzleExhaustVelocity/(m/s):8.0f} | {relErr(7549, nozzleExhaustVelocity/(m/s)):9.2f} |
| Thrust                      | kN   | 234.0 | {thrustNozzle/kN:8.1f} | {relErr(234.0, thrustNozzle/kN):9.2f} |
| Isp  (outlet plenum)        | s    |   770 | {IspNozzle/s:8.1f} | {relErr(770, IspNozzle/s):9.2f} |
""")

print(f"""
Power:
    Expected          = {nominalPower      :8.6e} W
    Neutro integrated = {powers['neutro']  :8.6e} W ({powers['neutro']  /nominalPower*100:7.3f} %)
    Fluid  integrated = {powers['fluid']   :8.6e} W ({powers['fluid']   /nominalPower*100:7.3f} %)
    int Ah(Tf-Tco) dV = {powers['AhdeltaT']:8.6e} W ({powers['AhdeltaT']/nominalPower*100:7.3f} %)
Temperature:
    Plenums          : {inlet['TBulk']       :.1f} -> {outlet['TBulk']       :6.1f} K  (83.3 -> 1972 K expected ({relErr(1972, outlet['TBulk']):.2f} %))
    Fuel element     : {inletFuel['TBulk']   :.1f} -> {outletFuel['TBulk']   :6.1f} K
    Central unloaded : {inletCentral['TBulk']:.1f} -> {outletCentral['TBulk']:6.1f} K
    Extrema          : {inlet['Tmin']:.1f} -> {outlet['Tmax']:6.1f} K  (83.3 -> 2356 K expected ({relErr(2356, outlet['Tmax']):.2f} %))
Velocity:
    Plenums          : {inlet['U']       :5.2f} -> {outlet['U']       :6.2f} m/s
    Fuel element     : {inletFuel['U']   :5.2f} -> {outletFuel['U']   :6.2f} m/s
    Central unloaded : {inletCentral['U']:5.2f} -> {outletCentral['U']:6.2f} m/s
Density:
    Plenums          : {inlet['rho']       :5.2f} -> {outlet['rho']       :5.2f} kg/m3
    Fuel element     : {inletFuel['rho']   :5.2f} -> {outletFuel['rho']   :5.2f} kg/m3
    Central unloaded : {inletCentral['rho']:5.2f} -> {outletCentral['rho']:5.2f} kg/m3
Pressure:
    Plenums          : {inlet['p']/MPa       :.3f} -> {outlet['p']/MPa       :.3f} MPa  (3.861 -> 3.427 MPa expected)
    Fuel element     : {inletFuel['p']/MPa   :.3f} -> {outletFuel['p']/MPa   :.3f} MPa
    Central unloaded : {inletCentral['p']/MPa:.3f} -> {outletCentral['p']/MPa:.3f} MPa
Mass flowrate:
    Fuel element     : {inletFuel['massFlow']   :6.3f} -> {outletFuel['massFlow']   :6.3f} kg/s
    Central unloaded : {inletCentral['massFlow']:6.3f} -> {outletCentral['massFlow']:6.3f} kg/s
    Plenums          : {inlet['massFlow']       :6.3f} -> {outlet['massFlow']       :6.3f} kg/s
    Plenums (core)   : {nFuelAssembly*inlet['massFlow']:6.3f} -> {nFuelAssembly*outlet['massFlow']:6.3f} kg/s  (31.0 kg/s expected ({relErr(31.0, nFuelAssembly*outlet['massFlow']):.2f} %)))
From mass conservation:
    rho_i S_i u_i = rho_o S_o u_o => {
        inlet['rho']  * inlet['S']  * inlet['U']:.3f} = {
        outlet['rho'] * outlet['S'] * outlet['U']:.3f} kg/s
From Perfect Gas law:
    r = P_i/(rho_i T_i) = P_o/(rho_o T_o) => {
        inlet['p'] /(inlet['rho'] *inlet['TBulk']) :.3f} = {
        outlet['p']/(outlet['rho']*outlet['TBulk']):.3f} m²/s²/K
Rocketery from Rocket propulsion elements 8th edition:
    Exhaust velocity: {nozzleExhaustVelocity/(m/s):.2f} m/s (v_2) (7549 m/s expected ({relErr(7549, nozzleExhaustVelocity/(m/s)):.2f} %))
    Ideal Exhaust v : {idealNozzleExhaustVelocity/(m/s):.2f} m/s (v_2)
    c               : {effectiveExhaustVelocity/(m/s):.2f} m/s (Effective exhaust velocity (Isp*g))
    c*              : {characteristicVelocity/(m/s):.2f} m/s (Characteristic velocity)
    thrust          : {thrustNozzle/N:.1f} N (m*c) (234000 N expected ({relErr(234000, thrustNozzle/N):.2f} %))
    thrust (ideal)  : {idealThrustNozzle/N:.1f} N (m*c)
    Isp             : {IspNozzle/s:.1f} s (v2/g) (770.3 s expected ({relErr(770.3, IspNozzle/s):.2f} %))
    Isp (ideal Tmax): {idealIspNozzle/s:.1f} s (v2/g)  (834 s expected ({relErr(834, idealIspNozzle/s):.2f} %))
""")

#==============================================================================*
