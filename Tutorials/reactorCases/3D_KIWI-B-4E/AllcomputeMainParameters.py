"""
Script to analyse the last time in a log.GeN-Foam file.

Author: Thomas Guilbaud, 03/12/2022
"""

#==============================================================================*
# Imports
#==============================================================================*

from numpy import pi, sqrt
import sys
from unitsSI import *


#==============================================================================*
# Constants
#==============================================================================*

gConst = 9.8 * m/(s**2)
Rconst = 8.314 * kg*m2/s**2 / mol/K # J/mol/K
gamma_H2 = 1.4
M_H2 = 2.014 * g/mol


#==============================================================================*
# Reactor Parameters
#==============================================================================*

nFuelAssembly = 1 # Number of fuel assembly in core 1500/6
fuelElementStructureFraction = 0.696911
# SfuelElement = 0.00190587*m2 * (1-fuelElementStructureFraction) # Fluid Sp
# SPlenum = 0.00222351*m2

outletNozzlePres = 0 * MPa
nozzleMinArea = pi * (8.720 * inch/m /2)**2
nozzleOutletArea = pi * (30.370 * inch/m /2)**2

# Theoretical values from LA-3185-MS page 110
nominalPower = 905 * MW
massFlowTh = 31.083 * kg/s
inletTemperatureTh = 93.9 * K
outletTemperatureTh = 1956 * K
inletPressureTh = 3.951 * MPa
outletPressureTh = 3.372 * MPa


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
            self['massFlow'] = self['U'] * self['rho'] * self['S']

    def findInLine(self, line: str, pattern: str, key: str=None, pos: int=4) -> None:
        if (pattern in line):
            if (not key):
                key = pattern.split()[-1]
            self[key] = float(line.split()[pos])


#==============================================================================*

# User information on usage
if (len(sys.argv) != 2):
    print(f"\n    Usage: python3 {sys.argv[0]} path/to/log.GeN-Foam\n")
    sys.exit(0)

#==============================================================================*

# Open the log file
filename = sys.argv[1]
logfile = open(filename, 'r')

# Data
time, executionTime, clockTime = 0, 0, 0
powers = SuperDict()
inlet,        outlet        = SuperDict(), SuperDict()
inletFuel,    outletFuel    = SuperDict(), SuperDict()

#inletFuel,    outletFuel    = SuperDict(), SuperDict()
inletCentral, outletCentral = SuperDict(), SuperDict()

# Read the file in reverse to extract the last results
for line in logfile.readlines()[::-1]:
    powers.findInLine(line, "volIntegrate(fluidRegion) of powerDensityNeutronics", "fluid")
    powers.findInLine(line, "volIntegrate(neutroRegion) of powerDensity", "neutro")
    powers.findInLine(line, "volIntegrate(fluidRegion) of hdeltaT", "AhdeltaT")

    inlet.findInLine(line, "areaAverage(inlet) of magU", 'U')
    inlet.findInLine(line, "areaAverage(inlet) of p")
    inlet.findInLine(line, "areaAverage(inlet) of T")
    inlet.findInLine(line, "areaAverage(inlet) of alphaRhoPhi")
    inlet.findInLine(line, "areaAverage(inlet) of alphaPhi")
    inlet.findInLine(line, "min(T) = ", 'Tmin', pos=2)

    outlet.findInLine(line, "areaAverage(outlet) of magU", 'U')
    outlet.findInLine(line, "areaAverage(outlet) of p")
    outlet.findInLine(line, "areaAverage(outlet) of T")
    outlet.findInLine(line, "areaAverage(outlet) of alphaRhoPhi")
    outlet.findInLine(line, "areaAverage(outlet) of alphaPhi")
    outlet.findInLine(line, "max(T) = ", 'Tmax', pos=2)
    outlet.findInLine(line, "max(Tsurface.lumpedNuclearStructure) = ", 'Tsurfacemax', pos=2)
    outlet.findInLine(line, "max(Tmatrix.lumpedNuclearStructure) = ", 'Tmatrixmax', pos=2)
    outlet.findInLine(line, "max(Tmax.lumpedNuclearStructure) = ", 'Tfuelmax', pos=2)
    outlet.findInLine(line, "max(powerDensityNeutronics) = ", 'powerDensityMax', pos=2)

    inletFuel.findInLine(line, "areaAverage(inletFuelElement) of magU", 'U')
    inletFuel.findInLine(line, "areaAverage(inletFuelElement) of p")
    inletFuel.findInLine(line, "areaAverage(inletFuelElement) of T")
    inletFuel.findInLine(line, "areaAverage(inletFuelElement) of alphaRhoPhi")
    inletFuel.findInLine(line, "areaAverage(inletFuelElement) of alphaPhi")

    outletFuel.findInLine(line, "areaAverage(outletFuelElement) of magU", 'U')
    outletFuel.findInLine(line, "areaAverage(outletFuelElement) of p")
    outletFuel.findInLine(line, "areaAverage(outletFuelElement) of T")
    outletFuel.findInLine(line, "areaAverage(outletFuelElement) of alphaRhoPhi")
    outletFuel.findInLine(line, "areaAverage(outletFuelElement) of alphaPhi")

    inletCentral.findInLine(line, "areaAverage(inletCentralUnloaded) of magU", 'U')
    inletCentral.findInLine(line, "areaAverage(inletCentralUnloaded) of p")
    inletCentral.findInLine(line, "areaAverage(inletCentralUnloaded) of T")
    inletCentral.findInLine(line, "areaAverage(inletCentralUnloaded) of alphaRhoPhi")
    inletCentral.findInLine(line, "areaAverage(inletCentralUnloaded) of alphaPhi")

    outletCentral.findInLine(line, "areaAverage(outletCentralUnloaded) of magU", 'U')
    outletCentral.findInLine(line, "areaAverage(outletCentralUnloaded) of p")
    outletCentral.findInLine(line, "areaAverage(outletCentralUnloaded) of T")
    outletCentral.findInLine(line, "areaAverage(outletCentralUnloaded) of alphaRhoPhi")
    outletCentral.findInLine(line, "areaAverage(outletCentralUnloaded) of alphaPhi")


    if (powers.any() and inlet.any() and outlet.any() and "Time =" == line[:6]):
        time = float(line.split()[2])
        break
    if (powers.any() and inlet.any() and outlet.any() and "ExecutionTime" == line[:13]):
        executionTime = float(line.split()[2])
        clockTime = float(line.split()[6])

print(f"Time = {time} s   Execution Time = {executionTime} s   Clock Time = {clockTime} s")

inlet['S'] = 0.56602329 # from log file
outlet['S'] = inlet['S']
# inletFuel['S'] = 0.47142417662379
# outletFuel['S'] = inletFuel['S']
# inletCentral['S'] = 0.081399241627666
# outletCentral['S'] = inletCentral['S']


# Compute the densities
inlet.computeDensity()
outlet.computeDensity()
inletFuel.computeDensity()
outletFuel.computeDensity()
inletCentral.computeDensity()
outletCentral.computeDensity()


# Compute nozzle exhaust velocity
nozzleExhaustVelocity = computeNozzleExhaustVelocity(
    outlet['T'], outlet['p'], outletNozzlePres, gamma_H2, M_H2
)

idealNozzleExhaustVelocity = computeNozzleExhaustVelocity(
    outlet['Tmax'], outlet['p'], outletNozzlePres, gamma_H2, M_H2
)

# Compute the Isp
IspNozzle = computeIsp(nozzleExhaustVelocity)
idealIspNozzle = computeIsp(idealNozzleExhaustVelocity)

# Compute the thrust
thrustNozzle = computeThrust(inlet['massFlow'], nozzleExhaustVelocity)
idealThrustNozzle = computeThrust(inlet['massFlow'], idealNozzleExhaustVelocity)

# effectiveExhaustVelocity = IspCore * gConst
effectiveExhaustVelocity = nozzleExhaustVelocity + outletNozzlePres * nozzleOutletArea/inlet['massFlow']

# thrustCore2 = effectiveExhaustVelocity2*inletMFlow[-1]

characteristicVelocity = outlet['p'] * nozzleMinArea / inlet['massFlow']

# Expected         : 83.3 -> 1833.3 K (exit fuel: 2222 K)
# The expected values are taken from Finseth 1991 Final Report, Table C-1: Reactor Test Summary
def relErr(th: float, val: float) -> float:
    return(100 * abs(th-val)/th)

# Theoretical values
nozzleExhaustVelocityTh = computeNozzleExhaustVelocity(T_chamber=outletTemperatureTh, P_chamber=outletPressureTh, P_outletNozzle=outletNozzlePres, gamma=gamma_H2, M_gas=M_H2)
IspNozzleTh = computeIsp(nozzleExhaustVelocityTh)
thrustNozzleTh = computeThrust(massFlowTh, nozzleExhaustVelocityTh)

print(f"""Power:
    Expected          = {nominalPower      :8.6e} W
    Neutro integrated = {powers['neutro']  :8.6e} W ({relErr(nominalPower, powers['neutro']):7.3f} %)
    Fluid  integrated = {powers['fluid']   :8.6e} W ({relErr(nominalPower, powers['fluid']):7.3f} %)
    int Ah(Tf-Tco) dV = {powers['AhdeltaT']:8.6e} W ({relErr(nominalPower, powers['AhdeltaT']):7.3f} %)
Temperature:
    Plenums          : {inlet['T']       :.1f} -> {outlet['T']       :6.1f} K  ({inletTemperatureTh:.1f} -> {outletTemperatureTh:.1f} K expected ({relErr(outletTemperatureTh, outlet['T']):.2f} %))
    Fuel element     : {inletFuel['T']   :.1f} -> {outletFuel['T']   :6.1f} K
    Central unloaded : {inletCentral['T']:.1f} -> {outletCentral['T']:6.1f} K
    Extrema          : {inlet['Tmin']:.1f} -> {outlet['Tmax']:6.1f} K  ({inletTemperatureTh:.1f} -> 2356 K expected ({relErr(2356, outlet['Tmax']):.2f} %))
Velocity:
    Plenums          : {inlet['U']       :5.2f} -> {outlet['U']       :6.2f} m/s
    Fuel element     : {inletFuel['U']   :5.2f} -> {outletFuel['U']   :6.2f} m/s
    Central unloaded : {inletCentral['U']:5.2f} -> {outletCentral['U']:6.2f} m/s
Density:
    Plenums          : {inlet['rho']       :5.2f} -> {outlet['rho']       :5.2f} kg/m3
    Fuel element     : {inletFuel['rho']   :5.2f} -> {outletFuel['rho']   :5.2f} kg/m3
    Central unloaded : {inletCentral['rho']:5.2f} -> {outletCentral['rho']:5.2f} kg/m3
Pressure:
    Plenums          : {inlet['p']/MPa       :.3f} -> {outlet['p']/MPa       :.3f} MPa  ({inletPressureTh/MPa:.3f} -> {outletPressureTh/MPa:.3f} MPa expected)
    Fuel element     : {inletFuel['p']/MPa   :.3f} -> {outletFuel['p']/MPa   :.3f} MPa
    Central unloaded : {inletCentral['p']/MPa:.3f} -> {outletCentral['p']/MPa:.3f} MPa
Mass flowrate:
    Fuel element     : {inletFuel['massFlow']   :6.3f} -> {outletFuel['massFlow']   :6.3f} kg/s
    Central unloaded : {inletCentral['massFlow']:6.3f} -> {outletCentral['massFlow']:6.3f} kg/s
    Plenums          : {inlet['massFlow']       :6.3f} -> {outlet['massFlow']       :6.3f} kg/s  ({massFlowTh:.1f} kg/s expected ({relErr(massFlowTh, inlet['massFlow']):.2f} %)))
From mass conservation:
    rho_i S_i u_i = rho_o S_o u_o => {
        inlet['rho']  * inlet['S']  * inlet['U']:.3f} = {
        outlet['rho'] * outlet['S'] * outlet['U']:.3f} kg/s
From Perfect Gas law:
    r = P_i/(rho_i T_i) = P_o/(rho_o T_o) => {
        inlet['p'] /(inlet['rho'] *inlet['T']) :.3f} = {
        outlet['p']/(outlet['rho']*outlet['T']):.3f} m²/s²/K
Rocketery from Rocket propulsion elements 8th edition:
    Exhaust velocity: {nozzleExhaustVelocity/(m/s):.2f} m/s (v_2) ({nozzleExhaustVelocityTh/(m/s):.0f} m/s expected ({relErr(nozzleExhaustVelocityTh/(m/s), nozzleExhaustVelocity/(m/s)):.2f} %))
    Ideal Exhaust v : {idealNozzleExhaustVelocity/(m/s):.2f} m/s (v_2)
    c               : {effectiveExhaustVelocity/(m/s):.2f} m/s (Effective exhaust velocity (Isp*g))
    c*              : {characteristicVelocity/(m/s):.2f} m/s (Characteristic velocity)
    thrust          : {thrustNozzle/N:.1f} N (m*c) ({thrustNozzleTh/N:.0f} N expected ({relErr(thrustNozzleTh/N, thrustNozzle/N):.2f} %))
    thrust (ideal)  : {idealThrustNozzle/N:.1f} N (m*c)
    Isp             : {IspNozzle/s:.1f} s (v2/g) ({IspNozzleTh/s:.1f} s expected ({relErr(IspNozzleTh/s, IspNozzle/s):.2f} %))
    Isp (ideal Tmax): {idealIspNozzle/s:.1f} s (v2/g)  (834 s expected ({relErr(834, idealIspNozzle/s):.2f} %))
""")


print(f"""Summary:

| Parameter                       | Unit | GeN-Foam | Reference | Diff [%] |
|:--------------------------------|:----:|:--------:|:---------:|:--------:|
| Power neutronics                | MW   | {powers['neutro']/MW:8.1f} | {nominalPower/MW:9.1f} | {relErr(nominalPower/MW, powers['neutro']/MW):8.2f} |
| Power $\int Ah(T_f-T_{{co}})dV$   | MW   | {powers['AhdeltaT']/MW:8.1f} | {nominalPower/MW:9.1f} | {relErr(nominalPower/MW, powers['AhdeltaT']/MW):8.2f} |
| Max power density               | MW/m3| {outlet['powerDensityMax'] * fuelElementStructureFraction/1e9:8.2f} | {2.90:9.2f} | {relErr(2.90, outlet['powerDensityMax'] * fuelElementStructureFraction/1e9):8.2f} |
| H$_2$ Inlet temperature         | K    | {inlet['T']:8.1f} | {inletTemperatureTh:9.1f} | {relErr(inletTemperatureTh, inlet['T']):8.2f} |
| H$_2$ Nozzle chamber temperature| K    | {outlet['T']:8.1f} | {outletTemperatureTh:9.1f} | {relErr(outletTemperatureTh, outlet['T']):8.2f} |
| Max H$_2$ temperature           | K    | {outlet['Tmax']:8.1f} | {2311:9} | {relErr(2311, outlet['Tmax']):8.2f} |
| Max surface temperature         | K    | {outlet['Tsurfacemax']:8.1f} | {'-':9} | {'-':8} |
| Max fuel temperature            | K    | {outlet['Tfuelmax']:8.1f} | {'-':9} | {'-':8} |
| Inlet pressure                  | MPa  | {inlet['p']/MPa:8.3f} | {inletPressureTh/MPa:9.3f} | {relErr(inletPressureTh/MPa, inlet['p']/MPa):8.2f} |
| Nozzle chamber pressure         | MPa  | {outlet['p']/MPa:8.3f} | {outletPressureTh/MPa:9} | {relErr(outletPressureTh/MPa, outlet['p']/MPa):8.2f} |
| Mass flow rate                  | kg/s | {inlet['massFlow']:8.1f} | {massFlowTh:9.1f} | {relErr(massFlowTh, inlet['massFlow']):8.2f} |
| Exhaust velocity*               | m/s  | {nozzleExhaustVelocity/(m/s):8.1f} | {nozzleExhaustVelocityTh/(m/s):9.1f} | {relErr(nozzleExhaustVelocityTh/(m/s), nozzleExhaustVelocity/(m/s)):8.2f} |
| Thrust*                         | kN   | {thrustNozzle/N/1e3:8.1f} | {thrustNozzleTh/N/1e3:9.1f} | {relErr(thrustNozzleTh/N, thrustNozzle/N):8.2f} |
| Specific impulse \isp*          | s    | {IspNozzle/s:8.1f} | {IspNozzleTh/s:9.1f} | {relErr(IspNozzleTh/s, IspNozzle/s):8.2f} |

""")

#==============================================================================*
