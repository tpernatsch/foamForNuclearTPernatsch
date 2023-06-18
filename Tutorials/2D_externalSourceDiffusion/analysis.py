"""
Script to analyse the source amplification behavior of the extended diffusion 
solver of GeN-Foam. The total neutron produced S+nuSigmaF*flux must be equal 
to S/(1-keff).

Author: Thomas Guilbaud, EPFL/Transmutex SA, 18/06/2023

Source:
[1] Augusto GANDINI & Massimo SALVATORES (2002) The Physics of Subcritical
    Multiplying Systems, Journal of Nuclear Science and Technology, 39:6, 
    673-686, DOI:10.1080/18811248.2002.9715249
"""

import sys

if (len(sys.argv) != 2):
    print(f"Usage: python3 {sys.argv[0]} path/to/case/timeStep")
    sys.exit(1)

#=============================================================================*
# Usefull function

def extractFromFunctionObjs(filename: str, pattern: str):
    with open(filename, 'r') as file:
        for line in file.readlines():
            if (pattern in line):
                return(float(line.split()[1].split(';')[0]))
    return(None)


#=============================================================================*
# Nuclear data

nuSigmaF = [1, 3]
SigmaDisp = [5, 4]
SigmaS = [
    [10,  5],
    [ 0, 10]
]
chiP = [1, 0]
chiD = [0.5, 0.5]
beta = [7.2315e-05, 0.000609661, 0.000471181, 0.00118907, 0.000445487, 9.58515e-05]


#=============================================================================*
# Extract results

filename = "steadyState/100/uniform/functionObjects/functionObjectProperties"

# Flux integrated over the mesh (neutrons.m/s)
flux0 = extractFromFunctionObjs(filename, "volIntegrate(zone0,flux0)")
flux1 = extractFromFunctionObjs(filename, "volIntegrate(zone0,flux1)")
flux = [
    extractFromFunctionObjs(filename, f"volIntegrate(zone0,flux{i})") for i in range(2)
]

# Source integrated over the mesh (neutrons/s)
S = extractFromFunctionObjs(filename, "volIntegrate(zone0,externalSourceFlux0)")

# Computed data
betat = sum(beta)
chi = [
    (1-betat)*chiP[0] + betat*chiD[0],
    (1-betat)*chiP[1] + betat*chiD[1]
]
nuSigmaFfluxt = sum(
    [nuS*ph for nuS, ph in zip(nuSigmaF, flux)]
)

# Steady-state equations for an infinite medium
eq0 = chi[0]*nuSigmaFfluxt - SigmaDisp[0]*flux[0] + SigmaS[1][0]*flux[1] + S
eq1 = chi[1]*nuSigmaFfluxt - SigmaDisp[1]*flux[1] + SigmaS[0][1]*flux[0] + S


#=============================================================================*
# Test results

print(f"Flux eq 0 = {eq0} must be close to 0")
print(f"Flux eq 1 = {eq1} must be close to 0")

NprodFiss = chi[0]*nuSigmaFfluxt + chi[1]*nuSigmaFfluxt
NprodScatt = SigmaS[1][0]*flux[1] + SigmaS[0][1]*flux[0]
Nabs = SigmaDisp[0]*flux[0] + SigmaDisp[1]*flux[1]
Nsrc = 2*S
Nprod = NprodFiss+NprodScatt

# Equations from [1], ksrc and keff must be equivalent
ksrc = NprodFiss / (NprodFiss + Nsrc)
keff = (NprodFiss) / (Nabs-NprodScatt)

# Match the keff only if eigenvalueNeutronics at true and externalSourceNeutronics at false
# According to [1], keff and ksrc must close in external source mode.
# keffEigenvalue = 0.949712

print( "| Parameter  | Value    |")
print( "|:-----------|:--------:|")
print(f"| ksrc       | {ksrc:.6f} |")
print(f"| keff       | {keff:.6f} |")
# print(f"| Eigenvalue | {keffEigenvalue:.6f} |")

mult = 1/(1-keff)
print(f"Source amplification = 1/(1-keff) = {mult}")

NsrcMult = Nsrc*mult
Ntot = NprodFiss+Nsrc

print("Total neutrons produced in [neutrons/s] integrated over the mesh:")
print(f"S/(1-keff)      = {NsrcMult}")
print(f"S+nuSigmaF*flux = {Ntot}")

relErr = abs(Ntot-NsrcMult)/((Ntot+NsrcMult)/2)
print(f"Relative error = {relErr}")

if (relErr > 1e-6):
    sys.exit(1)

sys.exit(0)

#=============================================================================*
