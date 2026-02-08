# Test the contact forces
import re
import numpy as np

A = 0.000111803
P0 = 1e7
F0 = P0*A
n0 = [0, 1, 0]

filename = "log.offbeat"
buffer = open(filename, "r").read()
patchData = {}
for patchName in ["top", "bottom", "contact1", "contact2"]:    
    m = re.findall(f"{patchName}.*:\\s+\\((.*)\\)", buffer)[-1]
    F = np.array([float(v) for v in m.strip().split()])
    magF = np.sqrt(F@F)
    n = F / magF
    print(f"{patchName}:\t{F}")
    if patchName in ["top", "contact1"]:
        assert np.abs(1 + n@n0) < 1e-4
        assert np.abs(magF / F0 - 1) < 1e-4
    elif patchName in ["bottom", "contact2"]:
        assert np.abs(1 - n@n0) < 1e-4
        assert np.abs(magF / F0 - 1) < 1e-4

print ("Test Successful")
