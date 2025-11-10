"""
Script to compute the Stodola coefficient for a turbine.

Author: Thomas Guilbaud, 2023/03/16
"""

#=============================================================================*

import numpy as np

#=============================================================================*

def StodolaCoeff(
    massFlow: float,
    rhoin: float,
    Pin: float,
    Pout: float
) -> float:
    return(
        massFlow/np.sqrt(rhoin*Pin*(1-(Pout/Pin)**2))
    )

#=============================================================================*


if __name__ == '__main__':
    # Use SI units
    massFlow = 192.8 # kg/s
    rhoin = 0.0718e3 # kg/m3
    maxP = 18.8e6 # Pa
    interP = 6e6 # Pa
    minP = 0.008e6 # Pa

    KtHP = StodolaCoeff(massFlow, rhoin, maxP, interP)
    KtLP = StodolaCoeff(massFlow, rhoin, interP, minP)

    print(f"Kt High Pressure = {KtHP}")
    print(f"Kt Low  Pressure = {KtLP}")


    print(StodolaCoeff(31, 11, 5.9e6, 3.861e6))

#=============================================================================*
