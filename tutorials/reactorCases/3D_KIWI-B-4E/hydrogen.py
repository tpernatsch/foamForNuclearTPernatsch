"""
Thermophysical properties of the hydrogen.

Author: Thomas Guilbaud, EPFL 
        Eymeric Simonnot, EPFL

- [x] Kp
- [x] Mole fraction x_H
- [x] Weight fraction
- [~] Pressure
- [x] Density
- [x] Specific Heat capacity
- [x] Specific enthalpy
- [x] Thermal conductivity
- [x] Dynamic viscosity
- [x] Diffusion coefficient
- [x] Compressibility psi (use thermodynamic definition)
- [~] Compression factor Z (to be verified)
- [x] Entropy
- [x] Specific internal energy (use thermodynamic definition)


[1] Fang et al, Study on high-temperature hydrogen dissociation for nuclear thermal 
    propulsion reactor, Nuclear Engineering and Design 392 (2022) 111753

[2] Reisfeld, THERMODYNAMIC AND TRANSPORT PROPERTIES OF DISSOCIATED HYDROGEN MIXTURES
    1957, LANL, LA-2123.

[3] Roder et al, Computer Programs for Thermodynamic and Transport Properties 
    of Hydrogen (Tabcode-ll), 1972, NBS Technical note 625.

[4] NIST Chemistry WebBook, Hydrogen atom - Gas phase thermochemestry data : 
    https://webbook.nist.gov/cgi/cbook.cgi?ID=C12385136&Mask=1

[5] Muzny, C.D., Huber, M.L., Kazakov, A.F., 2013. Correlation for the viscosity
    of normal hydrogen obtained from symbolic regression. J. Chem. Eng. Data 58 
    (4), 969–979. https://doi.org/10.1021/je301273j

[6] Vanderslice, J.T., Weissman, S., Mason, E.A., Fallon, R.J., 1962. High-temperature
    transport properties of dissociating hydrogen. Phys. Fluids 5 (2), 155–164. 
    https://doi.org/10.1063/1.1706590

[7] NIST Chemistry WebBook, Hydrogen molecule - Gas phase thermochemestry data : 
    https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=1

[8] Assael, M.J., Assael, J.-A.- M., Huber, M.L., Perkins, R.A., Takata, Y., 
    2011. Correlation of the Thermal Conductivity of Normal and Parahydrogen 
    from the Triple Point to 1000 K and up to 100 MPa. J. Phys. Chem. Ref. 
    Data 40 (3), 033101. 
    https://doi.org/10.1063/1.3606499

[9] Leachman, J. W., R. T Jacobsen, S. G. Penoncello, et E. W. Lemmon. 
    « Fundamental Equations of State for Parahydrogen, Normal Hydrogen, and Orthohydrogen ». 
    Journal of Physical and Chemical Reference Data 38, no 3 (1 septembre 2009): 721‑48. 
    https://doi.org/10.1063/1.3160306.

[10] Robert D.McCarty, Hydrogen Technological Survey, Thermophysical Properties
    Cryogenics Division, Institute for Basic Standards
    https://ntrs.nasa.gov/api/citations/19760004209/downloads/19760004209.pdf

[11] Mason, Edward A., et William E. Rice. « The Intermolecular Potentials of Helium and Hydrogen ».
     The Journal of Chemical Physics 22, no 3 (1 mars 1954): 522‑35.
     https://doi.org/10.1063/1.1740100.

[12] Mason, Edward A. « Transport Properties of Gases Obeying a Modified Buckingham (Exp-Six) Potential ».
     The Journal of Chemical Physics 22, no 2 (1 février 1954): 169‑86.
     https://doi.org/10.1063/1.1740026.

[13] Amdur, I., et E. A. Mason. « Properties of Gases at Very High Temperatures ».
     Physics of Fluids 1, no 5 (1958): 370. https://doi.org/10.1063/1.1724353.

"""

#=============================================================================*
# Imports

from numpy import sqrt, cos, arccos, exp, log, log10, pi, linspace, cbrt, array
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

#=============================================================================*
# Units

kg, m3, K, Pa = 1, 1, 1, 1
degR = 0.555556 * K
poundPerCubFoot = 16.0185 * kg/m3
psia = 6894.76 * Pa
atm = 101325 * Pa
MPa = 1e6 * Pa

#=============================================================================*
# Physical constants

# Convert units to OpenFOAM, 
W = lambda: 2.01588e-3 # [kg/mol]
W_H_ = 1.007647e-3 # [kg/mol]
M_12 = 2.0*W_H_*W()/(W_H_+W())
omega = lambda: -0.219
Tc = lambda: 33.145 * K
Pc = lambda: 1.2965e6 * Pa
Vc = lambda: 6.4481e-5 # [m3/mol]
R_ = 8.314472 # J/mol/K



def rms(l1, l2):
    return(
        sqrt(sum(
            [(e1-e2)**2 for e1, e2 in zip(l1, l2)]
        ))
    )

def sig(x, a, b):
    return(1/(1+exp(-a*(x-b))))

def func(x, a, b, c, d, e, f, g, h, i):
    sig_ = sig(x, h, i)
    return (a*x**2 + b*x + c) * sig_ + (1-sig_) * (d*sig(x, e, g) + f)

def computeRsquare(yData, yPredict):
    residuals = yData - yPredict
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((yData - np.mean(yData))**2)
    return(1 - (ss_res / ss_tot))


class PerfectGas:
    def __init__(self) -> None:
        pass

    def rho(self, T, p):
        return(p/(R_*T)*W())
    
    def H(self, T, p):
        return(0)
    
    def Cp(self, T, p):
        CpCoeff = [
            1.266445630e+04, 1.17094171e+01, -2.98179826e-02, 3.81902598e-05, 
            -2.54178633e-08, 9.28566400e-12, -1.77480852e-15, 1.39051574e-19
        ]
        return(sum(Cp_*pow(T, i) for i, Cp_ in enumerate(CpCoeff)))
    
    def mu(self, T, p):
        return(8.9385e-6)

    def kappa(self, T, p):
        return(0.01238494+0.00057219*T)

#=============================================================================*

def Kp(T: float) -> float:
    """ 
    Equilibrium constant Kp from [1] Eq 3 and [2].

    T: Temperature in K
    """
    return(
        10.0**(-2.37943e4/T + 6.33153)
    )


def x_H(T: float, p: float) -> float:
    """
    Mole fraction of atomic hydrogen x_H from [1] Eq 2 and [2].

    T: Temperature in K
    p: Pressure in Pa (originally in atm)
    """
    p = p/atm
    Kp_ = Kp(T)
    if (Kp_ == 0):
        return(0)
    return(
        2.0 / (1.0 + sqrt(1.0 + 4.0*p/Kp_))
    )


def w_H(T: float, p: float) -> float:
    """
    Weight fraction of atomic hydrogen w_H from [1] Eq 7'.

    T: Temperature in K
    p: Pressure in Pa
    """
    x_H_ = x_H(T, p)
    return(
        x_H_*W_H_ / (x_H_*W_H_ + (1.0-x_H_)*W())
    )


#=============================================================================*
# DEFINITION OF THE AVERAGE EFFECTIVE COLLISION INTEGRALS

def omega11_11(T: float) -> float:
    """
    Ω(1,1){1,1} defined in [1]
    Coeffs fitted with data from [6], table II
    Fitted on values from 1000 to 10000 K every 500 K and extended
    at cryogenic temperatures

    T: Temperature in K
    """

    if(T>1000):
        return(
            -1.173958E-19*pow(T, 5.0)
            +3.983904E-15*pow(T, 4.0)
            -5.443976E-11*pow(T, 3.0)
            +3.918410E-07*pow(T, 2.0)
            -1.704143E-03*T
            +6.552751
        )
    else:
        return 44.13*pow(T, -0.331)

    # return(
    #     -1.173958E-19*pow(T, 5.0)
    #     +3.983904E-15*pow(T, 4.0)
    #     -5.443976E-11*pow(T, 3.0)
    #     +3.918410E-07*pow(T, 2.0)
    #     -1.704143E-03*T
    #     +6.552751
    #     )

def omega22_11(T: float) -> float:
    """
    Ω(2,2){1,1} defined in [1]
    Coeffs fitted with data from [6], table II, for T>1000K
    Coeffs fitted with data from [12], using parameters from [6],
    [11], and Eq.(23) from [13] and extended to cryogenic temperatures

    T: Temperature in K
    """

    if(T>1000):
        return(
            5.570961E-20*pow(T, 5.0)
            -1.444695E-15*pow(T, 4.0)
            +1.024643E-11*pow(T, 3.0)
            +3.257067E-08*pow(T, 2.0)
            -8.379299E-04*T
            +6.786617
        )

    else:
        return 37.879*pow(T, -0.268)
    # return(
    #     5.570961E-20*pow(T, 5.0)
    #     -1.444695E-15*pow(T, 4.0)
    #     +1.024643E-11*pow(T, 3.0)
    #     +3.257067E-08*pow(T, 2.0)
    #     -8.379299E-04*T
    #     +6.786617
    #     )

def omega11_22(T: float) -> float:
    """ (not used)
    Ω(1,1){2,2} defined in [1]
    Coeffs fitted with data from [6], table II

    T: Temperature in K
    """
    return(5e-16*pow(T, 4.0) - 2e-11*pow(T, 3.0) + 2e-7*pow(T, 2.0) - 0.0015*T + 6.401)

def omega22_22(T: float) -> float:
    """
    Ω(2,2){2,2} defined in [1]
    Coeffs fitted with data from [6], table II

    T: Temperature in K
    """
    return(5e-16*pow(T, 4.0) - 2e-11*pow(T, 3.0) + 2e-7*pow(T, 2.0) - 0.0015*T + 7.1214);

def omega11_12(T: float) -> float:
    """
    Ω(1,1){1,2} defined in [1]
    Coeffs fitted with data from [6], table II
    Fitted on values from 1000 to 10000 K every 500 K and extended to 
    cryogenic temperatures

    T: Temperature in K
    """
    return(-1.024*log(T) + 11.054)
    return(
        -2.031857E-19*pow(T, 5.0)
        +6.646410E-15*pow(T, 4.0)
        -8.548538E-11*pow(T, 3.0)
        +5.556305E-07*pow(T, 2.0)
        -2.029612E-03*T
        +5.697181
    )

def omega22_12(T: float) -> float:
    """
    Ω(2,2){1,2} defined in [1]
    Coeffs fitted with data from [6], table II
    Fitted on values from 1000 to 10000 K every 500 K and extended to 
    cryogenic temperatures

    T: Temperature in K
    """
    return(-1.227*log(T) + 13.42)
    return(
        -2.324315E-19*pow(T, 5.0)
        +7.601576E-15*pow(T, 4.0)
        -9.779885E-11*pow(T, 3.0)
        +6.367384E-07*pow(T, 2.0)
        -2.343590E-03*T
        +6.916536
    )


#=============================================================================*

def A_12(T: float) -> float:
    """
    A{1,2} defined in [1]
    Coeffs fitted with data from [6], table II
    Fitted on values from 1000 to 10000 K every 500 K and extended to 
    cryogenic temperatures

    T: Temperature in K
    """
    return(0.0345*log(T) + 0.9917)
    return(
        -1.394617E-17*pow(T, 4.0)
        +4.049174E-13*pow(T, 3.0)
        -4.565141E-09*pow(T, 2.0)
        +2.901797E-05*T
        +1.210780
    )

def B_12(T: float) -> float:
    """
    B{1,2} defined in [1]
    Coeffs fitted with data from [6], table II
    Fitted on values from 1000 to 10000 K every 500 K and extended to 
    cryogenic temperatures

    T: Temperature in K
    """
    return(0.0462*log(T) + 0.8701)
    return(
        -1.483131E-17*pow(T, 4.0)
        +4.234837E-13*pow(T, 3.0)
        -4.759966E-09*pow(T, 2.0)
        +3.241425E-05*T
        +1.172436
    )

def C_12(T: float) -> float:
    """ (not used)
    C{1,2} defined in [1]
    Coeffs fitted with data from [6], table II
    Fitted on values from 1000 to 10000 K every 500 K and extended to 
    cryogenic temperatures

    T: Temperature in K
    """
    return(-0.031*log(T) + 1.1108)
    return(
         1.087058E-17*pow(T, 4.0) 
        -2.989146E-13*pow(T, 3.0) 
        +3.191297E-09*pow(T, 2.0) 
        -2.088554E-05*T 
        +9.090170E-01
    )


#=============================================================================*

def cpPerfectGas(T):
    """
    cp perfect gas of H2
    T: Temperature in K
    """
    A_, B_, C_, D_, E_, F_, G_, H_ = 0, 0, 0, 0, 0, 0, 0, 0

    if (T < 298):
        u_k = [1.616, -0.4117, -0.792, 0.758, 1.217]
        v_k = [531, 751, 1989, 2484, 6859] # 1/K

        sumTerms = sum([
            u_k_i*pow(v_k_i/T, 2.0)*exp(v_k_i/T) / pow(exp(v_k_i/T)-1.0, 2.0) 
            for u_k_i, v_k_i in zip(u_k, v_k)
        ])

        return((2.5*R_ + R_*sumTerms)/W())


    elif (298 <= T and T <= 1000.0):
        A_ = 33.066178
        B_ = -11.363417
        C_ = 11.432816
        D_ = -2.772874
        E_ = -0.158558
        # F_ = -9.980797
        # G_ = 172.707974
        # H_ = 0

    elif (1000.0 <= T and T <= 2500):
        A_ = 18.563083
        B_ = 12.257357
        C_ = -2.859786
        D_ = 0.268238
        E_ = 1.977990
        # F_ = -1.147438
        # G_ = 156.288133
        # H_ = 0
    
    elif (2500 <= T and T <= 6000):
        A_ = 43.413560
        B_ = -4.293079
        C_ = 1.272428
        D_ = -0.096876
        E_ = -20.533862
        # F_ = -38.515158
        # G_ = 162.081354
        # H_ = 0

    T_ = T/1000.0

    cp_Nist = A_ + B_*T_ + C_*pow(T_, 2.0) + D_*pow(T_, 3.0) + E_/pow(T_, 2.0)

    return (cp_Nist/W())


def cvPerfectGas(T):
    """
    cv perfect gas of H2 from 
    cp perfect gas and Mayer's relation for perfect gas
    T: Temperature in K
    """
    return cpPerfectGas(T)-R_/W()
    

def specificHeatCapacityH(T: float) -> float:
    """
    Coefficients from:
    - [4] for T in [298; 6000] K for monoatomic hydrogen heat capacity
    - Note : the equation of State for monoatomic hydrogen is PerfectGas, hence its capacity depends only on Temperature

    T: Temperature in K
    return heat capacity of monoatomic H in J/kg/K 
    """
    if (T < 298): # not super clean
        return(0)

    A_ = 20.78603
    B_ = 4.850638e-10
    C_ = -1.582916e-10
    D_ = 1.525102e-11
    E_ = 3.196347e-11

    T_ = T/1000.0

    Cp_NIST = A_ + B_*T_ + C_*pow(T_, 2.0) + D_*pow(T_, 3.0) + E_/pow(T_, 2.0)

    return(Cp_NIST/W_H_)


def specificHeatCapacityH2(T: float, p: float) -> float:
    """
    Coefficients from:
    - [7] for T in [298; 6000] K for hydrogen heat capacity
    u_k in [9] table 3
    v_k in [9] table 3

    T: Temperature in K
    return heat capacity of H2 in J/kg/K 
    """

    Rmass = R_/W()

    n_ = 0.4986 + 1.1735*omega() + 0.4754*pow(omega(), 2.0)
    alpha0_ = 0.42748*pow(Rmass*Tc(), 2.0) / Pc()
    alpha_ = alpha0_ / pow(T/Tc(), n_)
    b_ = 0.08664*Rmass*Tc() / Pc()
    c_ = Rmass*Tc()/(Pc() + alpha0_/(Vc() * (Vc()+b_))) + b_ - Vc()
    V = 1/densityH2(T,p)

    dPdT = Rmass/(V - b_ + c_) + n_*alpha_/(T*V*(V+b_))
    dPdV = -Rmass*T / pow(V-b_+c_, 2) + alpha_*(b_+ 2*V)/pow(V*(V+b_), 2)
    # dPdT = R_/(V - b_ + c_) 
    # dPdV = -R_*T / sqrt(V-b_+c_) + alpha_/(V*sqrt(V+b_))+ alpha_/((V+b_)*sqrt(V))

    dVdT = - dPdT / dPdV
    # a=0.2476
    # b=0.02661
    # dVdT = Rmass*pow(V,3)*(V-b)/(Rmass*T*pow(V,3)-2*a*pow(V-b,2))

    cp_dep = -p*dVdT + Rmass -n_*alpha_/T*(1+n_)/b_*log((V+b_)/V) - alpha_*(1+n_)*dVdT/(V*(V+b_))    

    return(cpPerfectGas(T) - cp_dep)


def specificHeatCapacity(T: float, p: float) -> float: 
    """
    Mass weighted heat capacity for H-H2 mix
    T: Temperature in K
    p: Pressure in Pa

    heatCapacity: heat capacity of H-H2 mix in J/mol/K 
    return specific heat capacity in J/kg/K 
    """
    # Quick fix to be verified
    # if (T < 298):
    #     CpCoeff = [
    #         1.266445630e+04, 1.17094171e+01, -2.98179826e-02, 3.81902598e-05, 
    #         -2.54178633e-08, 9.28566400e-12, -1.77480852e-15, 1.39051574e-19
    #     ]
    #     return(sum(Cp_*pow(T, i) for i, Cp_ in enumerate(CpCoeff)))

    w_H_ = w_H(T, p)
    
    return(w_H_*specificHeatCapacityH(T) + (1.0-w_H_)*specificHeatCapacityH2(T,p))


def diffusion(T: float, p: float) -> float:
    """
    D_H_H2 defined in [1]

    T: Temperature in K
    p: Pressure in Pa
    return diffusion coefficient of H-H2 in m2/s
    """
    p = p/atm
    return(2.628e-7/(p*omega11_12(T))*sqrt(pow(T, 3.0)/(M_12*1e3)))


def diffusionH2(T: float, p: float) -> float:
    """
    D_H2 defined in [1]

    T: Temperature in K
    """
    p = p/atm
    return(2.628e-7/(p*omega11_11(T))*sqrt(pow(T, 3.0)/(W()*1e3)))


def pressure(T: float, rho_: float) -> float:
    """
    Coefficients from:
    - [3] for T in [100; 400] K and p in [3; 10] MPa
    - [1] for the rest (Aungier-Redlich-Kwong)

    T: Temperature in K
    rho_: density in kg/m3 (originally in mol/L for T in [100; 400] K)
    return Pressure in Pa (orginally in atm for T in [100; 400] K)
    """

    p = 0

    # From [3]
    if (100*K*0 <= T and T <= 400*K or False):
        Rconst = 0.0820535 # L.atm/mol/K
        rho_ = rho_/1000.0/W() # Convert from kg/m3 to mol/L
        N1  =  1.1389049685e-3
        N2  =  1.8072093722e-1
        N3  = -5.3097164419e+1
        N4  =  2.3690885344e+3
        N5  = -1.6306395805e+6
        N6  =  3.2453595439e-8
        N7  = -6.3275640592e-7
        N8  =  2.3115255859e-3
        N9  =  3.4448764044e+1
        N10 = -3.7301781349e+3
        N11 =  1.0789473341e+5
        N12 =  8.0475741674e-7
        N13 =  1.6581404268e-1
        N14 = -1.8158230417e+1
        N15 =  4.9623738040e+2
        N16 =  5.262288563e-8
        delta = -0.0018 # N17
        p = (
            rho_*Rconst*T 
            + pow(rho_, 2.0) * ( N1*T + N2 + N3/T + N4/pow(T, 2.0) + N5/pow(T, 4.0) ) 
            + pow(rho_, 3.0) * ( N6*pow(T, 2.0) + N7*T + N8 ) 
            + pow(rho_, 3.0) * ( N9/pow(T, 2.0) + N10/pow(T, 3.0) + N11/pow(T, 4.0) ) * exp(delta*pow(rho_, 2.0))
            + pow(rho_, 4.0) * N12*T
            + pow(rho_, 5.0) * ( N13/pow(T, 2.0) + N14/pow(T, 3.0) + N15/pow(T, 4.0) ) * exp(delta*pow(rho_, 2.0))
            + pow(rho_, 6.0) * N16
        )
        p = p*atm

    # From [1], ARK model
    if (p <= 3e6*Pa or T >= 400*K) and True:
        # Molecular hydrogen
        n_ = 0.4986 + 1.1735*omega() + 0.4754*pow(omega(), 2.0)
        alpha0_ = 0.42748*pow(R_*Tc(), 2.0) / Pc()
        alpha_ = alpha0_ / pow(T/Tc(), n_)
        b_ = 0.08664*R_*Tc() / Pc()
        c_ = R_*Tc()/(Pc() + alpha0_/(Vc() * (Vc()+b_))) + b_ - Vc()
        V = W()/rho_
        p = R_*T/(V-b_+c_) - alpha_/(pow(V, 2.0)+b_*V)
        # p = R_*T/(W()/rho_-b_+c_) - alpha_/(pow(W()/rho_, 2.0)+b_*W()/rho_)

    return(p)


def density_H(T: float, p: float):
    """
    T: Temperature in K
    p: Pressure in Pa
    return density of H as a perfect gas in kg/m3
    """
    return(p*W_H_/(R_*T))


def densityH2(T: float, p: float):
    """
    Coefficients from:
    - [3] for T in [100; 400] K and p in [3; 10] MPa
    - [1] for the rest (Aungier-Redlich-Kwong)

    Current derivation from ARK model, solve directly the cubic function

    T: Temperature in K
    p: Pressure in Pa
    return density of H2 in kg/m3
    """

    # Molecular hydrogen
    n_ = 0.4986 + 1.1735*omega() + 0.4754*pow(omega(), 2.0)
    alpha0_ = 0.42748*pow(R_*Tc(), 2.0) / Pc()
    alpha_ = alpha0_ / pow(T/Tc(), n_)
    b_ = 0.08664*R_*Tc() / Pc()
    c_ = R_*Tc()/(Pc() + alpha0_/(Vc() * (Vc()+b_))) + b_ - Vc()

    # A*rho_^3 + B*rho_^2 + C*rho_ + D = 0
    # A_ = -alpha_ * (c_-b_)
    # B_ = W()*(R_*T*b_ - alpha_ - p*b_*(c_-b_))
    # C_ = pow(W(), 2.0)*(R_*T - p*c_)
    # D_ = -p * pow(W(), 3.0)

    A_ = alpha_ * (c_-b_)
    B_ = W()*(alpha_ + p*b_*(c_-b_) - R_*T*b_)
    C_ = pow(W(), 2.0)*(p*c_ - R_*T)
    D_ = pow(W(), 3.0) * p

    # Reduced form
    p_ = (3.0*A_*C_ - pow(B_, 2.0))/(3.0*pow(A_, 2.0))
    q_ = (2.0*pow(B_, 3.0) - 9.0*A_*B_*C_ + 27.0*pow(A_, 2.0)*D_)/(27.0*pow(A_, 3.0))

    # https://en.wikipedia.org/wiki/Cubic_equation#Trigonometric_and_hyperbolic_solutions
    discriminant = pow(p_, 3.0)/27.0 + pow(q_, 2.0)/4.0

    if (discriminant < 0):
        tk = lambda k: 2.0*sqrt(-p_/3.0) * cos(arccos(3.0*q_/(2.0*p_)*sqrt(-3.0/p_))/3.0 - 2*pi*k/3)
        rhoH2 = tk(0) - B_/(3.0*A_)

    else:
        u1 = -q_/2 + sqrt(discriminant)
        u2 = -q_/2 - sqrt(discriminant)
        rhoH2 = cbrt(u1) + cbrt(u2) - B_/(3.0*A_)

    return(rhoH2)


def rho(T: float, p: float):
    """ 
    T: Temperature in K
    p: Pressure in Pa
    return density of mix H-H2 in kg/m3
    """
    # Mix
    w_H_ = w_H(T, p)

    return(w_H_*density_H(T, p) + (1.0-w_H_)*densityH2(T, p))


def specificEnthalpyH(T: float) -> float:
    """
    From NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C12385136&Mask=1#Thermo-Gas
    T: Temperature in K for T in [298; 6000] K
    return specific enthalpy of the atomic hydrogen in J/kg like a perfect gas
    """

    A_ = 20.78603
    B_ = 4.850638e-10
    C_ = -1.582916e-10
    D_ = 1.525102e-11
    E_ = 3.196347e-11
    F_ = 211.8020
    G_ = 139.8711
    H_ = 217.9994

    T_ = T/1000.0

    enthalpyH_standard = 218
    enthalpyH = enthalpyH_standard + A_*T_ + B_*pow(T_, 2.0)/2.0 + C_*pow(T_, 3.0)/3.0 + D_*pow(T_, 4.0)/4.0 - E_/T_ + F_ - H_     #kJ/mol

    return(enthalpyH / W_H_ * 1e3)


def specificEnthalpyH2NIST(T: float, p: float) -> float:
    """
    From NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=1#Thermo-Gas
    Not as [1]
    T: Temperature in K
    p: Pressure in Pa
    return specific enthalpy of the molecular hydrogen in J/kg
    """
    A_, B_, C_, D_, E_, F_, G_, H_ = 0, 0, 0, 0, 0, 0, 0, 0

    if (0 <= T and T <= 1000.0):
        A_ = 33.066178
        B_ = -11.363417
        C_ = 11.432816
        D_ = -2.772874
        E_ = -0.158558
        F_ = -9.980797
        # G_ = 172.707974
        H_ = 0

    elif (1000.0 <= T and T <= 2500):
        A_ = 18.563083
        B_ = 12.257357
        C_ = -2.859786
        D_ = 0.268238
        E_ = 1.977990
        F_ = -1.147438
        # G_ = 156.288133
        H_ = 0
    
    elif (2500 <= T and T <= 6000):
        A_ = 43.413560
        B_ = -4.293079
        C_ = 1.272428
        D_ = -0.096876
        E_ = -20.533862
        F_ = -38.515158
        # G_ = 162.081354
        H_ = 0

    T_ = T/1000.0

    # In kJ/mol
    enthalpy_h2_nist = A_*T_ + B_*pow(T_, 2.0)/2.0 + C_*pow(T_, 3.0)/3.0 + D_*pow(T_, 4.0)/4.0 - E_/T_+F_-H_

    return(enthalpy_h2_nist * 1000.0 / W())


def specificEnthalpyH2(T: float, p: float) -> float:
    """
    T: Temperature in K
    p: Pressure in Pa
    return specific enthalpy of the molecular hydrogen in J/kg
    """
    # Molecular hydrogen
    n_ = 0.4986 + 1.1735*omega() + 0.4754*pow(omega(), 2.0)
    alpha0_ = 0.42748*pow(R_*Tc(), 2.0) / Pc() # J*m3/(mol^2)
    alpha_ = alpha0_ / pow(T/Tc(), n_)
    b_ = 0.08664*R_*Tc() / Pc() # m3/mol
    # -> alpha/b = J/mol

    V = W()/densityH2(T, p) # m3/mol

    # u_k = [1.616, -0.4117, -0.792, 0.758, 1.217]
    # v_k = [531, 751, 1989, 2484, 6859]
        
    # sumTerms = sum(u_k[i]*pow(v_k[i]/T , 2)*exp(v_k[i]/T)/pow(exp(v_k[i]/T)-1,2) for i in range(len(u_k)))
    # cpPerfectGas = (2.5*R_ + R_*sumTerms)/W()

    specific_enthalpyH2_pg = cpPerfectGas(T) * (T)        
    specific_enthalpyH2 = specific_enthalpyH2_pg + (p*V - R_*T - alpha_/b_*(1.0+n_)*log((V+b_)/V))/W()
    return(specific_enthalpyH2)


def specificEnthalpy(T: float, p: float) -> float:
    """
    Coefficients from:
    - [4] for T in [298; 6000] K for monoatomic hydrogen enthalpy
    - [1] for molecular hydrogen enthalpy - Eq. (9)
    - Note : the equation of State for monoatomic hydrogen is PerfectGas, hence its enthalpy depends only on Temperature
    
    T: Temperature in K
    p: Pressure in Pa
    return specific enthalpy in J/kg
    """
    specific_enthalpyH = specificEnthalpyH(T)

    specific_enthalpyH2 = specificEnthalpyH2(T, p)

    # Mix
    w_H_ = w_H(T, p)

    return(w_H_*specific_enthalpyH + (1.0-w_H_)*specific_enthalpyH2)


def viscosity0(T : float) -> float:
    """
    Coefficients from:
    - [5] zero_density_limit Eqs. (3) and (4), coefficients ai in Table 2. 

    T: Temperature in K
    viscosity: dynamic viscosity in Pa.s 
    return dynamic viscosity in Pa.s 
    """
    a_list = [2.09630e-1, -4.55274e-1, 1.43602e-1, -3.35325e-2, 2.76981e-3]
    epsilon_kb = 30.41 * K
    sigma = 0.297 # nm

    T_star = T/epsilon_kb

    S_star = exp(sum([a_i * pow(log(T_star), i) for i, a_i in enumerate(a_list)]))

    return(0.021357*sqrt(W()*1000.0*T) / (pow(sigma, 2.0) * S_star) * 1e-6)


def viscosity1(T : float) -> float:
    """
    Coefficients from:
    - [5] the excess contribution Eqs. (5) (6) and (7), coefficients bi in Table 3.

    T: Temperature in K
    viscosity: dynamic viscosity in Pa.s 
    return dynamic viscosity in Pa.s 
    """

    b_list = [-0.1870, 2.4871, 3.7151, -11.0972, 9.0965, -3.8292, 0.5166]

    epsilon_kb = 30.41 * K
    sigma = 0.297 # nm

    T_star = T/epsilon_kb

    B_star = sum([b_i/pow(T_star, i) for i, b_i in enumerate(b_list)])

    B_ = B_star*pow(sigma, 3.0)

    return(B_*viscosity0(T))


def viscosityH(T: float) -> float:
    """
    T: Temperature in K
    return the dynamic viscosity of the atomic hydrogen in Pa.s
    """
    return(2.6693e-6*sqrt(W_H_*1000.0*T)/omega22_22(T))


def dynamicViscosity(T: float, p: float) -> float:
    """
    Coefficients from:
    - [5] Eq. (9), coefficients Ci in Table 4.

    T: Temperature in K
    p: Pressure in Pa
    return dynamic viscosity in Pa.s 
    """
    # Molecular hydrogen 
    C1 = 6.43449673 
    C2 = 4.56334068e-2
    C3 = 2.32797868e-1
    C4 = 9.58326120e-1
    C5 = 1.27941189e-1
    C6 = 3.63576595e-1

    rho_sc = 90.5 * kg/m3
    rhoH2 = rho(T, p)
    rhoR = rhoH2/rho_sc
    Tr = T/Tc()
    viscosityH2_ = viscosity0(T) + viscosity1(T)*rhoH2 + 1e-6*(
        C1*pow(rhoR, 2.0) * exp(C2*Tr + C3/Tr + C4*pow(rhoR, 2.0)/(C5+Tr) + C6*pow(rhoR, 6.0)))

    # Atomic hydrogen 
    viscosityH_ = viscosityH(T)

    # # MIX
    # x_H_ = x_H(T, p)
    # D_H_H2 = diffusion(T, p)
    # A_ = A_12(T)

    # H_11 = pow(1.0-x_H_, 2.0)/viscosityH2_ + 2.0*(1.0-x_H_)*x_H_/(W()+W_H_) * R_*T/(p*D_H_H2) * (1.0 + 3.0*W_H_*A_/(5.0*W()))
    # H_22 = pow(x_H_, 2.0)/viscosityH_ + 2.0*(1.0-x_H_)*x_H_/(W()+W_H_) * R_*T/(p*D_H_H2) * (1.0 + 3.0*W()*A_/(5.0*W_H_))
    # H_12 = -2.0*(1.0-x_H_)*x_H_/(W()+W_H_) * R_*T/(p*D_H_H2) * (1.0 - 3.0*A_/5.0)

    # return((pow(1.0-x_H_, 2.0)/H_11 + pow(x_H_, 2.0)/H_22 - 2.0*(1.0-x_H_)*x_H_*H_12/(H_11*H_22))/ (1.0 - pow(H_12, 2.0)/(H_11*H_22)))


 # MIX
    
    x_H_ = x_H(T, p)
    D_H_H2 = diffusion(T, p)
    A_ = A_12(T)

    if x_H_ == 0:
        return viscosityH2_
    else :
        H_11 = pow(1-x_H_, 2.0)/viscosityH2_ + 2*(1-x_H_)*x_H_/(W()+W_H_) * R_*T/(p*D_H_H2) * (1 + 3*W_H_*A_/(5*W()))
        H_22 = pow(x_H_, 2.0)/viscosityH_ + 2*(1-x_H_)*x_H_/(W()+W_H_) * R_*T/(p*D_H_H2) * (1 + 3*W()*A_/(5*W_H_))
        H_12 = -2*(1-x_H_)*x_H_/(W()+W_H_) * R_*T/(p*D_H_H2) * (1 - 3*A_/5)

        return(((1-x_H_)**2/H_11 + x_H_**2/H_22 - 2*(1-x_H_)*x_H_*H_12/(H_11*H_22)) / (1 - pow(H_12, 2.0)/(H_11*H_22)))


def thermalConductivityH2(T: float, p: float) -> float:
    """
    Expressions from:
    - [8] 
    - lambda_0 expression from Eq. (2), coefficients Aij found in Table 2.
    - lambda_excess expression from Eq. (3), coefficients Bij found in Table 2.
    - lambda_critical expression from Eq. (8), coefficients Ci found in the text following Eq. (8).

    Note :  lambda_critical used here is obtained via empirical correlation,
    the deviation for the Roder and Diller data is -8% to +15%, while 
    those for the simplified Olchowy-Sengers enhancement vary from -9% to +9%.


    T: Temperature in K
    p: Pressure in Pa
    return thermal conductivity of H2 in W/m/K
    """

    a_1_list = [-3.40976e-1, 4.58820, -1.45080, 3.26394e-1, 3.16939e-3, 1.90592e-4, -1.139e-6]
    a_2_list = [1.38497e2, -2.21878e1, 4.57151, 1.0]

    Tc_ = 33.45 * K

    lambda_0_num = sum([a_1i * pow(T/Tc_, i) for i, a_1i in enumerate(a_1_list)])
    lambda_0_denom = sum([a_2i * pow(T/Tc_, i) for i, a_2i in enumerate(a_2_list)])

    lambda_0 = lambda_0_num/lambda_0_denom

    b_1_list = [3.63081e-2, -2.07629e-2, 3.14810e-2, -1.43097e-2, 1.74980e-3]
    b_2_list = [1.83370e-3, -8.86716e-3, 1.58260e-2, -1.06283e-2, 2.80673e-3]

    # Molecular hydrogen
    rhoH2 = rho(T, p)
    rho_c = 31.262 * kg/m3

    lambda_excess = sum(
        (b_1i + b_2i*T/Tc_)*pow(rhoH2/rho_c, i+1.0) for i, b_1i, b_2i in zip(
            range(0, len(b_1_list)), b_1_list, b_2_list
        )
    )

    c1 = 6.24e-4 # W/m/K
    c2 = -2.58e-7
    c3 = 0.837

    lambda_critical = c1/(c2 + T/Tc_ -1.0) * exp(-pow(c3*(rhoH2/rho_c-1.0), 2.0))   #empirical, a theoretical model exists but more complex

    return lambda_0 + lambda_excess + lambda_critical


def thermalConductivity(T: float, p: float) -> float: 
    """
    Expressions from:
    - [1] 

    T: Temperature in K
    p: Pressure in Pa
    return thermal conductivity of the mix H-H2 in W/m/K
    """

    x_H_ = x_H(T, p)
    W_H_H2_ = x_H_*W_H_ + (1.0-x_H_)*W()
    A12 = A_12(T)
    B12 = B_12(T)
    lambdaH2 = thermalConductivityH2(T, p)

    lambda_H = 15.0/4.0 * R_*viscosityH(T)/W_H_

    lambda_12 = 15.0/4.0 * R_/W_H_H2_ * 2.6693e-4*sqrt(M_12*T)/omega22_12(T)
    U_y = 4.0/15.0 * A12 * W_H_*W()/pow(M_12, 2.0) * pow(lambda_12, 2.0)/(lambda_H*lambdaH2) - (B12/5.0+1.0/12.0) - 1.0/A12*(3*B12/8-25/32) * pow(W_H_-W(), 2.0)/(W_H_*W())
    U_z = 4.0/15.0 * A12 * (lambda_12*(1.0/lambda_H+1.0/lambdaH2)*W_H_*W()/pow(M_12, 2.0)-1.0) - (B12/5.0+1.0/12.0)
    U_1 = 4.0/15.0 * A12 - (B12/5.0+1.0/12.0)*W_H_/W() + pow(W_H_-W(), 2.0)/(2.0*W_H_*W())
    U_2 = 4.0/15.0 * A12 - (B12/5.0+1.0/12.0)*W()/W_H_ + pow(W()-W_H_, 2.0)/(2.0*W_H_*W())

    X_lambda = pow(x_H_, 2.0)/lambda_H + 2.0*x_H_*(1.0-x_H_)/lambda_12 + pow(1.0-x_H_, 2.0)/lambdaH2
    Y_lambda = pow(x_H_, 2.0)/lambda_H*U_1 + 2.0*x_H_*(1.0-x_H_)/lambda_12*U_y + pow(1.0-x_H_, 2.0)/lambdaH2*U_2
    Z_lambda = pow(x_H_, 2.0)*U_1 + 2.0*x_H_*(1.0-x_H_)*U_z + pow(1.0-x_H_, 2.0)*U_2

    lambda_f = (1.0+Z_lambda)/(X_lambda+Y_lambda)

    delta_H = 2.0*W_H_*specificEnthalpyH(T) - W()*specificEnthalpyH2(T, p)

    lambda_r = p*diffusion(T, p)/T * pow(delta_H/(R_*T), 2.0) * x_H_*(1.0-x_H_)/pow(2.0-x_H_, 2.0)

    return(lambda_f + lambda_r)


def psi(T: float, p: float) -> float:
    """
    T: Temperature in K
    p: Pressure in Pa
    return compressibility factor psi in s2/m2
    """
    return(rho(T, p)/p)


def Z(T: float, p: float) -> float:
    """
    To be verified
    T: Temperature in K
    p: Pressure in Pa
    return compression factor of mix H-H2, no unit
    """
    x_H_ = x_H(T, p)
    W_H_H2_ = x_H_*W_H_ + (1.0-x_H_)*W()
    return(p*W_H_H2_/(rho(T, p)*R_*T))


def bulkModulus(T: float, p: float) -> float:
    """
    To be verified
    T: Temperature in K
    p: Pressure in Pa
    return bulk modulus as rho * dP/drho
    Use the ARK model to derive the bulk modulus of mix H-H2
    """
    # Is it correct to use the mixed molar mass? If yes, should we update it elsewhere? 
    x_H_ = x_H(T, p)
    W_H_H2_ = x_H_*W_H_ + (1.0-x_H_)*W()

    n_ = 0.4986 + 1.1735*omega() + 0.4754*pow(omega(), 2.0)
    alpha0_ = 0.42748*pow(R_*Tc(), 2.0) / Pc()
    alpha_ = alpha0_ / pow(T/Tc(), n_)
    b_ = 0.08664*R_*Tc() / Pc()
    c_ = R_*Tc()/(Pc() + alpha0_/(Vc() * (Vc()+b_))) + b_ - Vc()

    rho_ = rho(T, p)

    dP_drho = W_H_H2_*R_*T/pow(W_H_H2_-rho_*(b_-c_), 2.0) - alpha_*rho_*W_H_H2_*(2.0*W_H_H2_+b_*rho_)/pow(pow(W_H_H2_, 2.0)+b_*W_H_H2_*rho_, 2.0)
    return(rho_*dP_drho)


def CpMCv(T: float, p: float) -> float:
    """
    T: Temperature in K
    p: Pressure in Pa
    return cp-cv in J/kg/K
    https://en.wikipedia.org/wiki/Specific_heat_capacity#Relation_between_specific_heat_capacities
    """
    Rmass=R_/W()

    n_ = 0.4986 + 1.1735*omega() + 0.4754*pow(omega(), 2.0)
    alpha0_ = 0.42748*pow(Rmass*Tc(), 2.0) / Pc()
    alpha_ = alpha0_ / pow(T/Tc(), n_)
    b_ = 0.08664*Rmass*Tc() / Pc()
    c_ = Rmass*Tc()/(Pc() + alpha0_/(Vc() * (Vc()+b_))) + b_ - Vc()
    V = 1.0/densityH2(T,p)

    dPdT = Rmass/(V - b_ + c_) + n_*alpha_/(T*V*(V+b_))
    dPdV = -Rmass*T / pow((V-b_+c_), 2) + alpha_*(b_+ 2*V)/pow(V*(V+b_), 2)
        
    # dPdT = R_/(V - b_ + c_) 
    # dPdV = -R_*T / sqrt(V-b_+c_) + alpha_/(V*sqrt(V+b_))+ alpha_/((V+b_)*sqrt(V))

    dVdT = - dPdT / dPdV   #only if eq of state can be written as f(x,y,z)=0, true for cubic equations

    return(T*dPdT*dVdT)


def Cv(T: float, p: float) -> float:
    """
    T: Temperature in K
    p: Pressure in Pa
    return isochoric specific heat capacity in J/kg/K
    """
    return(specificHeatCapacity(T, p) - CpMCv(T, p))


# def specificInternalEnergyDeparture(T: float, p: float) -> float:
#     """
#     T: Temperature in K
#     p: Pressure in Pa
#     return the specific internal energy departure in J/kg
#     """
#     n_ = 0.4986 + 1.1735*omega() + 0.4754*pow(omega(), 2.0)
#     alpha0_ = 0.42748*pow(R_*Tc(), 2.0) / Pc()
#     alpha_ = alpha0_ / pow(T/Tc(), n_)
#     b_ = 0.08664*R_*Tc() / Pc()

#     x_H_ = x_H(T, p)
#     W_H_H2_ = x_H_*W_H_ + (1.0-x_H_)*W()
#     V = W_H_H2_/rho_(T, p)
#     return(alpha_/b_*(1.0+n_)*log((V+b_)/V) / W_H_H2_)

# def idealSpecificInternalEnergy(T):
#     u_ref = 0
#     return cvPerfectGas(T)*T+u_ref

# def specificInternalEnergy(T: float, p: float) -> float:
#     """
#     To be verified
#     T: Temperature in K
#     p: Pressure in Pa
#     return the internal energy in J/kg
#     """
#     # return(alpha_/b_*(1.0+n_)*log((V+b_)/V) + Cv(T, p)*T)
#     return(idealSpecificInternalEnergy(T)-specificInternalEnergyDeparture(T,p))


def specificInternalEnergy(T, p):
    """
    By definition H = U+pV

    T: Temperature in K
    p: Pressure in Pa
    return specificInternalEnergy of H-H2 mix in J/kg
    """
    # x_H_ = x_H(T, p)
    # W_H_H2_ = x_H_*W_H_ + (1.0-x_H_)*W()
    V = 1/rho(T, p)

    return(specificEnthalpy(T, p) - p*V)


def idealEntropyH2(T, p):
    """
    From thermodynamic course, also on https://fr.wikipedia.org/wiki/Entropie_d%27un_corps_pur
    Entropy ref from NIST
    T: Temperature in K
    p: Pressure in Pa
    return ideal gas entropy of H2 in J/kg/K
    """

    Tref = 273.15*K
    Pref = 101325*Pa
    s_ref = 130.68 # J/mol/K

    return (s_ref + cpPerfectGas(T)*W()*log(T/Tref) - R_*log(p/Pref)) / W()


def idealEntropyH(T, p):
    """
    from thermodynamic course, also on https://fr.wikipedia.org/wiki/Entropie_d%27un_corps_pur
    T: Temperature in K
    p: Pressure in Pa
    return ideal gas entropy of H in J/kg/K
    """

    Tref = 273.15*K
    Pref = 101325*Pa
    s_ref = 114.72 # J/mol/K

    return (s_ref + specificHeatCapacityH(T)*W_H_*log(T/Tref) - R_*log(p/Pref)) / W_H_


def entropyH2(T: float, p: float) -> float:
    """
    https://www.afs.enea.it/project/neptunius/docs/fluent/html/ug/node335.htm for entropy departure
    From NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=1#Thermo-Gas
    T: Temperature in K
    p: Pressure in Pa
    return the real gas entropy of H2 in J/kg/K
    """
    # Molecular hydrogen
    n_ = 0.4986 + 1.1735*omega() + 0.4754*pow(omega(), 2.0)
    alpha0_ = 0.42748*pow(R_*Tc(), 2.0) / Pc() # J*m3/(mol^2)
    alpha_ = alpha0_ / pow(T/Tc(), n_)
    b_ = 0.08664*R_*Tc() / Pc() # m3/mol
    # -> alpha/b = J/mol
    c_ = R_*Tc()/(Pc() + alpha0_/(Vc() * (Vc()+b_))) + b_ - Vc()

    V = W()/densityH2(T, p) # m3/mol

    entropyDeparture = (-R_*log((V-b_+c_)*p/(R_*T)) + n_*alpha_/(T*b_)*log((V+b_)/V))/W()

    return(idealEntropyH2(T, p) - entropyDeparture)


def entropy(T, p):
    """
    T: Temperature in K
    p: Pressure in Pa
    return real gas entropy of H-H2 mix in J/kg
    """
    # Mix
    w_H_ = w_H(T, p)

    return(w_H_*idealEntropyH(T, p) + (1.0-w_H_)*entropyH2(T, p))


def secondVirialCoeff(T):
    """
    Garberoglio, Giovanni, Piotr Jankowski,
    Krzysztof Szalewicz, et Allan H. Harvey.
    « Second Virial Coefficients of H 2 and 
    Its Isotopologues from a Six-Dimensional Potential ».
    The Journal of Chemical Physics 137, nᵒ 15 (21 octobre 2012): 154308. 
    https://doi.org/10.1063/1.4757565.

    T: Temperature in K
    """

    x_H_ = x_H(T, p)
    W_H_H2_ = x_H_*W_H_ + (1.0-x_H_)*W()

    if (T < 50):
        coeffList=[-5.66667e-7, 1.26295e-4, -1.16747e-2, 0.578078, -16.37, 258.895, -1933.3 ]
    
    # elif(T>50 and T<=400):
    #     coeffList=[-1.25632e-12, 1.6811e-9, -9.37815e-7, 2.72506e-4, -4.45883e-2, 4.03708, -155.946]
    # elif(T>400 and T<=700):
    #     coeffList=[-1.42222e-12, 3.75556-9, -3.89111e-6, 1.96061e-3, -4.72363e-1, 58.14]
    
    else:
        # coeffList=[3.06277e-15, -2.02294e-11, 5.2562e-8, -6.73655e-5, 4.09781e-2, 6.832 ] 
        return(15*1e-6/W_H_H2_)

    # Reverse the list
    B = sum([a * pow(T, i) for i, a in enumerate(coeffList[::-1])]) # cm3/mol

    return(B*1e-6/W_H_H2_) # m3/kg


def Prandtl(T, p):
    """
    T: Temperature in K
    p: Pressure in Pa
    return Prandtl number
    """
    return(specificHeatCapacity(T, p)*dynamicViscosity(T, p) / thermalConductivity(T, p))


if (__name__ == '__main__'):

    # T, p =  1933.1986 * K, 3416788.4 * Pa
    # print(f"T = {T:.3f}; p = {p:.3f} Pa")
    # print(f"rho                  = {rho(T, p):.7f} kg/m3")
    # print(f"viscosity            = {dynamicViscosity(T, p):.8e} Pa.s")
    # print(f"specificHeatCapacity = {specificHeatCapacity(T, p):.7f} J/kg")
    # print(f"thermalConductivity  = {thermalConductivity(T, p):.7f} W/m/K")
    # print(f"Prandtl              = {Prandtl(T, p):.7f}")

    # sys.exit(0)

    Trange = linspace(10, 6000, 10000)
    cpPerfectGasList = [cpPerfectGas(T) for T in Trange]
    plt.plot(Trange, cpPerfectGasList, label='Reference')
        
    guess = array([
        -9.87035830e-05, 1.70699738e+00, 1.41705707e+04, 4.70123961e+03, 
        2.32642127e-02, 9.70407310e+03, 1.37874475e+02, 2.99984374e-03, 
        1.10241106e+03
    ])
    popt, pcov = curve_fit(func, Trange, cpPerfectGasList, guess)
    # a, b, c, d, e, f, g, h, i = popt

    r_squared = computeRsquare(cpPerfectGasList, func(Trange, *popt))

    plt.plot(Trange, func(Trange, *popt), '--', label=f'Fit R$^2$ = {r_squared:.6f}')
    plt.xlabel("Temperature [K]")
    plt.ylabel(f"Specific heat capacity ideal gas H$_2$ [J/kg/K]")
    plt.grid(True)
    plt.legend()
    plt.savefig("images/specificHeatCapacityIdealH2fit.png")


    Tlist = [
        83.33334, 111.11112, 222.22224, 333.33336, 444.44448, 555.5556, 
        833.3334, 1111.1112, 1388.889, 1666.6668, 1944.4446, 2222.2224, 
        2500.0002, 2777.778, 3055.5558, 3333.3336
    ]

    plt.figure("Kp")
    # Seems correct
    TrangeKp = linspace(1111, 5000, 1000)
    plt.plot([1e4/T for T in TrangeKp], [log10(K_) for K_ in Kp(TrangeKp)])
    plt.xlabel("1e4/T [1/K]")
    plt.ylabel("log Kp")
    # plt.legend()
    plt.grid(True)
    plt.savefig("images/Kp.png")


    PG = PerfectGas()


    plt.figure("x_H")
    TrangeXH = linspace(80, 3500, 1000)
    plt.plot(TrangeXH, [x_H(T=T, p=0.1*atm) for T in TrangeXH], label="0.1 atm")
    plt.plot(TrangeXH, [x_H(T=T, p=1.0*atm) for T in TrangeXH], label="1 atm")
    plt.plot(TrangeXH, [x_H(T=T, p=10*atm) for T in TrangeXH], label="10 atm")
    plt.plot(TrangeXH, [x_H(T=T, p=100*atm) for T in TrangeXH], label="100 atm")
    plt.xlabel("Temperature [K]")
    plt.ylabel("Mole fraction of atomic hydrogen")
    plt.legend()
    plt.grid(True)
    plt.savefig("images/x_H.png")


    plt.figure("Specific heat capacity")
    TrangeHeatCap = linspace(80, 3500, 1000)
    plt.plot(TrangeHeatCap, [specificHeatCapacity(T=T, p=0.01*MPa) for T in TrangeHeatCap], label="0.01 MPa")
    plt.plot(TrangeHeatCap, [specificHeatCapacity(T=T, p=0.1*MPa) for T in TrangeHeatCap], label="0.1 MPa")
    plt.plot(TrangeHeatCap, [specificHeatCapacity(T=T, p=1*MPa) for T in TrangeHeatCap], label="1 MPa")
    plt.plot(TrangeHeatCap, [specificHeatCapacity(T=T, p=5*MPa) for T in TrangeHeatCap], label="5 MPa")
    plt.xlabel("Temperature [K]")
    plt.ylabel("Specific heat capacity [J/kg/K]")
    plt.grid(True)
    plt.legend()
    plt.savefig("images/specificHeatCapacity.png")


    plt.figure("Thermal conductivity")
    TrangeThermCond = linspace(80, 3500, 1000)
    # plt.plot(TrangeThermCond, [thermalConductivity(T=T, p=0.01*MPa) for T in TrangeThermCond], label="0.01 MPa")
    # plt.plot(TrangeThermCond, [thermalConductivity(T=T, p=0.1*MPa) for T in TrangeThermCond], label="0.1 MPa")
    # plt.plot(TrangeThermCond, [thermalConductivity(T=T, p=1*MPa) for T in TrangeThermCond], label="1 MPa")
    # plt.plot(TrangeThermCond, [thermalConductivity(T=T, p=10*MPa) for T in TrangeThermCond], label="10 MPa")
    # plt.plot(TrangeThermCond, [thermalConductivity(T=T, p=1*psia) for T in TrangeThermCond], label="1 psia")

    plt.plot(TrangeThermCond, [thermalConductivity(T=T, p=1*psia) for T in TrangeThermCond], label="1 psia")
    plt.plot(TrangeThermCond, [thermalConductivity(T=T, p=15*psia) for T in TrangeThermCond], label="15 psia")
    plt.plot(TrangeThermCond, [thermalConductivity(T=T, p=50*psia) for T in TrangeThermCond], label="50 psia")
    plt.plot(TrangeThermCond, [thermalConductivity(T=T, p=150*psia) for T in TrangeThermCond], label="150 psia")
    plt.plot(TrangeThermCond, [thermalConductivity(T=T, p=500*psia) for T in TrangeThermCond], label="500 psia")
    plt.plot(TrangeThermCond, [thermalConductivity(T=T, p=1500*psia) for T in TrangeThermCond], label="1500 psia")
    plt.plot(TrangeThermCond, [thermalConductivity(T=T, p=4500*psia) for T in TrangeThermCond], label="4500 psia")

    thermalConductivityList1 = [
        0.078090915, 0.087935136, 0.146328954, 0.197889969, 0.242677023, 
        0.287464078, 0.383685806, 0.48061539, 0.58168827, 0.71079849, 
        1.05088104, 2.1495294, 5.29853805, 11.2997403, 17.8487091, 
        17.687754
    ]
    thermalConductivityList15 = [
        0.078360904, 0.088161858, 0.146427604, 0.197912468, 0.242654524, 
        0.28739485, 0.383543888, 0.48061539, 0.581681347, 0.68951088, 
        0.854403323, 1.225422135, 2.190633525, 4.1969475, 7.733200275, 
        12.28364152
    ]
    thermalConductivityList50 = [
        0.079034146, 0.088727797, 0.146671633, 0.197971312, 0.242597411, 
        0.287223511, 0.383190826, 0.48061539, 0.581667502, 0.68587641, 
        0.81896724, 1.067755365, 1.627593548, 2.80113795, 4.9117266, 
        8.063763975
    ]
    thermalConductivityList150 = [
        0.080956954, 0.090346001, 0.147370836, 0.198135728, 0.242434725, 
        0.286733723, 0.382181828, 0.48061539, 0.581627696, 0.68414571, 
        0.803967263, 0.992527028, 1.367685675, 2.117222423, 3.461110973, 
        5.5503549
    ]
    thermalConductivityList500 = [
        0.087687646, 0.096007121, 0.149816315, 0.198713782, 0.241865325, 
        0.285016868, 0.378647738, 0.48061539, 0.581485778, 0.665684333, 
        0.79300674, 0.93734712, 1.175894693, 1.606551697, 2.353175677, 
        3.531782377
    ]
    thermalConductivityList1500 = [
        0.109581001, 0.112189166, 0.156619697, 0.200363139, 0.240238467, 
        0.280113795, 0.368552565, 0.48061539, 0.581082525, 0.68241501, 
        0.78816078, 0.91311732, 1.09155249, 1.38040632, 1.8535797, 
        2.5891272
    ]
    thermalConductivityList4500 = [
        0.160783761, 0.155295711, 0.174672628, 0.205312941, 0.235357893, 
        0.265402845, 0.338265315, 0.48061539, 0.579871035, 0.68206887, 
        0.78435324, 0.89338734, 1.02370905, 1.19868282, 1.4520573, 
        1.817235
    ]
    plt.scatter(Tlist, thermalConductivityList1) # 1 psia
    plt.scatter(Tlist, thermalConductivityList15) # 15 psia
    plt.scatter(Tlist, thermalConductivityList50) # 50 psia
    plt.scatter(Tlist, thermalConductivityList150) # 150 psia
    plt.scatter(Tlist, thermalConductivityList500) # 500 psia
    plt.scatter(Tlist, thermalConductivityList1500) # 1500 psia
    plt.scatter(Tlist, thermalConductivityList4500) # 4500 psia
    plt.xlabel("Temperature [K]")
    plt.ylabel("Thermal conductivity [W/m/K]")
    plt.yscale('log')
    plt.legend()
    plt.grid(True)
    plt.savefig("images/thermalConductivity.png")


    plt.figure("Diffusion")
    TrangeDiff = linspace(80, 3500, 1000)
    plt.plot(
        TrangeDiff, 
        [15*psia*diffusion(T=T, p=15*psia) for T in TrangeDiff], 
        label=f"15 psia = {15*psia/MPa:.3f} MPa"
    )
    p = 1*atm
    plt.scatter( #independent of pressure
        [1000, 1500, 2000, 2500, 3000, 3500], 
        [p*D_H_H2/100 for D_H_H2 in [0.172, 0.363, 0.620, 0.942, 1.33, 1.78]]
    )
    plt.xlabel("Temperature [K]")
    plt.ylabel("p * Diffusion [Pa.m2/s]")
    plt.legend()
    plt.grid(True)
    plt.savefig("images/diffusion.png")
    # print(f"D(83.3 K, 3.861 MPa) = {diffusion(83.3*K, 3.861*MPa)} m2/s")

    
    Trange = linspace(60, 3500, 400)

    plt.figure("Specific enthalpy")
    # plt.plot(Trange, [specificEnthalpy(T=T, p=0.01*MPa) for T in Trange], label="0.01 MPa")
    # plt.plot(Trange, [specificEnthalpy(T=T, p=0.1*MPa) for T in Trange], label="0.1 MPa")
    # plt.plot(Trange, [specificEnthalpy(T=T, p=1*MPa) for T in Trange], label="1 MPa")
    # plt.plot(Trange, [specificEnthalpy(T=T, p=10*MPa) for T in Trange], label="10 MPa")

    # Data from https://ntrs.nasa.gov/api/citations/19720005127/downloads/19720005127.pdf

    # plt.scatter( #0.01 MPa
    #     [T             for T   in [100, 500, 1000, 1500, 2000, 2500, 3000, 3400]], 
    #     [h*1000 for h in [2160, 7790, 15130, 22950, 36910, 123328, 260408, 284401]]
    # )


    # plt.scatter( #0.01 MPa
    #     [80, 120, 200, 300, 400, 500], 
    #     [h*1000 for h in [1354.3, 1803.6, 2825.9, 4226.9, 5668.5, 7118.4]]
    # )
    # plt.scatter( #0.1 MPa
    #     [80, 120, 200, 300, 400, 500], 
    #     [h*1000 for h in [1354.3, 1803.6, 2825.9, 4226.9, 5668.5, 7118.4]]
    # )
    # plt.scatter( #1 MPa
    #     [80, 120, 200, 300, 400, 500], 
    #     [h*1000 for h in [1333.7, 1795, 2825.6, 4230.7, 5673.9, 7124.6]]
    # )
    # plt.scatter( #10 MPa
    #     [80, 120, 200, 300, 400, 500], 
    #     [h*1000 for h in [1192.2, 1739.4, 2834.5, 4273.2, 5730.7, 7188.8]]
    # )
    # plt.scatter( #30 MPa
    #     [80, 120, 200, 300, 400, 500], 
    #     [h*1000 for h in [1232.2, 1773, 2915.7, 4391.5, 5869.6, 7338.9]]
    # )
    # plt.scatter( #100 MPa
    #     [80, 120, 200, 300, 400, 500], 
    #     [h*1000 for h in [1795.8, 2287.3, 3403.4, 4898.1, 6404.6, 7893.4]]
    # )

    plt.plot(Trange, [specificEnthalpy(T=T, p=1*psia) for T in Trange], label="1 psia")
    plt.plot(Trange, [specificEnthalpy(T=T, p=15*psia) for T in Trange], label="15 psia")
    plt.plot(Trange, [specificEnthalpy(T=T, p=50*psia) for T in Trange], label="50 psia")
    plt.plot(Trange, [specificEnthalpy(T=T, p=150*psia) for T in Trange], label="150 psia")
    plt.plot(Trange, [specificEnthalpy(T=T, p=500*psia) for T in Trange], label="500 psia")
    plt.plot(Trange, [specificEnthalpy(T=T, p=1500*psia) for T in Trange], label="1500 psia")
    plt.plot(Trange, [specificEnthalpy(T=T, p=4500*psia) for T in Trange], label="4500 psia")

    hlist1 = [ 
        1245433.44, 1615755.9, 3128749.12, 4712987.72, 6318858.12,
        7925147.2, 11981272.52, 16130856.52, 20429304.52, 24934720,
        30005400, 36983400, 50253230, 76664960, 126371556.7, 
        192871873.5 
    ]
    hlist15 = [ 
        1243386.56, 1614941.8, 3129051.5, 4713569.22, 6319532.66, 
        7925496.1, 11982063.36, 16131647.36, 20430072.1, 24934720, 
        29633240, 34971410, 42106415, 52928130, 71227935, 
        100192450 
    ]
    hlist50 = [ 
        1238292.62, 1612964.7, 3129795.82, 4715057.86, 6321230.64, 
        7927403.42, 11984040.46, 16133624.46, 20432002.68, 24934720, 
        29557645, 34645770, 40693370, 48863445, 61022610, 
        79560830 
    ]
    hlist150 = [
        1223731.86, 1607312.52, 3131889.22, 4719314.44, 6326091.98, 
        7932846.26, 11989669.38, 16139253.38, 20437492.04, 24934720, 
        29544084.42, 34490695.58, 40053720, 46942564.42, 56242680, 
        46202105.58 
    ]##################TO CHEEEECK
    hlist500 = [ 
        1176141.9, 1587518.26, 3139216.12, 4734200.84, 6343048.52, 
        7951896.2, 12009370.6, 16158954.6, 20456704.8, 24957980, 
        29540200, 34378280, 39604034.42, 45535334.42, 52753680, 
        61995645.58
    ]
    hlist1500 = [
        1089754.26, 1535439.12, 3160754.88, 4776720.12, 6391522.36,
        8006324.6, 12065659.8, 16215243.8, 20511598.4, 25004500,
        29586720, 34378280, 39448960, 44961580, 51265040,
        58754760 
    ]
    hlist4500 = [
        1136227.74, 1574143.76, 3257423.44, 4904277.96, 6536897.36,
        8169609.8, 12234527.4, 16384111.4, 20676279.2, 25144060,
        29726280, 34517840, 39448960, 44635940, 50195080, 
        56289200 
    ]
    plt.scatter(Tlist, hlist1) #, label="1 psia")
    plt.scatter(Tlist, hlist15) #, label="15 psia")
    plt.scatter(Tlist, hlist50) #, label="50 psia")
    plt.scatter(Tlist, hlist150) #, label="150 psia")
    plt.scatter(Tlist, hlist500) #, label="500 psia")
    plt.scatter(Tlist, hlist1500) #, label="1500 psia")
    plt.scatter(Tlist, hlist4500) #, label="4500 psia")

    plt.xlabel("Temperature [K]")
    plt.ylabel("Specific enthalpy [J/kg]")
    plt.yscale('log')
    plt.legend()
    plt.grid(True)
    plt.savefig("images/specificEnthalpy.png")


    plt.figure("Specific enthalpy low temperature")
    TrangeSpecificEnthalpy = linspace(80, 500, 10)
    plt.plot(TrangeSpecificEnthalpy, [specificEnthalpy(T=T, p=0.1*MPa) for T in TrangeSpecificEnthalpy], label="0.1 MPa")
    plt.plot(TrangeSpecificEnthalpy, [specificEnthalpy(T=T, p=1*MPa) for T in TrangeSpecificEnthalpy], label="1 MPa")
    plt.plot(TrangeSpecificEnthalpy, [specificEnthalpy(T=T, p=10*MPa) for T in TrangeSpecificEnthalpy], label="10 MPa")
    plt.plot(TrangeSpecificEnthalpy, [specificEnthalpy(T=T, p=30*MPa) for T in TrangeSpecificEnthalpy], label="30 MPa")
    plt.plot(TrangeSpecificEnthalpy, [specificEnthalpy(T=T, p=100*MPa) for T in TrangeSpecificEnthalpy], label="100 MPa")
    # plt.plot(TrangeSpecificEnthalpy, [specificEnthalpyH2NIST(T=T, p=1*MPa) for T in TrangeSpecificEnthalpy], ':', label="1 MPa NIST")

    TrangeRef = [80, 120, 200, 300, 400, 500]
    plt.scatter( # 0.1 MPa
        TrangeRef, 
        [h*1000 for h in [1354.3, 1803.6, 2825.9, 4226.9, 5668.5, 7118.4]]
    )
    plt.scatter( # 1 MPa
        TrangeRef, 
        [h*1000 for h in [1333.7, 1795, 2825.6, 4230.7, 5673.9, 7124.6]]
    )
    plt.scatter( # 10 MPa
        TrangeRef, 
        [h*1000 for h in [1192.2, 1739.4, 2834.5, 4273.2, 5730.7, 7188.8]]
    )
    plt.scatter( # 30 MPa
        TrangeRef, 
        [h*1000 for h in [1232.2, 1773, 2915.7, 4391.5, 5869.6, 7338.9]]
    )
    plt.scatter( # 100 MPa
        TrangeRef, 
        [h*1000 for h in [1795.8, 2287.3, 3403.4, 4898.1, 6404.6, 7893.4]]
    )
    plt.xlabel("Temperature [K]")
    plt.ylabel("Specific enthalpy [J/kg]")
    plt.legend()
    plt.grid(True)
    plt.savefig("images/specificEnthalpyLowTemperature_Th.png")


    plt.figure("Specific enthalpy low temperature (p-h)")
    pRange = linspace(0.1*MPa, 100*MPa, 10)
    plt.plot(pRange, [specificEnthalpy(T=70*K, p=P) for P in pRange], label="70 K")
    plt.plot(pRange, [specificEnthalpy(T=120*K, p=P) for P in pRange], label="120 K")
    plt.plot(pRange, [specificEnthalpy(T=160*K, p=P) for P in pRange], label="160 K")
    plt.plot(pRange, [specificEnthalpy(T=220*K, p=P) for P in pRange], label="220 K")
    plt.plot(pRange, [specificEnthalpy(T=280*K, p=P) for P in pRange], label="280 K")
    plt.plot(pRange, [specificEnthalpy(T=300*K, p=P) for P in pRange], label="300 K")
    plt.plot(pRange, [specificEnthalpy(T=350*K, p=P) for P in pRange], label="350 K")
    plt.plot(pRange, [specificEnthalpy(T=400*K, p=P) for P in pRange], label="400 K")

    ## Data from [10]
    pRange = [p*MPa for p in [1, 1.5, 2, 20, 50, 80]]
    plt.scatter( #70 K
        pRange,
        [h*1000 for h in [1221.7, 1207.3, 1193.2, 1043.6, 1260.7, 1514.1]]
    )
    plt.scatter( #120 K
        pRange,
        [h*1000 for h in [1795, 1790.3, 1785.9, 1735.9, 1897.4, 2127.2]]
    )
    plt.scatter( #160 K vsdv
        pRange,
        [h*1000 for h in [2294.2, 2292.4, 2290.8, 2296.6, 2455.8, 2673.8]]
    )
    plt.scatter( #220 K
        pRange,
        [h*1000 for h in [3099.9, 3100.4, 3101.0, 3157.9, 3334, 3551.9]]
    )
    plt.scatter( #280 K
        pRange,
        [h*1000 for h in [3944.5, 3946.2, 3948.0, 4031.5, 4234.8, 4457.2]]
    )
    plt.scatter( #300 K odf
        pRange,
        [h*1000 for h in [4230.7, 4232.8, 4234.9, 4329.2, 4528, 4747.6]]
    )
    plt.scatter( #350 K gref
        pRange,
        [h*1000 for h in [4950.5, 4953.2, 4955.8, 5064.3, 5275, 5500.5]]
    )
    plt.scatter( #400 K
        pRange,
        [h*1000 for h in [5673.9, 5677.0, 5680, 5798.4, 6018.3, 6249.0]]
    )
    plt.xlabel("Pressure [Pa]")
    plt.ylabel("Specific enthalpy [J/kg]")
    plt.legend()
    plt.grid(True)
    plt.savefig("images/specificEnthalpyLowTemperature_ph.png")


    plt.figure("Viscosity")
    # plt.plot(Trange, [dynamicViscosity(T=T, p=0.01*MPa) for T in Trange], label="0.01 MPa")
    # plt.plot(Trange, [dynamicViscosity(T=T, p=0.1*MPa) for T in Trange], label="0.1 MPa")
    # plt.plot(Trange, [dynamicViscosity(T=T, p=1*MPa) for T in Trange], label="1 MPa")
    # plt.plot(Trange, [dynamicViscosity(T=T, p=10*MPa) for T in Trange], label="10 MPa")

    plt.plot(Trange, [dynamicViscosity(T=T, p=1*psia) for T in Trange], label="1 psia")
    plt.plot(Trange, [dynamicViscosity(T=T, p=15*psia) for T in Trange], label="15 psia")
    plt.plot(Trange, [dynamicViscosity(T=T, p=50*psia) for T in Trange], label="50 psia")
    plt.plot(Trange, [dynamicViscosity(T=T, p=150*psia) for T in Trange], label="150 psia")
    plt.plot(Trange, [dynamicViscosity(T=T, p=500*psia) for T in Trange], label="500 psia")
    plt.plot(Trange, [dynamicViscosity(T=T, p=1500*psia) for T in Trange], label="1500 psia")
    plt.plot(Trange, [dynamicViscosity(T=T, p=4500*psia) for T in Trange], label="4500 psia")
    
    # plt.scatter( #0.01 MPa
    #     [T             for T   in [80, 120, 200, 300, 400, 500]], 
    #     [mu*pow(10,-7) for mu in [35.77, 54.09, 69.01, 81.39, 92.95, 104.23]]
    # )
    # plt.scatter( #0.1 MPa
    #     [T             for T   in [80, 120, 200, 300, 400, 500]], 
    #     [mu*pow(10,-7) for mu in [35.85, 54.08, 69.01, 81.41, 92.97, 104.26]]
    # )
    # plt.scatter( #1 MPa
    #     [T             for T   in [80, 120, 200, 300, 400, 500]], 
    #     [mu*pow(10,-7) for mu in [36.64, 54.01, 69.09, 81.56, 93.20, 104.56]]
    # )
    # plt.scatter( #10 MPa
    #     [T             for T   in [50, 80, 120, 200, 300, 400, 500]], 
    #     [mu*pow(10,-7) for mu in [67.02, 50.28, 57.97, 71.14, 83.49, 95.56, 107.5]]
    # )
    muList1 = [
        3.68525E-06, 4.37989E-06, 7.24294E-06, 9.47167E-06, 1.15039E-05, 
        1.35344E-05, 1.78057E-05, 2.16495E-05, 2.52176E-05, 2.85615E-05, 
        3.17676E-05, 3.49392E-05, 3.83176E-05, 4.20752E-05, 4.51089E-05, 
        4.5626E-05
    ]
    muList15 = [
        3.69387E-06, 4.39024E-06, 7.24639E-06, 9.47339E-06, 1.15022E-05, 
        1.3531E-05, 1.77971E-05, 2.16392E-05, 2.52038E-05, 2.85529E-05, 
        3.17417E-05, 3.48237E-05, 3.79091E-05, 4.11927E-05, 4.48159E-05, 
        4.85081E-05
    ]
    muList50 = [
        3.71627E-06, 4.41609E-06, 7.25673E-06, 9.47512E-06, 1.1497E-05, 
        1.35206E-05, 1.77799E-05, 2.16133E-05, 2.51728E-05, 2.85081E-05, 
        3.169E-05, 3.47496E-05, 3.77902E-05, 4.09307E-05, 4.43522E-05, 
        4.80944E-05
    ]
    muList150 = [
        3.77833E-06, 4.48849E-06, 7.28603E-06, 9.48029E-06, 1.14867E-05, 
        1.34913E-05, 1.77264E-05, 2.15409E-05, 2.50814E-05, 2.84012E-05, 
        3.15676E-05, 3.46134E-05, 3.76333E-05, 4.07428E-05, 4.41144E-05, 
        4.77893E-05
    ]
    muList500 = [
        4.02481E-06, 4.74187E-06, 7.38773E-06, 9.49925E-06, 1.14453E-05, 
        1.33931E-05, 1.75385E-05, 2.12876E-05, 2.47642E-05, 2.80272E-05, 
        3.11419E-05, 3.41342E-05, 3.70817E-05, 4.00878E-05, 4.32818E-05, 
        4.6724E-05
    ]
    muList1500 = [
        5.1228E-06, 5.49167E-06, 7.67386E-06, 9.55096E-06, 1.13384E-05, 
        1.31259E-05, 1.703E-05, 2.05636E-05, 2.38559E-05, 2.69585E-05, 
        2.99232E-05, 3.27673E-05, 3.5508E-05, 3.82142E-05, 4.09031E-05, 
        4.36783E-05
    ]
    muList4500 = [
        8.70808E-06, 7.83244E-06, 8.51847E-06, 9.73884E-06, 1.11006E-05, 
        1.2464E-05, 1.57321E-05, 1.87365E-05, 2.15461E-05, 2.42351E-05, 
        2.67861E-05, 2.9251E-05, 3.16469E-05, 3.39739E-05, 3.62492E-05, 
        3.85417E-05
    ]
    plt.scatter(Tlist, muList1) # 1 psia
    plt.scatter(Tlist, muList15) # 15 psia
    plt.scatter(Tlist, muList50) # 50 psia
    plt.scatter(Tlist, muList150) # 150 psia
    plt.scatter(Tlist, muList500) # 500 psia
    plt.scatter(Tlist, muList1500) # 1500 psia
    plt.scatter(Tlist, muList4500) # 4500 psia
    plt.xlabel("Temperature [K]")
    plt.ylabel("Dynamic viscosity [Pa.s]")
    plt.yscale('log')
    plt.legend()
    plt.grid(True)
    plt.savefig("images/dynamicViscosity.png")


    plt.figure("Density")
    plt.plot(Trange, [rho(T=T, p=1*psia) for T in Trange], label="1 psia")
    plt.plot(Trange, [rho(T=T, p=15*psia) for T in Trange], label="15 psia")
    plt.plot(Trange, [rho(T=T, p=50*psia) for T in Trange], label="50 psia")
    plt.plot(Trange, [rho(T=T, p=150*psia) for T in Trange], label="150 psia")
    plt.plot(Trange, [rho(T=T, p=500*psia) for T in Trange], label="500 psia")
    plt.plot(Trange, [rho(T=T, p=1500*psia) for T in Trange], label="1500 psia")
    plt.plot(Trange, [rho(T=T, p=4500*psia) for T in Trange], label="4500 psia")
    densityList1 = [
        0.020101565, 0.015004491, 0.007325242, 0.005005769, 0.003756329, 
        0.003006665, 0.002005511, 0.001508939, 0.001204588, 0.001004357, 
        0.000858589, 0.000744858, 0.000639137, 0.000527007, 0.000410073, 
        0.00031236
    ]
    densityList15 = [
        0.301552315, 0.225059363, 0.109877025, 0.075086531, 0.056344933, 
        0.045099974, 0.030082668, 0.022637288, 0.018030379, 0.015031723, 
        0.012885249, 0.011251366, 0.00992664, 0.008754088, 0.007627991, 
        0.006519513
    ]
    densityList50 = [
        1.009060462, 0.750197343, 0.366255681, 0.250288438, 0.187816444,
        0.150333247, 0.10027556, 0.075459761, 0.059846568, 0.049935947,
        0.04283176, 0.037473585, 0.033223888, 0.029674197, 0.026539385,
        0.02365446
    ]
    densityList150 = [
        2.667643847, 2.25059363, 1.098767042, 0.747941944, 0.561847485,
        0.449998588, 0.299314536, 0.225530306, 0.179317049, 0.14965887,
        0.128397568, 0.112390321, 0.099785395, 0.089443877, 0.08056965,
        0.072608475
    ]
    densityList500 = [
        10.29643702, 7.501978233, 3.662553601, 2.446899857, 1.846127515,
        1.482588565, 0.990653649, 0.747693657, 0.597462928, 0.498691502,
        0.427878696, 0.374599696, 0.332754673, 0.298638557, 0.269673978,
        0.243949933
    ]
    densityList1500 = [
        29.33204765, 21.23380307, 10.4762955, 7.084964858, 5.387008098,
        4.34901189, 2.936984641, 2.22292174, 1.78117266, 1.488595488,
        1.27763237, 1.119530169, 0.995307012, 0.894310622, 0.809252599,
        0.734686668
    ]
    densityList4500 = [
        58.11917933, 47.59262431, 27.13438061, 19.77198564, 15.27520346,
        12.46156096, 8.541042872, 6.514307221, 5.242841958, 4.394664501,
        3.785963021, 3.319585558, 2.958129008, 2.662428237, 2.41462266,
        2.203259081
    ]
    plt.scatter(Tlist, densityList1) # 1 psia
    plt.scatter(Tlist, densityList15) # 15 psia
    plt.scatter(Tlist, densityList50) # 50 psia
    plt.scatter(Tlist, densityList150) # 150 psia
    plt.scatter(Tlist, densityList500) # 500 psia
    plt.scatter(Tlist, densityList1500) # 1500 psia
    plt.scatter(Tlist, densityList4500) # 4500 psia
    plt.xlabel("Temperature [K]")
    plt.ylabel("Density [kg/m3]")
    plt.yscale('log')
    plt.legend()
    plt.grid(True)
    plt.savefig("images/density.png")


    plt.figure("Specific heat capacity H2")
    # https://www.engineeringtoolbox.com/hydrogen-d_976.html data at high temp diverges with nasa data
    Trange = linspace(60, 3000, 2000)
    plt.plot(Trange, [specificHeatCapacityH2(T=T, p=0.1*MPa) for T in Trange], label="0.1 MPa")
    plt.plot(Trange, [specificHeatCapacityH2(T=T, p=1.5*MPa) for T in Trange], label="1.5 MPa")
    plt.plot(Trange, [specificHeatCapacityH2(T=T, p=3.5*MPa) for T in Trange], label="3.5 MPa")
    plt.plot(Trange, [specificHeatCapacityH2(T=T, p=4.0*MPa) for T in Trange], label="4.0 MPa")

    TrangeRef = [25, 30, 50, 90, 260, 280, 500, 800, 1200, 1600, 2000, 2500, 3000]
    plt.scatter( # 0.1 MPa 
        TrangeRef, 
        [CpNasa*1000 for CpNasa in [11.17, 10.83, 10.48, 10.96, 14.12, 14.23, 14.51, 14.69, 15.36, 16.30, 17.01, 17.80, 18.39]]
    )
    plt.scatter( # 1.5 MPa
        TrangeRef, 
        [CpNasa*1000 for CpNasa in [12.05, 18.68, 13.31, 11.51, 14.17, 14.28, 14.52, 14.70, 15.36, 16.25, 17.33, 17.80, 18.39]]
    )
    plt.scatter( # 3.5 MPa
        TrangeRef, 
        [CpNasa*1000 for CpNasa in [10.74, 13.73, 19.27, 12.26, 14.24, 14.33, 14.54, 14.70, 15.37, 16.25, 17.22, 17.80, 18.39 ]]
    )
    plt.scatter( # 4.0 MPa
        TrangeRef, 
        [CpNasa*1000 for CpNasa in [10.52, 13.22, 20.27, 12.44, 14.26, 14.35, 14.54, 14.70, 15.37, 16.25, 17.21, 17.80, 18.39]]
    )
    
    plt.xlabel("Temperature [K]")
    plt.ylabel(f"Specific heat capacity $c_p$ of $H_2$ [J/kg/K]")
    plt.grid(True)
    plt.legend()
    plt.savefig("images/specificHeatCapacityH2.png")


    plt.figure("CvH2")
    # https://www.engineeringtoolbox.com/hydrogen-d_976.html data at high temp diverges with nasa data
    plt.plot(Trange, [Cv(T=T, p=0.1*MPa) for T in Trange], label="0.1 MPa")
    plt.plot(Trange, [Cv(T=T, p=1.5*MPa) for T in Trange], label="1.5 MPa")
    plt.plot(Trange, [Cv(T=T, p=3.5*MPa) for T in Trange], label="3.5 MPa")
    plt.plot(Trange, [Cv(T=T, p=4.0*MPa) for T in Trange], label="4.0 MPa")
    
    TrangeRef = [25, 30, 50, 90, 260, 280, 500, 800, 1200, 1600]
    plt.scatter( #0.1 MPa 
        TrangeRef, 
        [CvNasa*1000 for CvNasa in [6.29, 6.25, 6.22, 6.8, 9.99, 10.10, 10.4, 10.57, 11.24, 12.17]]
    )
    plt.scatter( #1.5 MPa
        TrangeRef, 
        [CvNasa*1000 for CvNasa in [6.12, 6.38, 6.43, 6.84, 10, 10.11, 10.39, 10.57, 11.24, 12.13]]
    )
    plt.scatter( #3.5 MPa
        TrangeRef, 
        [CvNasa*1000 for CvNasa in [6.10, 6.33, 6.70, 6.9, 10.02, 10.13, 10.4, 10.57, 11.24, 12.12 ]]
    )
    plt.scatter( #4.0 MPa
        TrangeRef, 
        [CvNasa*1000 for CvNasa in [6.09, 6.33, 6.73, 6.92, 10.02, 10.13, 10.4, 10.57, 11.24, 12.12]]
    )
    plt.xlabel("Temperature [K]")
    plt.ylabel("Specific heat capacity $c_V$ [J/kg/K]")
    plt.grid(True)
    plt.legend()
    plt.savefig("images/CvH2.png")


    plt.figure("Prandtl_P")
    Prange = linspace(0.1*MPa, 4*MPa, 10)
    plt.plot(Prange, [Prandtl(T=70*K, p=P) for P in Prange], label="70 K")
    plt.plot(Prange, [Prandtl(T=120*K, p=P) for P in Prange], label="120 K")
    plt.plot(Prange, [Prandtl(T=160*K, p=P) for P in Prange], label="160 K")
    plt.plot(Prange, [Prandtl(T=220*K, p=P) for P in Prange], label="220 K")
    plt.plot(Prange, [Prandtl(T=280*K, p=P) for P in Prange], label="280 K")
    plt.plot(Prange, [Prandtl(T=300*K, p=P) for P in Prange], label="300 K")
    plt.plot(Prange, [Prandtl(T=350*K, p=P) for P in Prange], label="350 K")
    plt.plot(Prange, [Prandtl(T=400*K, p=P) for P in Prange], label="400 K")
    plt.xlabel("Pressure [Pa]")
    plt.ylabel("Prandtl number at constant temperature [-]")
    plt.legend()
    plt.grid(True)
    plt.savefig("images/Prandtl_P.png")



    plt.figure("Prandtl_T")
    TrangePrandtl = linspace(80, 3500, 2000)
    plt.plot(TrangePrandtl, [Prandtl(T, p=0.01*MPa) for T in TrangePrandtl], label="0.01 MPa")
    plt.plot(TrangePrandtl, [Prandtl(T, p=0.1*MPa) for T in TrangePrandtl], label="0.1 MPa")
    plt.plot(TrangePrandtl, [Prandtl(T, p=0.2*MPa) for T in TrangePrandtl], label="0.2 MPa")
    plt.plot(TrangePrandtl, [Prandtl(T, p=0.3*MPa) for T in TrangePrandtl], label="0.3 MPa")
    plt.plot(TrangePrandtl, [Prandtl(T, p=0.4*MPa) for T in TrangePrandtl], label="0.4 MPa")
    plt.plot(TrangePrandtl, [Prandtl(T, p=0.5*MPa) for T in TrangePrandtl], label="0.5 MPa")
    plt.plot(TrangePrandtl, [Prandtl(T, p=1*MPa) for T in TrangePrandtl], label="1 MPa")
    plt.plot(TrangePrandtl, [Prandtl(T, p=10*MPa) for T in TrangePrandtl], label="10 MPa")
    plt.plot(TrangePrandtl, [Prandtl(T, p=16*MPa) for T in TrangePrandtl], label="16 MPa")
    plt.xlabel("Temperature [K]")
    plt.ylabel("Prandtl number at constant pressure [-]")
    plt.grid(True)
    plt.legend()
    plt.savefig("images/Prandtl_T.png")


    plt.figure("Entropy")
    Trange = linspace(80, 3500, 2000)
    plt.plot(Trange, [entropy(T=T, p=1*psia) for T in Trange], label="1 psia")
    plt.plot(Trange, [entropy(T=T, p=15*psia) for T in Trange], label="15 psia")
    plt.plot(Trange, [entropy(T=T, p=50*psia) for T in Trange], label="50 psia")
    plt.plot(Trange, [entropy(T=T, p=150*psia) for T in Trange], label="150 psia")
    plt.plot(Trange, [entropy(T=T, p=500*psia) for T in Trange], label="500 psia")
    plt.plot(Trange, [entropy(T=T, p=1500*psia) for T in Trange], label="1500 psia")
    plt.plot(Trange, [entropy(T=T, p=4500*psia) for T in Trange], label="4500 psia")
    
    sList1 = [
        63953.37, 67960.1376, 77392.998, 82814.904, 87336.648, 90405.5724, 
        96451.3116, 100776.276, 104209.452, 107182.08, 109987.236, 113336.676, 
        118821.384, 128911.572, 145784.376, 166760.244
    ]
    sList15 = [
        52862.5368, 57020.0292, 66448.7028, 71887.356, 76409.1, 79478.0244, 
        85507.0164, 89848.728, 93051.63, 96254.532, 98850.348, 101446.164, 
        105017.5044, 108584.658, 116259.0624, 123929.28
    ]
    sList50 = [
        47876.058, 52062.858, 61529.2128, 66967.866, 71489.61, 74575.2816, 
        80587.5264, 84929.238, 88123.7664, 91314.108, 93847.122, 96380.136, 
        99227.16, 102074.184, 107119.278, 112164.372
    ]
    sList150 = [
        43174.2816, 47407.1364, 56940.48, 62383.32, 66905.064, 69990.7356, 
        76019.7276, 80344.692, 83547.594, 86750.496, 89233.2684, 91711.854, 
        94194.6264, 96673.212, 100307.3544, 103937.31
    ]
    sList500 = [
        37748.1888, 41997.7908, 51983.3088, 57451.2696, 61985.574, 65071.2456, 
        71104.4244, 75446.136, 78640.6644, 81831.006, 84280.284, 86729.562, 
        88977.8736, 91230.372, 94006.2204, 96777.882
    ]
    sList1500 = [
        32213.2392, 36726.6096, 47302.4664, 52837.416, 57401.028, 60499.26, 
        66519.8784, 70861.59, 74056.1184, 77246.46, 79687.3644, 82124.082, 
        84259.35, 86394.618, 88718.292, 91041.966
    ]
    sList4500 = [
        27230.9472, 31514.0436, 42202.944, 47143.368, 51790.716, 54888.948, 
        60905.3796, 65230.344, 68433.246, 71636.148, 74064.492, 76492.836, 
        78556.9284, 80616.834, 82576.2564, 84531.492
    ]
    # print(idealEntropyH(p=3.10264e+07*Pa, T=3333.33*K))
    # print(idealEntropyH2(p=3.10264e+07*Pa, T=3333.33*K))
    # print(entropyH2(p=3.10264e+07*Pa, T=3333.33*K))
    # print(entropy(p=3.10264e+07*Pa, T=3333.33*K), 84531.5)
    plt.scatter(Tlist, sList1) #1 psia
    plt.scatter(Tlist, sList15) # 15 psia
    plt.scatter(Tlist, sList50) # 50 psia
    plt.scatter(Tlist, sList150) # 150 psia
    plt.scatter(Tlist, sList500) # 500 psia
    plt.scatter(Tlist, sList1500) # 1500 psia
    plt.scatter(Tlist, sList4500) # 4500 psia
    plt.xlabel("Temperature [K]")
    plt.ylabel("Specific entropy [J/kg/K]")
    plt.yscale('log')
    plt.legend()
    plt.grid(True)
    plt.savefig("images/specificEntropy.png")


    plt.figure("Specific internal energy")
    plt.plot(Trange, [specificInternalEnergy(T=T, p=1*psia) for T in Trange], label="1 psia")
    plt.plot(Trange, [specificInternalEnergy(T=T, p=15*psia) for T in Trange], label="15 psia")
    plt.plot(Trange, [specificInternalEnergy(T=T, p=50*psia) for T in Trange], label="50 psia")
    plt.plot(Trange, [specificInternalEnergy(T=T, p=150*psia) for T in Trange], label="150 psia")
    plt.plot(Trange, [specificInternalEnergy(T=T, p=500*psia) for T in Trange], label="500 psia")
    plt.plot(Trange, [specificInternalEnergy(T=T, p=1500*psia) for T in Trange], label="1500 psia")
    plt.plot(Trange, [specificInternalEnergy(T=T, p=4500*psia) for T in Trange], label="4500 psia")
    plt.xlabel("Temperature [K]")
    plt.ylabel("Specific internal energy [J/kg]")
    plt.yscale('log')
    plt.legend()
    plt.grid(True)
    plt.savefig("images/specificInternalEnergy.png")


    plt.figure("Virial Coefficient")
    plt.plot(Trange, [secondVirialCoeff(T=T) for T in Trange])
    plt.xlabel("Temperature [K]")
    plt.ylabel("Virial coefficient [m3/kg]")
    # plt.legend()
    plt.grid(True)
    plt.savefig("images/virialCoefficient.png")

    # plt.show()
