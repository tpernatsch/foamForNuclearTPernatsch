/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2312                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2022 OpenCFD Ltd.         |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/


inline Foam::scalar Foam::H2::Kp(scalar T) const
{
    return pow(10.0, -2.37943e4/T + 6.33153);
}


inline Foam::scalar Foam::H2::x_H(scalar p, scalar T) const
{
    return 2.0 / (1.0 + sqrt(1.0 + 4.0*(p/101325.0)/Kp(T)));
}


inline Foam::scalar Foam::H2::w_H(scalar p, scalar T) const
{
    const scalar x_H_(x_H(p, T));
    return x_H_*W_H_ / (x_H_*W_H_ + (1.0-x_H_)*W());
}


//---------------------------------------------------------------------------//
// Definition of the average effective collision integrals

inline Foam::scalar Foam::H2::omega11_11(scalar T) const
{
    if (T <= 1000)
    {
        return(44.13*pow(T, -0.331));
    }
    return(
        -1.173958e-19*pow(T, 5.0)
        +3.983904e-15*pow(T, 4.0)
        -5.443976e-11*pow(T, 3.0)
        +3.918410e-07*sqr(T)
        -1.704143e-03*T
        +6.552751
    );
}


inline Foam::scalar Foam::H2::omega22_11(scalar T) const
{
    if (T <= 1000)
    {
        return(37.879*pow(T, -0.268));
    }
    return(
         5.570961e-20*pow(T, 5.0)
        -1.444695e-15*pow(T, 4.0)
        +1.024643e-11*pow(T, 3.0)
        +3.257067e-08*sqr(T)
        -8.379299e-04*T
        +6.786617
    );
}


// inline Foam::scalar Foam::H2::omega11_22(scalar T) const
// { 
//     // (not used)
//     // Ω(1,1){2,2} defined in [1]
//     // Coeffs fitted with data from [6], table II
//     return(
//         5e-16*pow(T, 4.0) - 2e-11*pow(T, 3.0) 
//         + 2e-7*pow(T, 2.0) - 0.0015*T + 6.401
//     );
// }


inline Foam::scalar Foam::H2::omega22_22(scalar T) const
{
    return(
        5e-16*pow(T, 4.0)
        -2e-11*pow(T, 3.0)
        +2e-7*sqr(T)
        -0.0015*T
        +7.1214
    );
}


inline Foam::scalar Foam::H2::omega11_12(scalar T) const
{
    return(-1.024*log(T) + 11.054);
}


inline Foam::scalar Foam::H2::omega22_12(scalar T) const
{
    return(-1.227*log(T) + 13.42);
}


//---------------------------------------------------------------------------//

inline Foam::scalar Foam::H2::A_12(scalar T) const
{
    return(0.0345*log(T) + 0.9917);
}


inline Foam::scalar Foam::H2::B_12(scalar T) const
{
    return(0.0462*log(T) + 0.8701);
}


inline Foam::scalar Foam::H2::C_12(scalar T) const
{
    return(-0.031*log(T) + 1.1108);
}


//---------------------------------------------------------------------------//

inline Foam::scalar Foam::H2::densityH(scalar p, scalar T) const
{
    return p * W_H_ / (R_*T);
}


inline Foam::scalar Foam::H2::densityH2(scalar p, scalar T) const
{
    // Molecular hydrogen
    // const scalar n_(0.4986 + 1.1735*omega() + 0.4754*pow(omega(), 2.0));
    // const scalar alpha0_(0.42748*pow(R_*Tc(), 2.0) / Pc());
    const scalar alpha_(alpha0_ / pow(T/Tc(), n_));
    // const scalar b_(0.08664*R_*Tc() / Pc());
    // const scalar c_(R_*Tc() / (Pc() + alpha0_ / (Vc() * (Vc()+b_))) + b_ - Vc());

    // A*rho**3 + B*rho**2 + C*rho + D = 0
    const scalar A_(alpha_ * (c_-b_));
    const scalar B_(W()*(alpha_ + p*b_*(c_-b_) - R_*T*b_ ));
    const scalar C_(sqr(W()) * (p*c_ - R_*T));
    const scalar D_(pow(W(), 3.0) * p);

    // Reduced form
    const scalar p_((3.0*A_*C_ - sqr(B_)) / (3.0*sqr(A_)));
    const scalar q_
    (
        (2.0*pow(B_, 3.0) - 9.0*A_*B_*C_ + 27.0*sqr(A_)*D_) / (27.0*pow(A_, 3.0))
    );

    // https://en.wikipedia.org/wiki/Cubic_equation//Trigonometric_and_hyperbolic_solutions
    const scalar discriminant(pow(p_, 3.0)/27.0 + sqr(q_)/4.0);

    if (discriminant < 0)
    {
        return(
            2.0*sqrt(-p_/3.0) * cos(acos(3.0*q_/(2.0*p_)*sqrt(-3.0/p_))/3.0) 
            - B_/(3.0*A_)
        );
    }
    else
    {
        const scalar u1(-q_/2.0 + sqrt(discriminant));
        const scalar u2(-q_/2.0 - sqrt(discriminant));
        return cbrt(u1) + cbrt(u2) - B_/(3.0*A_);
    }
}


inline Foam::scalar Foam::H2::rho(scalar p, scalar T) const
{
    const scalar w_H_(w_H(p, T));

    return(w_H_*densityH(p, T) + (1.0-w_H_)*densityH2(p, T));
}


inline Foam::scalar Foam::H2::psi(scalar p, scalar T) const
{
    return rho(p, T)/p;
}


inline Foam::scalar Foam::H2::CpMCv(scalar p, scalar T) const
{
    const scalar Rmass(R_/W());

    // const scalar n_(0.4986 + 1.1735*omega() + 0.4754*pow(omega(), 2.0));
    const scalar alpha0_spec(alpha0_ / sqr(W())); // 0.42748*pow(Rmass*Tc(), 2.0) / Pc());
    const scalar alpha_(alpha0_spec / pow(T/Tc(), n_));
    const scalar b_spec(b_ / W()); // (0.08664*Rmass*Tc() / Pc());
    const scalar Vc_(Vc()*1e-3);
    const scalar c_(Rmass*Tc()/(Pc() + alpha0_spec/(Vc_ * (Vc_+b_spec))) + b_spec - Vc_);
    const scalar V(1.0/densityH2(p, T));

    const scalar dPdT(Rmass/(V - b_spec + c_) + n_*alpha_/(T*V*(V+b_spec)));
    const scalar dPdV
    (
        -Rmass*T / sqr(V-b_spec+c_) + alpha_*(b_spec+ 2.0*V)/sqr(V*(V+b_spec))
    );

    // Only if eq of state can be written as f(x,y,z)=0, true for cubic 
    // equations
    const scalar dVdT(-dPdT / dPdV);

    return(T*dPdT*dVdT);
}


inline Foam::scalar Foam::H2::Hs(const scalar p, const scalar T) const
{
    return h(p, T);
}


inline Foam::scalar Foam::H2::Hc() const
{
    return 0;
}


inline Foam::scalar Foam::H2::Ha(const scalar p, const scalar T) const
{
    return h(p, T);
}


inline Foam::scalar Foam::H2::idealEntropyH2(const scalar p, const scalar T) const
{
    // const scalar Tref(273.15); // K
    // const scalar Pref(101325.0); // Pa
    // const scalar s_ref(130.68e3); // J/kmol/K

    return
    (
        // (s_ref + cpPerfectGasH2(T)*W()*log(T/Tref) - R_*log(p/Pref)) / W()
        (130.68e3 + cpPerfectGasH2(T)*W()*log(T/273.15) - R_*log(p/101325.0)) / W()
    );
}


inline Foam::scalar Foam::H2::idealEntropyH(const scalar p, const scalar T) const
{   
    // const scalar Tref(273.15); // K
    // const scalar Pref(101325.0); // Pa
    // const scalar s_ref(114.72e3); // J/kmol/K

    return
    (
        // (s_ref + heatCapacityH(T)*W_H_*log(T/Tref) - R_*log(p/Pref)) / W_H_
        (114.72e3 + heatCapacityH(T)*W_H_*log(T/273.15) - R_*log(p/101325.0)) / W_H_
    );
}


inline Foam::scalar Foam::H2::entropyH2(const scalar p, const scalar T) const
{
    // Molecular hydrogen
    // const scalar n_(0.4986 + 1.1735*omega() + 0.4754*pow(omega(), 2.0));
    // const scalar alpha0_(0.42748*pow(R_*Tc(), 2.0) / Pc()); // J*m3/(kmol^2)
    const scalar alpha_(alpha0_ / pow(T/Tc(), n_));
    // const scalar b_(0.08664*R_*Tc() / Pc()); // m3/kmol
    // -> alpha/b = J/kmol
    // const scalar c_(R_*Tc() / (Pc() + alpha0_ / (Vc() * (Vc()+b_))) + b_ - Vc());

    const scalar V(W()/densityH2(p, T)); // m3/kmol

    const scalar entropyDeparture
    (
        (-R_*log((V-b_+c_)*p/(R_*T)) + n_*alpha_/(T*b_)*log((V+b_)/V))/W()
    );

    return(idealEntropyH2(p, T) - entropyDeparture);
}


inline Foam::scalar Foam::H2::S(const scalar p, const scalar T) const
{
    const scalar w_H_(w_H(p, T));

    return(w_H_*idealEntropyH(p, T) + (1.0-w_H_)*entropyH2(p, T));
}


inline Foam::scalar Foam::H2::specificEnthalpyH(scalar T) const
{
    // const scalar A_(20.78603);
    // const scalar B_(4.850638e-10);
    // const scalar C_(-1.582916e-10);
    // const scalar D_(1.525102e-11);
    // const scalar E_(3.196347e-11);
    // const scalar F_(211.8020);
    // const scalar G_(139.8711);
    // const scalar H_(217.9994);

    const scalar T_(T/1000.0);

    /*const scalar enthalpyH(
        218.0 
        + A_*T_ 
        + B_*pow(T_, 2.0)/2.0 
        + C_*pow(T_, 3.0)/3.0 
        + D_*pow(T_, 4.0)/4.0 
        - E_/T_ + F_ - H_
    ); // kJ/mol */
    const scalar enthalpyH(
        218.0 
        + 20.78603*T_ 
        + 4.850638e-10*sqr(T_)/2.0 
        - 1.582916e-10*pow(T_, 3.0)/3.0 
        + 1.525102e-11*pow(T_, 4.0)/4.0 
        - 3.196347e-11/T_ + 211.8020 - 217.9994
    ); // kJ/mol

    return(enthalpyH / W_H_ * 1e6);
}


inline Foam::scalar Foam::H2::specificEnthalpyH2(scalar p, scalar T) const
{ 
    // Molecular hydrogen
    // const scalar n_(0.4986 + 1.1735*omega() + 0.4754*pow(omega(), 2.0));
    // const scalar alpha0_(0.42748*pow(R_*Tc(), 2.0) / Pc()); // J*m3/(mol^2)
    const scalar alpha_(alpha0_ / pow(T/Tc(), n_));
    // const scalar b_(0.08664*R_*Tc() / Pc()); // m3/mol
    // -> alpha/b = J/mol

    const scalar V(W()/densityH2(p, T)); // m3/mol

    const scalar specificEnthalpyH2_pg(cpPerfectGasH2(T) * T);
    return(
        (specificEnthalpyH2_pg + (p*V - R_*T - alpha_/b_*(1.0+n_)*log((V+b_)/V))/W())
    );
}


inline Foam::scalar Foam::H2::h(scalar p, scalar T) const
{
    const scalar w_H_(w_H(p, T));
    return(w_H_*specificEnthalpyH(T) + (1.0-w_H_)*specificEnthalpyH2(p, T));
}


inline Foam::scalar Foam::H2::cvPerfectGasH2(scalar T) const
{
    // Mayer's relation for perfect gas
    return cpPerfectGasH2(T) - R_/W();
}


inline Foam::scalar Foam::H2::cpPerfectGasH2(scalar T) const
{
    // Blending function
    const scalar sig_(sigmoid(T, 2.999842e-03, 1102.411));
    return
    (
        (1.0-sig_) * (4701.236*sigmoid(T, 2.326426e-02, 137.8745) + 9704.076)
        + sig_ * (-9.870353e-05*sqr(T) + 1.706997*T + 14170.57) 
    );
}


inline Foam::scalar Foam::H2::heatCapacityH(scalar T) const
{
    if (T < 298)
    {
        return 0.0;
    }

    // const scalar A_(20.78603);
    // const scalar B_(4.850638e-10);
    // const scalar C_(-1.582916e-10);
    // const scalar D_(1.525102e-11);
    // const scalar E_(3.196347e-11);

    const scalar T_(T/1000.0);

    const scalar cpNIST
    (
        // A_ + B_*T_ + C_*pow(T_, 2.0) + D_*pow(T_, 3.0) + E_/pow(T_, 2.0)
        20.78603 
        + 4.850638e-10*T_ 
        - 1.582916e-10*sqr(T_) 
        + 1.525102e-11*pow(T_, 3.0) 
        + 3.196347e-11/sqr(T_)
    );

    return(cpNIST / W_H_*1.0e3);
}


inline Foam::scalar Foam::H2::heatCapacityH2(scalar p, scalar T) const
{   
    const scalar Rmass(R_/W());

    // const scalar n_(0.4986 + 1.1735*omega() + 0.4754*sqr(omega()));
    const scalar alpha0_spec(alpha0_ / sqr(W())); // 0.42748*pow(Rmass*Tc(), 2.0) / Pc());
    const scalar alpha_(alpha0_spec / pow(T/Tc(), n_));
    const scalar b_spec(b_ / W()); // (0.08664*Rmass*Tc() / Pc());
    const scalar Vc_(Vc()*1e-3);
    const scalar c_(Rmass*Tc()/(Pc() + alpha0_spec/(Vc_ * (Vc_+b_spec))) + b_spec - Vc_);
    const scalar V(1.0/densityH2(p, T));

    const scalar dPdT(Rmass/(V - b_spec + c_) + n_*alpha_/(T*V*(V+b_spec)));
    const scalar dPdV
    (
        -Rmass*T / sqr(V-b_spec+c_) + alpha_*(b_spec+ 2.0*V) / sqr(V*(V+b_spec))
    );
    
    const scalar dVdT(-dPdT / dPdV);
    
    const scalar cpDeparture
    (
        -p*dVdT + Rmass - (n_*alpha_/T)*((1.0+n_)/b_spec)*log((V+b_spec)/V) 
        - alpha_*(1.0+n_)*dVdT / (V*(V+b_spec))
    );

    return(cpPerfectGasH2(T) - cpDeparture);
}


inline Foam::scalar Foam::H2::Cp(scalar p, scalar T) const
{
    // Mass weighted heat capacity for H-H2 mix
    const scalar w_H_(w_H(p, T));
    
    return(w_H_*heatCapacityH(T) + (1.0-w_H_)*heatCapacityH2(p, T));
}


inline Foam::scalar Foam::H2::B(scalar p, scalar T) const
{
    const scalar x_H_(x_H(T, p));
    const scalar W_H_H2_(x_H_*W_H_ + (1.0-x_H_)*W());

    if (T >= 50)
    {
        return(15.0 * 1e-6 / W_H_H2_);
    }

    const scalar B_
    (
        - 5.66667e-7 * pow(T, 6.0)
        + 1.26295e-4 * pow(T, 5.0)
        - 1.16747e-2 * pow(T, 4.0)
        + 0.578078 * pow(T, 3.0)
        - 16.37 * sqr(T)
        + 258.895 * T
        - 1933.3
    ); // cm3/mol

    return(B_*1e-6/W_H_H2_); // m3/kg
}


inline Foam::scalar Foam::H2::viscosity0(scalar T) const
{
    // const scalar epsilon_kb(30.41); // K
    // const scalar sigma(0.297); // nm

    const scalar T_star(T/30.41);

    const scalar S_star(exp(
        2.09630e-1
        - 4.55274e-1 * log(T_star)
        + 1.43602e-1 * sqr(log(T_star))
        - 3.35325e-2 * pow(log(T_star), 3.0)
        + 2.76981e-3 * pow(log(T_star), 4.0)
    ));

    return
    (
        // 0.021357 * sqrt(W()*T) / (pow(sigma, 2.0) * S_star) * 1e-6
        0.021357 * sqrt(W()*T) / (sqr(0.297) * S_star) * 1e-6
    );
}


inline Foam::scalar Foam::H2::viscosity1(scalar T) const
{
    // const scalar epsilon_kb(30.41); // K
    // const scalar sigma(0.297); // nm

    const scalar T_star(T/30.41);

    const scalar B_star(
        -0.1870 
        + 2.4871/T_star 
        + 3.7151/sqr(T_star)
        -11.0972/pow(T_star, 3.0)
        + 9.0965/pow(T_star, 4.0)
        - 3.8292/pow(T_star, 5.0)
        + 0.5166/pow(T_star, 6.0)
    );

    return
    (
        B_star * pow(0.297, 3.0) * viscosity0(T)
    );
}


inline Foam::scalar Foam::H2::viscosityH(scalar T) const
{
    return 2.6693e-6 * sqrt(W_H_*T) / omega22_22(T);
}


inline Foam::scalar Foam::H2::mu(scalar p, scalar T) const
{
    // Molecular hydrogen 
    // const scalar C1(6.43449673);
    // const scalar C2(4.56334068e-2);
    // const scalar C3(2.32797868e-1);
    // const scalar C4(9.58326120e-1);
    // const scalar C5(1.27941189e-1);
    // const scalar C6(3.63576595e-1);

    // const scalar rho_sc(90.5); // kg/m3
    const scalar rhoH2(rho(p, T));
    const scalar rhoR(rhoH2/90.5);
    const scalar Tr(T/Tc());
    /*const scalar viscosityH2_(
        viscosity0(T) + viscosity1(T)*rhoH2 + 1e-6*(
            C1*sqr(rhoR) * exp(
                C2*Tr + C3/Tr + C4*sqr(rhoR)/(C5+Tr) + C6*pow(rhoR, 6.0)
            )
        )
    );*/
    const scalar viscosityH2_(
        viscosity0(T) + viscosity1(T)*rhoH2 + 1e-6*(
            6.43449673*sqr(rhoR) * exp(
                4.56334068e-2*Tr 
                + 2.32797868e-1/Tr 
                + 9.58326120e-1*sqr(rhoR)/(1.27941189e-1+Tr) 
                + 3.63576595e-1*pow(rhoR, 6.0)
            )
        )
    );

    // Atomic hydrogen 
    const scalar viscosityH_(viscosityH(T));

    // Mix
    const scalar x_H_(x_H(p, T));

    // Pure molecular hydrogen
    if (x_H_ == 0) 
    {
        return viscosityH2_;
    }

    const scalar D_H_H2(D(p, T));
    const scalar A_(A_12(T));

    const scalar H_11(
        sqr(1.0-x_H_)/viscosityH2_
        + 2.0*(1.0-x_H_)*x_H_/(W()+W_H_) * R_*T/(p*D_H_H2) * (1.0 + 3.0*W_H_*A_/(5.0*W()))
    );
    const scalar H_22(
        sqr(x_H_)/viscosityH_ 
        + 2.0*(1.0-x_H_)*x_H_/(W()+W_H_) * R_*T/(p*D_H_H2) * (1.0 + 3.0*W()*A_/(5.0*W_H_))
    );
    const scalar H_12(
        -2.0*(1.0-x_H_)*x_H_/(W()+W_H_) * R_*T/(p*D_H_H2) * (1.0 - 3.0*A_/5.0)
    );

    return(
        (sqr(1.0-x_H_)/H_11 + sqr(x_H_)/H_22 - 2.0*(1.0-x_H_)*x_H_*H_12/(H_11*H_22)) 
        / (1.0 - sqr(H_12)/(H_11*H_22))
    );
}


inline Foam::scalar Foam::H2::thermalConductivityH2(scalar p, scalar T) const
{
    // Reduced temperature (divided by critical temperature)
    const scalar T_(T/33.45);

    const scalar lambda_0_num
    (
        -3.40976e-1
        +4.58820 * T_
        -1.45080 * sqr(T_)
        +3.26394e-1 * pow(T_, 3.0)
        +3.16939e-3 * pow(T_, 4.0)
        +1.90592e-4 * pow(T_, 5.0)
        -1.139e-6 * pow(T_, 6.0)
    );
    const scalar lambda_0_denom
    (
        1.38497e2
        - 2.21878e1 * T_
        + 4.57151 * sqr(T_)
        + 1.0 * pow(T_, 3.0)
    );
    const scalar lambda_0(lambda_0_num/lambda_0_denom);

    // Molecular hydrogen
    const scalar rhoH2(rho(p, T));
    const scalar rho_c(31.262); // kg/m3

    const scalar rho_r(rhoH2/rho_c);
    const scalar lambda_excess
    (
          ( 3.63081e-2 + 1.83370e-3*T_) * rho_r
        + (-2.07629e-2 - 8.86716e-3*T_) * sqr(rho_r)
        + ( 3.14810e-2 + 1.58260e-2*T_) * pow(rho_r, 3.0)
        + (-1.43097e-2 - 1.06283e-2*T_) * pow(rho_r, 4.0)
        + ( 1.74980e-3 + 2.80673e-3*T_) * pow(rho_r, 5.0)
    );

    // const scalar c1(6.24e-4); // W/m/K
    // const scalar c2(-2.58e-7);
    // const scalar c3(0.837);

    // Empirical, a theoretical model exists but more complex
    const scalar lambda_critical
    (
        // c1/(c2 + T_ -1.0) * exp(-pow(c3*(rho_r-1.0), 2.0))
        6.24e-4/(-2.58e-7 + T_ -1.0) * exp(-sqr(0.837*(rho_r-1.0)))
    );

    return lambda_0 + lambda_excess + lambda_critical;
}


inline Foam::scalar Foam::H2::kappa(scalar p, scalar T) const
{
    const scalar x_H_(x_H(p, T));
    const scalar W_H_H2_(x_H_*W_H_ + (1.0-x_H_)*W());
    const scalar A12(A_12(T));
    const scalar B12(B_12(T));
    const scalar M_12(2.0*W_H_*W()/(W_H_+W()));
    const scalar lambdaH2(thermalConductivityH2(p, T));

    const scalar lambda_H(15.0/4.0 * R_*viscosityH(T)/W_H_);

    const scalar lambda_12(15.0/4.0 * R_/W_H_H2_ * 2.6693e-4*sqrt(M_12*1e-3*T)/omega22_12(T));
    const scalar U_y
    (
        4.0/15.0 * A12 * W_H_*W()/sqr(M_12) * sqr(lambda_12)/(lambda_H*lambdaH2) 
        - (B12/5.0+1.0/12.0) 
        - 1.0/A12*(3.0*B12/8.0-25.0/32.0) * sqr(W_H_-W())/(W_H_*W())
    );
    const scalar U_z
    (
        4.0/15.0 * A12 * (
            lambda_12*(1.0/lambda_H + 1.0/lambdaH2)*W_H_*W()/sqr(M_12) -1.0
        ) 
        - (B12/5.0+1.0/12.0)
    );
    const scalar U_1
    (
        4.0/15.0 * A12 - (B12/5.0+1.0/12.0)*W_H_/W() 
        + sqr(W_H_-W())/(2.0*W_H_*W())
    );
    const scalar U_2
    (
        4.0/15.0 * A12 - (B12/5.0+1.0/12.0)*W()/W_H_ 
        + sqr(W()-W_H_)/(2.0*W_H_*W())
    );

    const scalar X_lambda
    (
        sqr(x_H_)/lambda_H + 2.0*x_H_*(1.0-x_H_)/lambda_12 
        + sqr(1.0-x_H_)/lambdaH2
    );
    const scalar Y_lambda
    (
        sqr(x_H_)/lambda_H*U_1 + 2.0*x_H_*(1.0-x_H_)/lambda_12*U_y 
        + sqr(1.0-x_H_)/lambdaH2*U_2
    );
    const scalar Z_lambda
    (
        sqr(x_H_)*U_1 + 2.0*x_H_*(1.0-x_H_)*U_z + sqr(1.0-x_H_)*U_2
    );

    const scalar lambda_f((1.0+Z_lambda)/(X_lambda+Y_lambda));

    const scalar delta_H
    (
        2.0*W_H_*specificEnthalpyH(T) - W()*specificEnthalpyH2(p, T)
    );

    const scalar lambda_r
    (
        p*D(p, T)/T * sqr(delta_H/(R_*T)) 
        * x_H_*(1.0-x_H_)/sqr(2.0-x_H_)
    );

    return lambda_f + lambda_r;
}


inline Foam::scalar Foam::H2::diffusionH2(scalar p, scalar T) const
{
    // Convert pressure from Pa to atm
    // const scalar p_(p/101325.0);

    return 2.628e-7 / (p/101325.0*omega11_11(T)) * sqrt(pow(T, 3.0)/(W()*1.0e3));
}


inline Foam::scalar Foam::H2::D(scalar p, scalar T) const
{
    // const scalar p_(p/101325.0);
    const scalar M_12(2.0*W_H_*W()/(W_H_+W()));

    return 2.628e-7 / (p/101325.0*omega11_12(T)) * sqrt(pow(T, 3.0)/M_12);
}


inline Foam::scalar Foam::H2::alphah(const scalar p, const scalar T) const
{
    return kappa(p, T)/Cp(p, T);
}


inline Foam::scalar Foam::H2::Pr(const scalar p, const scalar T) const
{
    return Cp(p, T) * mu(p, T) / kappa(p, T);
}


inline Foam::scalar Foam::H2::sigmoid(scalar T, scalar a, scalar b) const
{
    return 1.0 / (1.0 + exp( -a * (T - b) ));
}


// ************************************************************************* //
