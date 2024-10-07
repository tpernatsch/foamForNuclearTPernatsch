/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2406                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2024 OpenCFD Ltd.          |
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

Class
    Foam::H2

Description
    Hydrogen thermophysical properties. Takes into account the dissociation of
    molecular hydrogen at high temperature.

    Author of this file (and associated .C or included .H files):
    Thomas Guilbaud <thomas.guilbaud@epfl.ch>, EPFL (Switzerland)
    Eymeric Simonnot <eymeric.simonnot@epfl.ch>, EPFL (Switzerland)

References
    [1] Fang et al, Study on high-temperature hydrogen dissociation for nuclear
        thermal propulsion reactor, Nuclear Engineering and Design 392 (2022) 
        111753

    [2] Reisfeld, THERMODYNAMIC AND TRANSPORT PROPERTIES OF DISSOCIATED 
        HYDROGEN MIXTURES, 1957, LANL, LA-2123.

    [3] Roder et al, Computer Programs for Thermodynamic and Transport 
        Properties of Hydrogen (Tabcode-ll), 1972, NBS Technical note 625.

    [4] NIST Chemistry WebBook, Hydrogen atom - Gas phase thermochemestry data: 
        https://webbook.nist.gov/cgi/cbook.cgi?ID=C12385136&Mask=1

    [5] Muzny, C.D., Huber, M.L., Kazakov, A.F., 2013. Correlation for the 
        viscosity of normal hydrogen obtained from symbolic regression. 
        J. Chem. Eng. Data 58 (4), 969–979. https://doi.org/10.1021/je301273j

    [6] Vanderslice, J.T., Weissman, S., Mason, E.A., Fallon, R.J., 1962. 
        High-temperature transport properties of dissociating hydrogen. 
        Phys. Fluids 5 (2), 155–164. https://doi.org/10.1063/1.1706590

    [7] NIST Chemistry WebBook, Hydrogen molecule - Gas phase thermochemestry 
        data: https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=1

    [8] Assael, M.J., Assael, J.-A.- M., Huber, M.L., Perkins, R.A., 
        Takata, Y., 2011. Correlation of the Thermal Conductivity of Normal and
        Parahydrogen from the Triple Point to 1000 K and up to 100 MPa. 
        J. Phys. Chem. Ref. Data 40 (3), 033101. 
        https://doi.org/10.1063/1.3606499

    [9] Leachman, J. W., R. T Jacobsen, S. G. Penoncello, et E. W. Lemmon. 
        « Fundamental Equations of State for Parahydrogen, Normal Hydrogen, and 
        Orthohydrogen ». Journal of Physical and Chemical Reference Data 38, 
        no 3 (1 septembre 2009): 721‑48. https://doi.org/10.1063/1.3160306.

    [10] Robert D.McCarty, Hydrogen Technological Survey, Thermophysical 
        Properties Cryogenics Division, Institute for Basic Standards
        https://ntrs.nasa.gov/api/citations/19760004209/downloads/19760004209.pdf

    [11] Mason, Edward A., et William E. Rice. « The Intermolecular Potentials 
        of Helium and Hydrogen ». The Journal of Chemical Physics 22, no 3 
        (1 mars 1954): 522‑35. https://doi.org/10.1063/1.1740100.

    [12] Mason, Edward A. « Transport Properties of Gases Obeying a Modified 
        Buckingham (Exp-Six) Potential ». The Journal of Chemical Physics 22, 
        no 2 (1 février 1954): 169‑86. https://doi.org/10.1063/1.1740026.

    [13] Amdur, I., et E. A. Mason. « Properties of Gases at Very High 
        Temperatures ». Physics of Fluids 1, no 5 (1958): 370. 
        https://doi.org/10.1063/1.1724353.

SourceFiles
    H2.C

\*---------------------------------------------------------------------------*/

#ifndef H2_H
#define H2_H

#include "liquidProperties.H"
#include "NSRDSfunc0.H"
#include "NSRDSfunc1.H"
#include "NSRDSfunc2.H"
#include "NSRDSfunc3.H"
#include "NSRDSfunc4.H"
#include "NSRDSfunc5.H"
#include "NSRDSfunc6.H"
#include "NSRDSfunc7.H"
#include "APIdiffCoefFunc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class H2 Declaration
\*---------------------------------------------------------------------------*/

class H2
:
    public liquidProperties
{
    // Private data

        //- Hydrogen atomic weight [kg/kmol] 
        scalar W_H_;

        //- Perfect gas constant
        scalar R_;

        //- 
        scalar n_;

        //-
        scalar alpha0_;

        //-
        scalar b_;

        //-
        scalar c_;


        //- Useless but needed because liquidProperties is over constrainded
        NSRDSfunc5 rho_;
        NSRDSfunc1 pv_;
        NSRDSfunc6 hl_;
        NSRDSfunc0 Cp_;
        NSRDSfunc0 h_;
        NSRDSfunc7 Cpg_;
        NSRDSfunc4 B_;
        NSRDSfunc1 mu_;
        NSRDSfunc2 mug_;
        NSRDSfunc0 kappa_;
        NSRDSfunc2 kappag_;
        NSRDSfunc6 sigma_;
        APIdiffCoefFunc D_;

public:

    friend class liquidProperties;

    //- Runtime type information
    TypeName("H2");


    // Constructors

        //- Construct null
        H2();

        //- Construct from components
        H2
        (
            const liquidProperties& l,
            const NSRDSfunc5& density,
            const NSRDSfunc1& vapourPressure,
            const NSRDSfunc6& heatOfVapourisation,
            const NSRDSfunc0& heatCapacity,
            const NSRDSfunc0& enthalpy,
            const NSRDSfunc7& idealGasHeatCapacity,
            const NSRDSfunc4& secondVirialCoeff,
            const NSRDSfunc1& dynamicViscosity,
            const NSRDSfunc2& vapourDynamicViscosity,
            const NSRDSfunc0& thermalConductivity,
            const NSRDSfunc2& vapourThermalConductivity,
            const NSRDSfunc6& surfaceTension,
            const APIdiffCoefFunc& vapourDiffussivity
        );

        //- Construct from dictionary
        H2(const dictionary& dict);

        //- Construct and return clone
        virtual autoPtr<liquidProperties> clone() const
        {
            return autoPtr<liquidProperties>::NewFrom<H2>(*this);
        }


    // Member Functions

        // Hydrogen dissociation

            //- Equilibrium constant from [1] Eq 3 and [2]
            inline scalar Kp(scalar T) const;
            
            //- Mole fraction of atomic hydrogen from [1] Eq 2 and [2]
            inline scalar x_H(scalar p, scalar T) const;
            
            //- Weight fraction of atomic hydrogen from [1] Eq 2 and [2]
            inline scalar w_H(scalar p, scalar T) const;


        // Average effective collision integrals, definition of the average 
        // effective collision integrals defined in [1] and [6]

            //- Ω(1,1){1,1} defined in [1]
            //  Coeffs fitted with data from [6], table II
            //  Fitted on values from 1000 to 10000 K every 500 K and extended
            //  at cryogenic temperatures
            inline scalar omega11_11(scalar T) const;

            //- Ω(2,2){1,1} defined in [1]
            //  Coeffs fitted with data from [6], table II, for T>1000K
            //  Coeffs fitted with data from [12], using parameters from [6], 
            //  [11], and Eq.(23) from [13]
            inline scalar omega22_11(scalar T) const;

            // inline scalar omega11_22(scalar T) const;
            
            //- Ω(2,2){2,2} defined in [1]
            //  Coeffs fitted with data from [6], table II
            inline scalar omega22_22(scalar T) const;

            //- Ω(1,1){1,2} defined in [1]
            //  Coeffs fitted with data from [6], table II
            //  Fitted on values from 1000 to 10000 K every 500 K and extended 
            //  to cryogenic temperatures
            inline scalar omega11_12(scalar T) const;
            
            //- Ω(2,2){1,2} defined in [1]
            //  Coeffs fitted with data from [6], table II
            //  Fitted on values from 1000 to 10000 K every 500 K and extended 
            //  to cryogenic temperatures
            inline scalar omega22_12(scalar T) const;

            //- A{1,2} defined in [1]
            //  Coeffs fitted with data from [6], table II
            //  Fitted on values from 1000 to 10000 K every 500 K and extended
            //  to cryogenic temperatures
            inline scalar A_12(scalar T) const;

            //- B{1,2} defined in [1]
            //  Coeffs fitted with data from [6], table II
            //  Fitted on values from 1000 to 10000 K every 500 K and extended
            // to cryogenic temperatures
            inline scalar B_12(scalar T) const;
            
            //- C{1,2} defined in [1]
            //  Coeffs fitted with data from [6], table II
            //  Fitted on values from 1000 to 10000 K every 500 K and extended 
            //  to cryogenic temperatures
            inline scalar C_12(scalar T) const;


        // Fundamental equation of state properties

            //- Density of atomic hydrogen [kg/m^3]
            inline scalar densityH(scalar p, scalar T) const;
            
            //- Density of molecular hydrogen [kg/m^3]
            inline scalar densityH2(scalar p, scalar T) const;
            
            //- Density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Liquid compressibility rho/p [s^2/m^2]
            //  Note: currently it is assumed the liquid is incompressible
            inline scalar psi(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K)]
            inline scalar CpMCv(scalar p, scalar T) const;


        // Fundamental thermodynamic properties

            //- Sensible enthalpy [J/kg]
            inline scalar Hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- Absolute Enthalpy [J/kg]
            inline scalar Ha(const scalar p, const scalar T) const;

            //- Ideal gas entropy of H2 in [J/kg/K]
            //  Entropy ref from NIST and https://fr.wikipedia.org/wiki/Entropie_d%27un_corps_pur
            inline scalar idealEntropyH2(const scalar p, const scalar T) const;

            //- Ideal gas entropy of H in [J/kg/K]
            //  https://fr.wikipedia.org/wiki/Entropie_d%27un_corps_pur
            inline scalar idealEntropyH(const scalar p, const scalar T) const;

            //- Real gas entropy of H2 in [J/kg/K]
            //  https://www.afs.enea.it/project/neptunius/docs/fluent/html/ug/node335.htm for entropy departure
            //  From NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=1//Thermo-Gas
            inline scalar entropyH2(const scalar p, const scalar T) const;
            
            //- Entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Specific enthalpy of atomic hydrogen [J/kg]
            //  From NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C12385136&Mask=1//Thermo-Gas
            inline scalar specificEnthalpyH(scalar T) const;

            //- Specific enthalpy of molecular hydrogen [J/kg]
            inline scalar specificEnthalpyH2(scalar p, scalar T) const;

            //- Liquid Enthalpy [J/(kg)]
            inline scalar h(scalar p, scalar T) const;

            //- Ideal gas isochoric heat capacity of molecular hydrogen [J/(kg K)]
            inline scalar cvPerfectGasH2(scalar T) const;

            //- Ideal gas isobaric heat capacity of molecular hydrogen [J/(kg K)]
            inline scalar cpPerfectGasH2(scalar T) const;

            //- Heat capacity of atomic hydrogen [J/(kg K)]
            //  Coefficients from [4] for T in [298; 6000] K for monoatomic 
            //  hydrogen heat capacity
            //  Note: the equation of State for monoatomic hydrogen is 
            //  PerfectGas, hence its capacity depends only on temperature
            inline scalar heatCapacityH(scalar T) const;

            //- Heat capacity of molecular hydrogen [J/(kg K)]
            //  Coefficients from [7] for T in [298; 6000] K
            inline scalar heatCapacityH2(scalar p, scalar T) const;

            //- Liquid heat capacity [J/(kg K)]
            inline scalar Cp(scalar p, scalar T) const;

            //- Second Virial Coefficient [m^3/kg]
            inline scalar B(scalar p, scalar T) const;


        // Physical properties

            //- Coefficients from [5] zero_density_limit Eqs. (3) and (4), 
            //  coefficients ai in Table 2. 
            inline scalar viscosity0(scalar T) const;

            //- Coefficients from [5] the excess contribution Eqs. (5), (6) and 
            //  (7), coefficients bi in Table 3.
            inline scalar viscosity1(scalar T) const;

            //- Viscosity of atomic hydrogen [Pa s]
            inline scalar viscosityH(scalar T) const;

            //- Viscosity [Pa s]
            //  Coefficients from [5] Eq. (9), coefficients Ci in Table 4.
            inline scalar mu(scalar p, scalar T) const;

            //- Thermal conductivity of molecular hydrogen [W/(m K)]
            // Expressions from:
            // - Ref. [8] 
            // - lambda_0 expression from Eq. (2), coefficients Aij found in Table 2.
            // - lambda_excess expression from Eq. (3), coefficients Bij found in Table 2.
            // - lambda_critical expression from Eq. (8), coefficients Ci found in the text following Eq. (8).
            // Note: lambda_critical used here is obtained via empirical 
            // correlation, the deviation for the Roder and Diller data is -8% 
            // to +15%, while those for the simplified Olchowy-Sengers 
            // enhancement vary from -9% to +9%.
            inline scalar thermalConductivityH2(scalar p, scalar T) const;
            
            //- Thermal conductivity  [W/(m K)]
            //  Expressions from [1]
            inline scalar kappa(scalar p, scalar T) const;

            //- Vapour diffusion of molecular hydrogen [m2/s]
            //  Expressions from [1]
            inline scalar diffusionH2(scalar p, scalar T) const;

            //- Vapour diffusivity [m2/s]
            //  Expressions from [1]
            inline scalar D(scalar p, scalar T) const;

            //- Liquid thermal diffusivity of enthalpy [kg/ms]
            inline scalar alphah(const scalar p, const scalar T) const;

            //- Prandtl number [-]
            inline scalar Pr(const scalar p, const scalar T) const;

                    
            //- Vapour pressure [Pa]
            inline scalar pv(scalar p, scalar T) const { return 0; }
                                    
            //- Heat of vapourisation [J/kg]
            inline scalar hl(scalar p, scalar T) const { return h(p, T); }
                                    
            //- Ideal gas heat capacity [J/(kg K)]
            inline scalar Cpg(scalar p, scalar T) const { return Cp(p, T); }
                                                                            
            //- Vapour viscosity [Pa s]
            inline scalar mug(scalar p, scalar T) const { return mu(p, T); }

            //- Vapour thermal conductivity [W/(m K)]      
            inline scalar kappag(scalar p, scalar T) const { return kappa(p, T); }
                                                                
            //- Surface tension [N/m]
            inline scalar sigma(scalar p, scalar T) const { return 0; }

            //- Vapour diffusivity [m2/s]
            inline scalar D(scalar p, scalar T, scalar Wb) const { return D(p, T); }


        // Functions

            //- Sigmoid function
            inline scalar sigmoid(scalar T, scalar a, scalar b) const;

    // I-O

        //- Write the function coefficients
        void writeData(Ostream& os) const;

        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const H2& l);
};


Ostream& operator<<(Ostream& os, const H2& l);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "H2I.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
