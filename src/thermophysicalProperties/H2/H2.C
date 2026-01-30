/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2512                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.          |
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

#include "H2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(H2, 0);
    addToRunTimeSelectionTable(liquidProperties, H2,);
    addToRunTimeSelectionTable(liquidProperties, H2, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::H2::H2()
:
    liquidProperties
    (
        2.01588, // kg/kmol
        33.145, // K
        1.2965e6, // Pa
        6.4481e-2, // m3/kmol
        0.305, // -
        13.81, // K
        7.042e3, // Pa
        20.28, // K
        0.0, // C.m
        -0.219, // -
         5.08e3 // Pa^(1/2) // https://www.hindawi.com/journals/jchem/2016/4701919/
    ),
    W_H_(1.007647), // kg/kmol
    R_(8.314472e3), // J/kmol/K
    n_(0.4986 + 1.1735*omega() + 0.4754*sqr(omega())),
    alpha0_(0.42748*sqr(R_*Tc()) / Pc()),
    b_(0.08664*R_*Tc() / Pc()),
    c_(R_*Tc() / (Pc() + alpha0_ / (Vc() * (Vc()+b_))) + b_ - Vc()),
    // Arbitrary values not used in the code
    rho_(98.343885, 0.30542, 647.13, 0.081),
    pv_(73.649, -7258.2, -7.3037, 4.1653e-06, 2),
    hl_(647.13, 2889425.47876769, 0.3199, -0.212, 0.25795, 0),
    Cp_
    (
        15341.1046350264,
        -116.019983347211,
        0.451013044684985,
        -0.000783569247849015,
        5.20127671384957e-07,
        0
    ),
    h_
    (
        -17957283.7993676,
        15341.1046350264,
        -58.0099916736053,
        0.150337681561662,
        -0.000195892311962254,
        1.04025534276991e-07
    ),
    Cpg_
    (
        1851.73466555648,
        1487.53816264224,
        2609.3,
        493.366638912018,
        1167.6
    ),
    B_
    (
        -0.0012789342214821,
        1.4909797391063,
        -1563696.91923397,
        1.85445462114904e+19,
        -7.68082153760755e+21
    ),
    mu_(-51.964, 3670.6, 5.7331, -5.3495e-29, 10),
    mug_(2.6986e-06, 0.498, 1257.7, -19570),
    kappa_(-0.4267, 0.0056903, -8.0065e-06, 1.815e-09, 0, 0),
    kappag_(6.977e-05, 1.1243, 844.9, -148850),
    sigma_(647.13, 0.18548, 2.717, -3.554, 2.047, 0),
    D_(15.0, 15.0, 18.015, 28)
{}


Foam::H2::H2
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
)
:
    liquidProperties(l),
    rho_(density),
    pv_(vapourPressure),
    hl_(heatOfVapourisation),
    Cp_(heatCapacity),
    h_(enthalpy),
    Cpg_(idealGasHeatCapacity),
    B_(secondVirialCoeff),
    mu_(dynamicViscosity),
    mug_(vapourDynamicViscosity),
    kappa_(thermalConductivity),
    kappag_(vapourThermalConductivity),
    sigma_(surfaceTension),
    D_(vapourDiffussivity)
{}


Foam::H2::H2(const dictionary& dict)
:
    H2()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::H2::writeData(Ostream& os) const
{
    liquidProperties::writeData(os);
    os  << token::SPACE
        << Tc() << token::SPACE
        << Pc() << token::SPACE
        << Vc() << token::SPACE
        << Zc() << token::SPACE
        << Tt() << token::SPACE
        << Pt() << token::SPACE
        << Tb() << token::SPACE
        << dipm() << token::SPACE
        << omega() << token::SPACE
        << delta();
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const H2& l)
{
    l.writeData(os);
    return os;
}


// ************************************************************************* //
