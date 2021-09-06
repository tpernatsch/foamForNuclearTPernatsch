/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "FSPair.H"
#include "NoKazimiFSDragCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FSDragCoefficientModels
{
    defineTypeNameAndDebug(NoKazimi, 0);
    addToRunTimeSelectionTable
    (
        FSDragCoefficientModel, 
        NoKazimi, 
        FSDragCoefficientModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSDragCoefficientModels::NoKazimi::NoKazimi
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    FSDragCoefficientModel
    (
        pair,
        dict,
        objReg
    )
{
    scalar Dp(dict.get<scalar>("pinDiameter"));
    scalar Dw(dict.get<scalar>("wireDiameter"));
    scalar H(dict.get<scalar>("wireLeadLen"));
    scalar Pt(Dp+1.0444*Dw);
    A_ = (32/sqrt(100.0*H))*pow(Pt/Dp, 1.5);
    B_ = 1.034/pow(Pt/Dp, 0.124);
    C_ = 29.7*pow(Pt/Dp, 6.9)/pow(H/(Dp+Dw), 2.239);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::FSDragCoefficientModels::NoKazimi::value
(
    const label& celli
) const
{
    const scalar& Rei(Re(celli));
    scalar fl(A_/Rei);              //- I am missing a *pow(Twall/Tbulk, 1.5)
                                    //  for simplicity, as Twall is 0 in those
                                    //  structure regions where there are no
                                    //  powerModels or passiveModels specified
                                    //  I could handle it via additional ifs,
                                    //  but eh
    scalar ft(0.316*(B_+C_*pow(Rei,0.086))/pow(Rei, 0.25));
    if (Rei > 400)
    {
        if (Rei < 2600)
        {
            scalar psi
            (
                min
                (
                    max
                    (
                        ((Rei-400.0)/2200.0), 
                        0.0
                    ), 
                    1.0
                )
            );
            return sqrt(psi)*ft+sqrt(1.0-psi)*fl;
        }
        else    
            return ft;
    }
    return fl;
}

// ************************************************************************* //
