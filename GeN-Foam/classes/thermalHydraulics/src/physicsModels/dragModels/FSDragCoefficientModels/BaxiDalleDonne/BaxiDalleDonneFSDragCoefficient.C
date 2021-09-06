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
#include "BaxiDalleDonneFSDragCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FSDragCoefficientModels
{
    defineTypeNameAndDebug(BaxiDalleDonne, 0);
    addToRunTimeSelectionTable
    (
        FSDragCoefficientModel, 
        BaxiDalleDonne, 
        FSDragCoefficientModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSDragCoefficientModels::BaxiDalleDonne::BaxiDalleDonne
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
    ),
    A_(0),
    B_(0),
    C_(0)
{
    scalar Dp(dict.get<scalar>("pinDiameter"));
    scalar Dw(dict.get<scalar>("wireDiameter"));
    scalar H(dict.get<scalar>("wireLeadLen"));
    scalar Pt(Dp+1.0444*Dw);
    A_ = (80.0/sqrt(100*H))*pow(Pt/Dp, 1.5); //- H provided in m but needs to
                                             //  be in cm
    B_ = 1.034/pow(Pt/Dp, 0.124);
    C_ = 29.7*pow(Pt/Dp, 6.9)/pow(H/(Dp+Dw), 2.239);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::FSDragCoefficientModels::BaxiDalleDonne::value
(
    const label& celli
) const
{
    const scalar& Rei(Re(celli));
    scalar fl(A_/Rei);  //- There should be a *pow(Twall/Tbulk, 1.5), but I am
                        //  lazy now (plus, what about cases where I do not 
                        //  solve for energy? I could handle this with some
                        //  flags and ifs)
    scalar ft(0.316*(B_+C_*pow(Rei,0.086))/pow(Rei, 0.25));
    if (Rei > 400)
    {
        if (Rei < 5000)
        {
            scalar psi
            (
                min
                (
                    max
                    (
                        ((Rei-400.0)/4600.0), 
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
