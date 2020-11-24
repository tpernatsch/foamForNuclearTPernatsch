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

#include "NoKazimiFS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(NoKazimiFS, 0);
    addToRunTimeSelectionTable(dragModel, NoKazimiFS, FSDragModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::NoKazimiFS::NoKazimiFS
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FSPair& FSPair
)
:
    dragModel
    (
        objReg,
        dict,
        FSPair
    ),
    A_(0),
    B_(0),
    C_(0)
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

Foam::scalar Foam::dragModels::NoKazimiFS::fd(const scalar& Re) const
{
    scalar fl(A_/Re);   //- There should be a *pow(Twall/Tbulk, 1.5), but yeah,
                        //  go access those fields...
    scalar ft(0.316*(B_+C_*pow(Re,0.086))/pow(Re, 0.25));
    scalar fd(fl);
    if (Re >= 400)
    {
        if (Re < 2600)
        {
            scalar psi
            (
                max
                (
                    min
                    (
                        ((Re-400.0)/2200.0), 
                        0.0
                    ), 
                    1.0
                )
            );
            fd = sqrt(psi)*ft+sqrt(1.0-psi)*fl;
        }
        else    fd = ft;
    }
    return fd;
}

void Foam::dragModels::NoKazimiFS::correctKd(volTensorField& Kd) const
{    
    #include "calcKdFromFd.H"
}


// ************************************************************************* //
