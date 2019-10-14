/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "BrowningPotter.H"
#include "fluid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(BrowningPotter, 0);
    addToRunTimeSelectionTable
    (
        saturationModel,
        BrowningPotter,
        saturationModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::BrowningPotter::
BrowningPotter
(
    const dictionary& dict, 
    const fluid& fluid1,
    const fluid& fluid2
)
:
    saturationModel
    (
        dict,
        fluid1,
        fluid2
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::saturationModels::BrowningPotter::pSat
(
    const volScalarField& T
) const
{
    return onePa_*Foam::exp(lnPSat(T));
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::BrowningPotter::pSatPrime
(
    const volScalarField& T
) const
{
    return tmp<volScalarField>
    (
        pSat(T)*(-0.4672/T + oneK_*12633.37/sqr(T))
    );
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::BrowningPotter::lnPSat
(
    const volScalarField& T
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
                (11.9463 - 12633.37/T - 0.4672*Foam::log(T/oneK_))
            +   Foam::log(1e6) // To have the ln of P in Pa rather than MPa
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::BrowningPotter::Tsat
(
    const volScalarField& p
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
                oneK_
            *
                923840.0
            /   (
                    -   11275 
                    +   Foam::sqrt
                    (
                        127125625 + 1847680 * 
                        (
                                7.8270 
                            -   Foam::log(p/1000000/onePa_) // p in MPa
                        ) 
                    )
                )
        )
    );
}


// ************************************************************************* //
