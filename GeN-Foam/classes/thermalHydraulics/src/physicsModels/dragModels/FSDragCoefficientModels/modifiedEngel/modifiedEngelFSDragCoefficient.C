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
#include "modifiedEngelFSDragCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FSDragCoefficientModels
{
    defineTypeNameAndDebug(modifiedEngel, 0);
    addToRunTimeSelectionTable
    (
        FSDragCoefficientModel, 
        modifiedEngel, 
        FSDragCoefficientModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSDragCoefficientModels::modifiedEngel::modifiedEngel
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
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::FSDragCoefficientModels::modifiedEngel::value
(
    const label& celli
) const
{
    const scalar& Rei(Re(celli));
    scalar fl(110.0/Rei);
    scalar ft(0.37/pow(Rei, 0.25));
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
