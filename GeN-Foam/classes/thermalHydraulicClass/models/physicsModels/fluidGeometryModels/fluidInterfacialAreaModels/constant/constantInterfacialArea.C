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

#include "constantInterfacialArea.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidInterfacialAreaModels
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable
    (
        fluidInterfacialAreaModel, 
        constant, 
        fluidInterfacialAreaModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidInterfacialAreaModels::constant::constant
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const fluid& fluid1,
    const fluid& fluid2
)
:
    fluidInterfacialAreaModel
    (
        objReg,
        dict,
        fluid1,
        fluid2
    ),
    iA_
    (
        IOobject
        (
            "",
            fluid1.mesh().time().timeName(),
            fluid1.mesh()
        ),
        fluid1.mesh(),
        dimensionedScalar("value", dimArea/dimVol, dict)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidInterfacialAreaModels::constant::~constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::fluidInterfacialAreaModels::constant::iA() const
{
    tmp<volScalarField> tiA
    (
        new volScalarField
        (
            iA_
        )
    );
    
    return tiA;
}


// ************************************************************************* //
