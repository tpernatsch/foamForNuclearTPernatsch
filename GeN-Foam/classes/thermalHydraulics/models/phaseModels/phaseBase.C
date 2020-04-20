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

#include "phaseBase.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseBase::phaseBase
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& name,
    bool readIfPresFlag,
    bool writeFlag
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", name),
            mesh.time().timeName(),
            mesh,
            (
                (readIfPresFlag) ?
                IOobject::READ_IF_PRESENT :
                IOobject::NO_READ
            ),
            (
                (writeFlag) ?
                IOobject::AUTO_WRITE :
                IOobject::NO_WRITE
            )
        ),
        mesh,
        dimensionedScalar("alpha", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    dict_
    (
        dict
    ),
    mesh_(mesh),
    name_(name),
    residualAlpha_
    (
        dimensionedScalar::lookupOrDefault
        (
            "residualAlpha", 
            dict, 
            dimless, 
            1e-6
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseBase::~phaseBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
