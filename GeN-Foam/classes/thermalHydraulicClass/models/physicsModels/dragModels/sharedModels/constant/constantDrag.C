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

#include "constantDrag.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable(dragModel, constant, FFDragModels);
    addToRunTimeSelectionTable(dragModel, constant, FSDragModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::constant::constant
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FFPair& FFPair
)
:
    dragModel
    (
        objReg,
        dict,
        FFPair
    ),
    Kd_
    (
        IOobject
        (
            "Kd."+this->pairName_,
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedTensor
        (
            "tensorKd", 
            this->dimK, 
            tensor
            (
                dict.lookupType<scalar>("value"), 0, 0,
                0, dict.lookupType<scalar>("value"), 0,
                0, 0, dict.lookupType<scalar>("value")
            )
        )
    )
{}

Foam::dragModels::constant::constant
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FSPair& FSPair,
    const word region
)
:
    dragModel
    (
        objReg,
        dict,
        FSPair,
        region
    ),
    Kd_
    (
        IOobject
        (
            "Kd."+this->pairName_,
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedTensor("value", this->dimK, dict)
    )
{
    //- Limit Kd to the selected structure region domain
    Kd_.primitiveFieldRef() *= FSPair.structure().cellFields()[region];
    Kd_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volTensorField> Foam::dragModels::constant::Kd() const
{    
    return tmp<volTensorField>(new volTensorField(Kd_));
}


// ************************************************************************* //
