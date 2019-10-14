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

#include "phaseConstant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(phaseConstant, 0);
    addToRunTimeSelectionTable(dragModel, phaseConstant, FFDragModels);
    addToRunTimeSelectionTable(dragModel, phaseConstant, FSDragModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::phaseConstant::phaseConstant
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
    alpha1_(FFPair.fluid1()),
    alpha2_(FFPair.fluid2()),
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

Foam::dragModels::phaseConstant::phaseConstant
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
    alpha1_(FSPair.fluidRef()),
    alpha2_(FSPair.structure().alphaFields()[region]),
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
    //- Kd is limited to structure region because of how alpha2 is defined
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volTensorField> Foam::dragModels::phaseConstant::Kd() const
{   
    return 
        tmp<volTensorField>
        (
            new volTensorField
            (
                    Kd_*4*alpha1_*alpha2_
                /   sqr(max(alpha1_+alpha2_, 1e-6))
            )
        );
}


// ************************************************************************* //
