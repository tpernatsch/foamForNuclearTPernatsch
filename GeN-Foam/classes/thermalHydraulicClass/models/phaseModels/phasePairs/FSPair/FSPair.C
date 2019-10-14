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
#include "fluid.H"
#include "structureModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FSPair, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSPair::FSPair
(
    const fluid& fluid,
    const structureModel& structure
)
:
    phasePair
    (
        fluid.mesh(),
        fluid.name(),
        structure.name()
    ),
    fluid_(fluid),
    structure_(structure),
    Re_
    (
        IOobject
        (
            IOobject::groupName("Re", this->name_),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 1),
        zeroGradientFvPatchScalarField::typeName
    ),
    lRe_
    (
        IOobject
        (
            IOobject::groupName("lRe", this->name_),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("", dimless, vector::one),
        zeroGradientFvPatchScalarField::typeName
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FSPair::correct()
{
    //- Init refs and dispersion marker fields
    const volVectorField& U(fluid_.U());

    //- Previously fluid_.Dh(). The difference is that fluid_.Dh() will be
    //  different than the structure Dh in regions where fluid_ is dispersed.
    //  Nonetheless, the Reynolds for the fluid-structure pair when the fluid
    //  is dispersed is fairly meaningless, so I changed it to utilize the
    //  structure_.Dh() which does not change in time.
    const volScalarField& Dh(structure_.Dh());
    const volVectorField& lDh(structure_.lDh());

    //- Correct
    Re_ = mag(U)*Dh/fluid_.thermo().nu()();

    //- Rotate U in the local frame, make an outer product with lDh so that
    //  we can a matrix whose diagonal elements are the components of lRe.
    //  Then directly set these diagonal elemets (I guess this is the most
    //  efficient approch)
    tmp<volTensorField> tlRet = 
        (structure_.Rg2l()&U)*lDh/fluid_.thermo().nu()();
    lRe_.component(0) = tlRet().component(0);
    lRe_.component(1) = tlRet().component(4);
    lRe_.component(2) = tlRet().component(8);

    Re_.correctBoundaryConditions();
    lRe_.correctBoundaryConditions();
}

const Foam::wordList& Foam::FSPair::structureRegions() const
{
    return structure_.regions();
}


// ************************************************************************* //
