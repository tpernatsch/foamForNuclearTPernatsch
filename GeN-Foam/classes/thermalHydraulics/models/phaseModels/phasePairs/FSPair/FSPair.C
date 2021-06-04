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
    minRe_
    (
        dimensionedScalar::lookupOrDefault
        (
            "residualFluidStructureRe",
            IOdictionary
            (
                IOobject
                (
                    "phaseProperties",
                    mesh_.time().constant(),
                    mesh_
                )
            ),
            dimless,
            10
        )
    ),
    Re_
    (
        IOobject
        (
            (
                (fluid.name() != "") ?
                IOobject::groupName("Re", this->name()) :
                "Re"
            ),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        minRe_,
        zeroGradientFvPatchScalarField::typeName
    ),
    lRe_
    (
        IOobject
        (
            (
                (fluid.name() != "") ?
                IOobject::groupName("lRe", this->name()) :
                "lRe"
            ),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("", dimless, vector::one),
        zeroGradientFvPatchScalarField::typeName
    )
{
    //- If a fluid-structure pair is being instantiated in a monoPhase solver
    //  (where the fluid has no name), rename the pair simply as "FSPair". This
    //  is required to make things work in the makeInTable functions in 
    //  dragModel, heatTransfer (called when the models are constructed in each
    //  regime class, see regime.C). In particular, this matters when reading
    //  structure region keywords (i.e. cellZones in the byZone structureModel)
    //  in the phaseProperties dict to define which drag/heatTransfer model 
    //  should be used in which region 
    if (fluid.name() == "")
    {
        this->rename("FSPair");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FSPair::correct()
{
    //- Init refs
    const volVectorField& U(fluid_.U());

    //- Previously fluid_.Dh(). The difference is that fluid_.Dh() will be
    //  different than the structure Dh in regions where fluid_ is dispersed.
    //  Nonetheless, the Reynolds for the fluid-structure pair when the fluid
    //  is dispersed is fairly meaningless, so I changed it to utilize the
    //  structure_.Dh() which does not change in time.
    const volScalarField& Dh(structure_.Dh());
    const volVectorField& lDh(structure_.lDh());
    
    //- Correct
    Re_ = fluid_.normalized()*mag(U)*Dh/fluid_.thermo().nu()();
    Re_ = max(Re_, minRe_);

    //- Rotate U in the local frame, make an outer product with lDh so that
    //  we get a matrix whose diagonal elements are the components of lRe.
    //  Then directly set these diagonal elemets (I guess this is the most
    //  efficient approch)
    tmp<volTensorField> tlRet = 
        fluid_.normalized()*(structure_.Rg2l()&U)*lDh/fluid_.thermo().nu()();
    lRe_.replace(0, max(tlRet().component(0), minRe_));
    lRe_.replace(1, max(tlRet().component(4), minRe_));
    lRe_.replace(2, max(tlRet().component(8), minRe_));

    Re_.correctBoundaryConditions();
    lRe_.correctBoundaryConditions();

    if (fPtr_.valid())
    {
        volScalarField& f(fPtr_());
        f.correctBoundaryConditions();
    }
}


// ************************************************************************* //
