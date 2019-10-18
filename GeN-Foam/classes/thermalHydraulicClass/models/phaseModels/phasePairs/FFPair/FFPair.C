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

#include "FFPair.H"
#include "fluid.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FFPair, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FFPair::FFPair
(
    const fluid& fluid1,
    const fluid& fluid2
)
:
    phasePair
    (
        fluid1.mesh(),
        fluid1.name(),
        fluid2.name()
    ),
    fluid1_(fluid1),
    fluid2_(fluid2),
    DhDispersed_
    (
        IOobject
        (
            "DhDispersed",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimLength, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    DhContinuous_
    (
        IOobject
        (
            "DhContinuous",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimLength, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    nuDispersed_
    (
        IOobject
        (
            "nuDispersed",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimArea/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
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
    magUr_
    (
        IOobject
        (
            "magUr",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimVelocity, 0),
        zeroGradientFvPatchScalarField::typeName
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FFPair::correct()
{
    //- Init refs and dispersion marker fields
    const volVectorField& U1(fluid1_.U());
    const volVectorField& U2(fluid2_.U());
    const volScalarField& Dh1(fluid1_.Dh());
    const volScalarField& Dh2(fluid2_.Dh());

    DhDispersed_.primitiveFieldRef() = 
            Dh1.primitiveField()*fluid1_.dispersion() 
        +   Dh2.primitiveField()*fluid2_.dispersion();

    DhContinuous_.primitiveFieldRef() = 
            Dh1.primitiveField()*(1.0-fluid1_.dispersion())
        +   Dh2.primitiveField()*(1.0-fluid2_.dispersion());

    nuDispersed_.primitiveFieldRef() = 
            fluid1_.thermo().nu()().primitiveField()*fluid1_.dispersion() 
        +   fluid2_.thermo().nu()().primitiveField()*fluid2_.dispersion();
    
    nuDispersed_ = 
        max
        (
            nuDispersed_, 
            dimensionedScalar("", dimArea/dimTime, VSMALL)
        );

    magUr_ = mag(U1-U2);

    Re_ = magUr_*DhDispersed_/nuDispersed_;

    DhDispersed_.correctBoundaryConditions();
    DhContinuous_.correctBoundaryConditions();
    nuDispersed_.correctBoundaryConditions();
    Re_.correctBoundaryConditions();
    magUr_.correctBoundaryConditions();   

    
}

Foam::tmp<Foam::volScalarField> Foam::FFPair::alphaSum() const
{
    tmp<volScalarField> tsum
    (
        new volScalarField
        (
            fluid1_+fluid2_
        )
    );
    volScalarField& sum = tsum.ref();
    sum.correctBoundaryConditions();
    return tsum;
        
}

// ************************************************************************* //
