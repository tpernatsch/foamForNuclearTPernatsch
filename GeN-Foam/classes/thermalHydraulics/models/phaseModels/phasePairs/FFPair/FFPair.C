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
    rhoContinuous_
    (
        IOobject
        (
            "rhoContinuous",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimDensity, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    nuContinuous_
    (
        IOobject
        (
            "nuContinuous",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimArea/dimTime, SMALL),
        zeroGradientFvPatchScalarField::typeName
    ),
    PrContinuous_
    (
        IOobject
        (
            "PrContinuous",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Re_
    (
        IOobject
        (
            IOobject::groupName("Re", this->name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 1),
        zeroGradientFvPatchScalarField::typeName
    ),
    minRe_
    (
        dimensionedScalar::lookupOrDefault
        (
            "residualFluidFluidRe",
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
            1e-3
        )
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
    scalarField continuity1(1.0-fluid1_.dispersion());
    scalarField continuity2(1.0-fluid2_.dispersion());

    DhDispersed_.primitiveFieldRef() = 
        Dh1.primitiveField()*fluid1_.dispersion() 
    +   Dh2.primitiveField()*fluid2_.dispersion();
    DhDispersed_.correctBoundaryConditions();

    DhContinuous_.primitiveFieldRef() = 
        Dh1.primitiveField()*continuity1
    +   Dh2.primitiveField()*continuity2;
    DhContinuous_.correctBoundaryConditions();

    rhoContinuous_.primitiveFieldRef() =
        fluid1_.thermo().rho()().primitiveField()*continuity1
    +   fluid2_.thermo().rho()().primitiveField()*continuity2;
    rhoContinuous_.correctBoundaryConditions();

    nuContinuous_.primitiveFieldRef() = 
        fluid1_.thermo().nu()()*continuity1
    +   fluid2_.thermo().nu()()*continuity2;
    nuContinuous_.correctBoundaryConditions();
    
    PrContinuous_.primitiveFieldRef() = 
        fluid1_.Pr()*continuity1
    +   fluid2_.Pr()*continuity2;
    PrContinuous_.correctBoundaryConditions();
    
    magUr_ = mag(U1-U2);
    
    Re_ = max(magUr_*DhDispersed_/nuContinuous_, minRe_);
}


// ************************************************************************* //
