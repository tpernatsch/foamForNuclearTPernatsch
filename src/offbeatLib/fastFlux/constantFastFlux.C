/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "constantFastFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"
#include "IndirectList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantFastFlux, 0);
    
    addToRunTimeSelectionTable
    (
        fastFlux, 
        constantFastFlux, 
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::constantFastFlux::advanceFluence()
{
    scalar Dt(mesh_.time().deltaT().value());

    scalarField& fastFluencei = fastFluence_.ref();  

    const scalarField& fastFluenceOldi = 
    fastFluence_.oldTime().internalField();
    
    const scalarField& fastFluxi = flux_.internalField();

    fastFluencei = fastFluenceOldi + fastFluxi*Dt;

    fastFluence_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantFastFlux::constantFastFlux
(
    const fvMesh& mesh,
    const materials& materials,
    const dictionary& fastFluxDict
)
:
    fastFlux(mesh, materials, fastFluxDict),
    fastFluence_
    (
        IOobject
        (
            "fastFluence",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("fastFluence", dimless, 0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    flux_
    (
        IOobject
        (
            "fastFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("fastFlux", dimless, 0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    currentTime_(-1)
{
    // Prepare old time field
    flux_.oldTime();
    fastFluence_.oldTime();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantFastFlux::~constantFastFlux()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::constantFastFlux::correct()
{
    if(!mesh_.foundObject<volScalarField>("fastFluxPrevIter"))
    {
        flux_.storePrevIter();
    }

    flux_.prevIter().storePrevIter();        
    flux_.storePrevIter();

    const scalarField& fastFluxPrev(flux_.prevIter());
    const scalarField& fastFluxPrevPrev(flux_.prevIter().prevIter());

    // Correct only once per time step or until flux changes
    if
    ( 
        gMax(mag(fastFluxPrevPrev - fastFluxPrev)) > VSMALL 
        or
        mesh_.time().value() > currentTime_
    )
    { 

        currentTime_ = mesh_.time().value();
        advanceFluence();
    }

}

// ************************************************************************* //
