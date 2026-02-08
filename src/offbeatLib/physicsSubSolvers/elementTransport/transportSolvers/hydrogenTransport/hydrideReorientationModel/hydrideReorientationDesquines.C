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

#include "hydrideReorientationDesquines.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"
#include "zeroGradientFvPatchFields.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * Static Data Member  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hydrideReorientationDesquines, 0);
    addToRunTimeSelectionTable
    (
        hydrideReorientationModel, 
        hydrideReorientationDesquines, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hydrideReorientationDesquines::hydrideReorientationDesquines
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    hydrideReorientationModel(mesh, lawDict),
    sigma0_
    (
        IOobject
        (
            "sigma0",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPressure, 0.0),
        calculatedFvPatchField<scalar>::typeName
    ),
    sigma100_
    (
        IOobject
        (
            "sigma100",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPressure, 0.0),
        calculatedFvPatchField<scalar>::typeName
    ),
    fractionalHydridesChange_
    (
        IOobject
        (
            "fractionalHydridesChange",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    T_tssd_
    (
        IOobject
        (
            "T_tssd",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0.0),
        calculatedFvPatchField<scalar>::typeName
    )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hydrideReorientationDesquines::~hydrideReorientationDesquines()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::hydrideReorientationDesquines::correctReorientation
(
    const labelList& addr
)
{
    // Using Kammenzind's coefficients
    const scalar TSSd_A = 66000;
    const scalar TSSd_Q = 32144;
    //const scalar mu_H0 = -656;

    //Using Une's coefficients
    //const scalar TSSd_A = 7.8029e4;
    //const scalar TSSd_Q = 32620;
    
    const volScalarField& T_(mesh_.lookupObject<volScalarField>("T"));
    const volScalarField& Ctot_(mesh_.lookupObject<volScalarField>("Ctot"));
    const volScalarField& Cpp_(mesh_.lookupObject<volScalarField>("Cpp"));
    // hydrideChange_ is Cpp_ - Cpp_.oldTime()
    const volScalarField& hydrideChange_(mesh_.lookupObject<volScalarField>("hydrideChange"));
    const volSymmTensorField& sigmaCyl_(mesh_.lookupObject<volSymmTensorField>("sigmaCyl"));
    // Take hoop component of the stress
    volScalarField sigmaHoop_(sigmaCyl_.component(symmTensor::YY));

    const volScalarField& OldRadialHydrides = radialHydrides_.oldTime();
    const volScalarField& OldFracRHy =  fracRHy_.oldTime();

    forAll(addr, addrI)
    {    
        const label cellI = addr[addrI];

        fractionalHydridesChange_[cellI] = 0;
        // reorientation only if tensile stress
        if (sigmaHoop_[cellI] < 0)
        {
            fractionalHydridesChange_[cellI] = 0;
        }

        else if (sigmaHoop_[cellI] > 0 && hydrideChange_[cellI] > 0)
        {
            //- Stress value for 100% radial hydrides (MPa)
            sigma100_[cellI] = 65*(1- exp(-Ctot_[cellI]/65)) + 110;

            //- Temperature required to dissolve Ctot during heatup (°C)
            T_tssd_[cellI] = -(TSSd_Q / 8.31446 / log(Ctot_[cellI]/TSSd_A)) - 273.15;

            //T_tssd_[cellI] = (-TSSd_Q + mu_H0) / 8.31446 / log(Ctot_[cellI]/TSSd_A) - 273.15;

            //- Threshold stress value for radial reorientation (MPa)
            sigma0_[cellI] = 0.02*Ctot_[cellI] - 0.3862*min(T_[cellI]-273.15, T_tssd_[cellI]) + 186.9;

            fractionalHydridesChange_[cellI] =  (sigmaHoop_[cellI]/1e6  - sigma0_[cellI]) / (sigma100_[cellI] - sigma0_[cellI]) * hydrideChange_[cellI];
            if(sigmaHoop_[cellI]/1e6 < sigma0_[cellI])
            {
                fractionalHydridesChange_[cellI] = 0;
            }
            else if(sigmaHoop_[cellI]/1e6 > sigma100_[cellI])
            {
                fractionalHydridesChange_[cellI] = hydrideChange_[cellI];
            }
        }

        else if (hydrideChange_[cellI] < 0)
        {
            fractionalHydridesChange_[cellI] = hydrideChange_[cellI]*OldFracRHy[cellI];
            //radialHydrides_[cellI] = OldRadialHydrides[cellI] - (hydrideChange_[cellI] * OldfracRHy_[cellI]);
        }
        
        radialHydrides_[cellI] = OldRadialHydrides[cellI] + fractionalHydridesChange_[cellI];
        if (radialHydrides_[cellI] <= 0)
        {
            radialHydrides_[cellI] = 0;
        }

        if (Cpp_[cellI] <= 0)
        {
            fracRHy_[cellI] = 0;
        } 
        else
        {
            fracRHy_[cellI] = radialHydrides_[cellI]/Cpp_[cellI];
        }

        if (fracRHy_[cellI] > 1)
        {
            fracRHy_[cellI] = 1;
        }
        else if (fracRHy_[cellI] <= 0)
        {
            fracRHy_[cellI] = 0;
        }
    }
    fracRHy_.correctBoundaryConditions();
    radialHydrides_.correctBoundaryConditions();
}

// ************************************************************************ //
