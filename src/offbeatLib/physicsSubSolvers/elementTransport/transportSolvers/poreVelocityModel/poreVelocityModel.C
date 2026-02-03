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

#include "poreVelocityModel.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"

#include "fvMatrix.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(poreVelocityModel, 0);
    defineRunTimeSelectionTable(poreVelocityModel, dictionary);
    addToRunTimeSelectionTable
    (
        poreVelocityModel, 
        poreVelocityModel, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::poreVelocityModel::poreVelocityModel
(
    const fvMesh& mesh, 
    const dictionary& dict,
    const word defaultModel
)
:
    mesh_(mesh),
    materialModelDict_(dict),
    poreVelocity_
    (
        createOrLookup<vector>
        (
            mesh, 
            "poreVelocity", 
            dimLength/dimTime, 
            vector::zero,
            "calculated"
        )
    ),
    alphaPoreVelocity_(dict.lookupOrDefault<scalar>("alphaPoreVelocity", 1.0)),
    regularizeGradT_(dict.lookupOrDefault<bool>("regularizeGradT", false)),
    regularizationScale_(dict.lookupOrDefault<scalar>("regularizationScale", 0)),
    regGradT_
    (
        IOobject
        (
            "regGradT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("", dimTemperature/dimLength, vector::zero),
        "zeroGradient"
    )
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::poreVelocityModel>
Foam::poreVelocityModel::New
(
    const fvMesh& mesh, 
    const dictionary& dict,
    const word defaultModel
)
{
    
    const word poreVelocityModelName =
        dict.lookupOrDefault<word>("poreVelocityModel", defaultModel);

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(poreVelocityModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("poreVelocityModel::New(const fvMesh&, const dictionary&)")
            << "Unknown poreVelocityModel type "
            << poreVelocityModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info << "Selecting poreVelocityModel type "
             << poreVelocityModelName << endl;
    }

    autoPtr<Foam::poreVelocityModel> 
    poreVelocityModelPtr(cstrIter()(mesh, dict, defaultModel));

    return poreVelocityModelPtr;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::poreVelocityModel::~poreVelocityModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volVectorField Foam::poreVelocityModel::regularizeGradT()
{

    // Initialize components 
    volScalarField GTR_x
    (
        IOobject
        (
            "GTR_x",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimTemperature/dimLength, 0.0),
        "zeroGradient"
    ); 
    volScalarField GTR_y
    (
        IOobject
        (
            "GTR_y",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimTemperature/dimLength, 0.0),
        "zeroGradient"
    ); 
    volScalarField GTR_z
    (
        IOobject
        (
            "GTR_z",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimTemperature/dimLength, 0.0),
        "zeroGradient"
    );

    // Assign components
    volVectorField GT = mesh_.lookupObject<volVectorField>("gradT");
    GTR_x = GT.component(vector::X);
    GTR_y = GT.component(vector::Y);
    GTR_z = GT.component(vector::Z);

    // Define length scale
    scalar L = regularizationScale_;

    // Create a diffusivity field that diffuses over the desired space scale
    scalar nuValue(pow(L, 2.0) / mesh_.time().deltaTValue());
    dimensionedScalar nu("", dimArea/dimTime, nuValue);

    fvScalarMatrix eq_x
    (
        fvm::ddt(GTR_x)
      - fvm::laplacian(nu, GTR_x)
    );
    fvScalarMatrix eq_y
    (
        fvm::ddt(GTR_y)
      - fvm::laplacian(nu, GTR_y)
    );
    fvScalarMatrix eq_z
    (
        fvm::ddt(GTR_z)
      - fvm::laplacian(nu, GTR_z)
    );

    eq_x.solve();
    eq_y.solve();
    eq_z.solve();

    // Re-pack components into a vector Field
    const vector i(1,0,0);
    const vector j(0,1,0);
    const vector k(0,0,1);

    regGradT_ = GTR_x * i + GTR_y * j + GTR_z * k;

    regGradT_.correctBoundaryConditions();

    return regGradT_;
}

// ************************************************************************* //
