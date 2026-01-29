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

#include "corrosion.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "mechanicsSubSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(corrosion, 0);
    defineRunTimeSelectionTable(corrosion, dictionary);
    addToRunTimeSelectionTable
    (
        corrosion, 
        corrosion, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::corrosion::corrosion
(
    fvMesh& mesh,
    const materials& mat,
    const dictionary& corrosionDict
)
:
    mesh_(mesh),
    mat_(mat),
    corrosionDict_(corrosionDict),
    oxideThickness_
    (
        IOobject
        (
            "oxideThickness",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("oxideThickness", dimLength, 0.0)
    ),
    DOxideThickness_
    (
        IOobject
        (
            "DOxideThickness",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("DOxideThickness", dimLength, 0.0)
    ),
    DMetalThickness_
    (
        IOobject
        (
            "DMetalThickness",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("DMetalThickness", dimLength, 0.0)
    ),
    updateMesh_(corrosionDict_.lookupOrDefault<bool>("updateMesh", false))
{
    oxideThickness_.storeOldTime();
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::corrosion>
Foam::corrosion::New
(
    fvMesh& mesh,
    const materials& mat,
    const dictionary& solverDict
)
{
    dictionary corrosionOptDict
    (
        solverDict.subOrEmptyDict("corrosionOptions")
    );

    word type = solverDict.lookupOrDefault<word>("corrosion", "constant");
    
    // Map deprecated `fromLatestTime` into new name `constant`
    if (type == "fromLatestTime")
    {
        WarningIn("elementTransport::New(fvMesh&, const materials&, const dictinoary&)")
            << "Type 'fromLatestTime' for `corrosion` is deprecated. " 
            << "Please use 'constant' instead." << endl;
        
        type = "constant";
    }

    Info << "Selecting corrosion: " << type << endl;
    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("corrosion::New(const fvMesh& mesh, volScalarField& Q)")
            << "Unknown corrosion type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting corrosion type "
            << type << endl;
    }
    
    return autoPtr<corrosion>(cstrIter()
        (
            mesh, 
            mat, 
            corrosionOptDict
        ));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::corrosion::~corrosion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::corrosion::updateMesh()
{
    if(updateMesh_)
    {
        // Calculate surface displacement proportional to oxide increment
        surfaceVectorField nf(mesh_.Sf()/mesh_.magSf());
        surfaceVectorField DDf = -DMetalThickness_*nf;

        tmp<volVectorField> tmetalThicknessDD
        (
            new volVectorField
            (
                IOobject
                (
                    "tmetalThicknessDD",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector("0", dimLength, vector::zero)
            )
        );

        volVectorField& metalThicknessDD(tmetalThicknessDD.ref());

        forAll(metalThicknessDD.boundaryField(), patchID)
        {
            labelList faceCells(mesh_.boundary()[patchID].faceCells());

            // Set the boundary displacement
            metalThicknessDD.boundaryFieldRef()[patchID] = 
            DDf.boundaryField()[patchID];

            // Set the inner field displacement
            forAll(faceCells, cellI)
            {
                metalThicknessDD[faceCells[cellI]] += DDf.boundaryField()[patchID][cellI];
            }
        }

        metalThicknessDD.correctBoundaryConditions();

        // Call mechanicsSubSolver static function
        mechanicsSubSolver::updateMesh(mesh_, metalThicknessDD);
    }      
}
// ************************************************************************* //
