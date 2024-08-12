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

#include "sliceMapper.H"
#include "fvm.H"
#include "fvc.H"
#include "calculatedFvPatchFields.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sliceMapper, 0);
    defineRunTimeSelectionTable(sliceMapper, dictionary);

    addToRunTimeSelectionTable
    (
        sliceMapper, 
        sliceMapper, 
        dictionary
    );
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::sliceMapper::createSliceID()
{
    // Create sliceID field if not already available
    if (!sliceID_.valid())
    {
#ifdef OPENFOAMFOUNDATION    
        sliceID_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "sliceID",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("sliceID", dimless, 0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
#elif OPENFOAMESI
        sliceID_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "sliceID",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("sliceID", dimless, 0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
#endif        
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sliceMapper::sliceMapper
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& mapperOptDict
)
:
    regIOobject
    (
        IOobject
        (
            "sliceMapper",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    mat_(mat),
    mapperOptDict_(mapperOptDict),
    sliceAddrList_(),
    sliceID_(),
    isFuel_(),
    nSlices_(0)
{}


Foam::autoPtr<Foam::sliceMapper>
Foam::sliceMapper::New
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& solverDict
)
{
    
    // Initialize default type for thermalSolver class
    word type("byMaterial");

    dictionary mapperOptDict
    (
        solverDict.subOrEmptyDict("mapperOptions")
    );

    solverDict.lookup("sliceMapper") >> type;
    Info << "Selecting sliceMapper model " << type << endl;

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("sliceMapper::New(const fvMesh& mesh)")
            << "Unknown sliceMapper dependence type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting sliceMapper type "
            << type << endl;
    }
    
    return autoPtr<sliceMapper>(cstrIter()
        (
            mesh, 
            mat,
            mapperOptDict
        ));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sliceMapper::~sliceMapper()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// ************************************************************************* //
