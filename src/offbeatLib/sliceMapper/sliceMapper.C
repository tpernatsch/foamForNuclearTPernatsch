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
    const dictionary& sliceMapperDict
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
    sliceMapperDict_(sliceMapperDict),
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
    // Default pick type
    word type("none");

    // Prepare options dict
    dictionary sliceMapperDict
    (
        solverDict.subOrEmptyDict("sliceMapperOptions")
    );

    if (solverDict.found("sliceMapper"))
    {
        if (solverDict.isDict("sliceMapper"))
        {
            // Case 1: sliceMapper is a dictionary
            sliceMapperDict = solverDict.subDict("sliceMapper");
            sliceMapperDict.lookup("type") >> type;
        }
        else
        {
            // Case 2: sliceMapper is a word (type)
            solverDict.lookup("sliceMapper") >> type;
            
            if(solverDict.found("mapperOptions"))
            {
                sliceMapperDict = solverDict.subDict("mapperOptions");
                    
                WarningIn("sliceMapper::New(const fvMesh&, materials&, const dictinoary&)")
                    << "Options dictionary 'mapperOptions' for `sliceMapper` class is deprecated. " 
                    << "Please use 'sliceMapperOptions' instead or even better create a 'sliceMapper' dict where the type of sliceMapper model is indicated by the keyword 'type'." << nl << endl;
            }
            // else, sliceMapperDict stays as sliceMapperOptions
        }
    }
    
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
            sliceMapperDict
        ));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sliceMapper::~sliceMapper()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// ************************************************************************* //
