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

#include "elementTransport.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(elementTransport, 0);
    defineRunTimeSelectionTable(elementTransport, dictionary);

    addToRunTimeSelectionTable
    (
        elementTransport, 
        elementTransport, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::elementTransport::elementTransport
(
    fvMesh& mesh,
    const materials& mat,
    const dictionary& elementTransportDict
)
:
    mesh_(mesh),
    mat_(mat),
    elementTransportDict_(elementTransportDict)
{
}


Foam::autoPtr<Foam::elementTransport>
Foam::elementTransport::New
(
    fvMesh& mesh,
    const materials& mat,
    const dictionary& solverDict
)
{
    word type;

    dictionary elementTransportDict
    (
        solverDict.subOrEmptyDict("elementTransportOptions")
    );

    solverDict.lookup("elementTransport") >> type;
    Info << "Selecting elementTransportSolver " << type << endl;
    
    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("elementTransport::New(const fvMesh&, const materials&, const dictionary&)")
            << "Unknown elementTransport type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<elementTransport>
    (
        cstrIter()
        (
            mesh, mat, elementTransportDict
        )
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::elementTransport::~elementTransport()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
