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

#include "UO2meltingMATPRO.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(UO2meltingMATPRO, 0);
    addToRunTimeSelectionTable(failureModel, UO2meltingMATPRO, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UO2meltingMATPRO::UO2meltingMATPRO
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    failureModel(mesh, dict),
    Tmelt_(failureModelDict_.lookupOrDefault<scalar>("Tmelt", 3113.0))
{    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::UO2meltingMATPRO::~UO2meltingMATPRO()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::UO2meltingMATPRO::isFailed
(
    const labelList& addr
)
{

    // Set failure switch to false
    failed_ = false;

    // Set counter of failed cells to 0
    scalar counter(0);

    // Constant reference to the fields
    const volScalarField& T(mesh_.lookupObject<volScalarField>("T"));

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];

        if ( T[cellI] >= Tmelt_)
        {
            if( !mesh_.foundObject<volScalarField>("failedMaterial"))
            {
                //- Initialize failedMaterial field (from parent class "failureModel") 
                initializeFailedMaterialField();
            }

            failedMaterial_()[cellI] = 1;
            failed_ = true;
            counter++;
        }
    }

    return failed_;

}


// ************************************************************************* //
