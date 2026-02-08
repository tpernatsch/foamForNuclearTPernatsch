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

#include "TSSpGeneric.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"
#include "Tuple2.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TSSpGeneric, 0);
    addToRunTimeSelectionTable
    (
        TSSpModel,
        TSSpGeneric, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TSSpGeneric::TSSpGeneric
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& lawDict
)
:
    TSSpModel(mesh, mat, lawDict),
    TSSp0_(readScalar(lawDict.lookup("TSSp0"))),
    Qp_(readScalar(lawDict.lookup("Qp")))
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TSSpGeneric::~TSSpGeneric()
{}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::TSSpGeneric::correctTSSp
(
    const labelList& addr
)
{
    const volScalarField& T_(mesh_.lookupObject<volScalarField>("T"));
    const scalarField& Ti = T_.internalField();

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];

        TSSp_[cellI] = TSSp0_ * exp(-Qp_ / (8.31446 * Ti[cellI]));
        
        const cell& c = mesh_.cells()[cellI];  
        forAll(c, faceI)
        {
            const label patchID = 
            mesh_.boundaryMesh().whichPatch(c[faceI]);

            if (patchID > -1 and TSSp_.boundaryField()[patchID].size())
            {
                const label faceID = 
                mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                scalarField& TSSpP(TSSp_.boundaryFieldRef()[patchID]);
                    
                const scalarField& Tp(T_.boundaryField()[patchID]);
                TSSpP[faceID] = TSSp0_ * exp(-Qp_ / (8.31446 * Tp[faceID]));
            }
        }
    }
    TSSp_.correctBoundaryConditions();
}
// ************************************************************************* //