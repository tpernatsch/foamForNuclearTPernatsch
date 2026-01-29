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

#include "TSSdMerlino.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"
#include "Tuple2.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TSSdMerlino, 0);
    addToRunTimeSelectionTable
    (
        TSSdModel, 
        TSSdMerlino, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TSSdMerlino::TSSdMerlino
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& lawDict
)
:
    TSSdModel(mesh, mat, lawDict)
{

    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TSSdMerlino::~TSSdMerlino()
{}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::TSSdMerlino::correctTSSd
(
    const labelList& addr
)
{
    const volScalarField& T_(mesh_.lookupObject<volScalarField>("T"));
    const scalarField& Ti = T_.internalField();

    const scalar TSSd0 = 106000;
    const scalar Qd = 36000;

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];

        TSSd_[cellI] = TSSd0 * exp(-Qd / (8.31446 * Ti[cellI]));

        const cell& c = mesh_.cells()[cellI];  
        forAll(c, faceI)
        {
            const label patchID = 
            mesh_.boundaryMesh().whichPatch(c[faceI]);

            if (patchID > -1 and TSSd_.boundaryField()[patchID].size())
            {
                const label faceID = 
                mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                scalarField& TSSdP(TSSd_.boundaryFieldRef()[patchID]);
                    
                const scalarField& Tp(T_.boundaryField()[patchID]);
                TSSdP[faceID] = TSSd0 * exp(-Qd / (8.31446 * Tp[faceID]));
            }
        }
    }
    TSSd_.correctBoundaryConditions();
}
// ************************************************************************* //