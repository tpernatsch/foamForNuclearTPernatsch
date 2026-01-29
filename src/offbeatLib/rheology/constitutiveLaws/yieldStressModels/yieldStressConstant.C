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

#include "yieldStressConstant.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(yieldStressConstant, 0);
    addToRunTimeSelectionTable
    (
        yieldStressModel, 
        yieldStressConstant, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::yieldStressConstant::yieldStressConstant
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    yieldStressModel(mesh, lawDict),
    sigmaYValue_(readScalar(lawDict.lookup("sigmaY")))
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::yieldStressConstant::~yieldStressConstant()
{}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::yieldStressConstant::correctYieldStress
(
    const labelList& addr
)
{
    forAll(addr, i)
    {  
        //- cells
        sigmaY_[addr[i]] = sigmaYValue_;

        //- boundary fields
        const cell& c = mesh_.cells()[addr[i]];  
        forAll(c, faceI)
        {
            const label patchID = 
            mesh_.boundaryMesh().whichPatch(c[faceI]);

            if (patchID > -1 and sigmaY_.boundaryField()[patchID].size())
            {
                
                const label faceID = 
                mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                scalarField& sigmaYP(sigmaY_.boundaryFieldRef()[patchID]);

                sigmaYP[faceID] = sigmaYValue_;
            }
        }
    }
}


// ************************************************************************* //


