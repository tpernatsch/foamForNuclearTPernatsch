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

#include "elasticity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(elasticity, 0);
    addToRunTimeSelectionTable
    (
        constitutiveLaw, 
        elasticity, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::elasticity::elasticity
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    constitutiveLaw(mesh, lawDict)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::elasticity::~elasticity()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::elasticity::correct
(
    volScalarField& sigmaHyd, 
    volSymmTensorField& sigmaDev, 
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
    const volScalarField& mu(mesh_.lookupObject<volScalarField>("mu"));
    const volScalarField& threeK(mesh_.lookupObject<volScalarField>("threeK"));

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        
        scalar trEpsilon = tr(epsilonEl[cellI]);
        scalar K = 1.0/3.0*threeK[cellI];
        
        sigmaHyd[cellI] = K*trEpsilon;
        sigmaDev[cellI] = 2*mu[cellI]*dev(epsilonEl[cellI]);

        const cell& c = mesh_.cells()[cellI];         

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);   

            if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
            {
                const label faceID = 
                mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                const scalarField& muF = mu.boundaryField()[patchID];
                const scalarField& threeKF = threeK.boundaryField()[patchID];

                const symmTensorField& epsilonElF = 
                epsilonEl.boundaryField()[patchID];
                
                scalarField& sigmaHydF = 
                sigmaHyd.boundaryFieldRef()[patchID];

                symmTensorField& sigmaDevF = 
                sigmaDev.boundaryFieldRef()[patchID];

                scalar trEpsilonF = tr(epsilonElF[faceID]);
                scalar KF = 1.0/3.0*threeKF[faceID];
                
                sigmaHydF[faceID] = KF*trEpsilonF;
                sigmaDevF[faceID] = 2*muF[faceID]*dev(epsilonElF[faceID]); 
            }
        }
    }
}

// ************************************************************************* //


