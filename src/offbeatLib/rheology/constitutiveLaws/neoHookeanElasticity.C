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

#include "neoHookeanElasticity.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanElasticity, 0);
    addToRunTimeSelectionTable
    (
        constitutiveLaw, 
        neoHookeanElasticity, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::neoHookeanElasticity::neoHookeanElasticity
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    constitutiveLaw(mesh, lawDict)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::neoHookeanElasticity::~neoHookeanElasticity()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::neoHookeanElasticity::correct
(
    volScalarField& sigmaHyd, 
    volSymmTensorField& sigmaDev, 
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
    const volScalarField& mu(mesh_.lookupObject<volScalarField>("mu"));
    const volScalarField& threeK(mesh_.lookupObject<volScalarField>("threeK"));
    const volTensorField& F = mesh_.lookupObject<volTensorField>("F");

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];       

        // Calculate the Jacobian of the deformation gradient
        scalar J = det(F[cellI]);

        // Calculate the volume preserving left Cauchy Green strain
        symmTensor bEbar = pow(J, -2.0/3.0)*symm(F[cellI] & F[cellI].T());

        // Calculate the deviatoric stress
        sigmaHyd[cellI] = (1.0/J)*(0.5*(threeK[cellI]/3.0)*(pow(J,2) - 1));

        // Calculate the Cauchy stress
        sigmaDev[cellI] = mu[cellI]*dev(bEbar);

        const cell& c = mesh_.cells()[cellI]; 

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);   

            if (patchID > -1 and F.boundaryField()[patchID].size())
            {
                const label faceID = 
                mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                const scalarField& muF = mu.boundaryField()[patchID];
                const scalarField& threeKF = threeK.boundaryField()[patchID];
                const tensorField& FF = F.boundaryField()[patchID];

                scalar J = det(FF[faceID]);

                symmTensor bEbar = 
                pow(J, -2.0/3.0)*symm(FF[faceID] & FF[faceID].T());

                sigmaHyd[faceID] = (1.0/J)*(0.5*(threeKF[faceID]/3.0)*(pow(J,2) - 1));
                sigmaDev[faceID] = muF[faceID]*dev(bEbar);
            }
        }
    }
}

// ************************************************************************* //


