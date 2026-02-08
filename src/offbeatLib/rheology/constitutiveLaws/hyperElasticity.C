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

#include "hyperElasticity.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hyperElasticity, 0);
    addToRunTimeSelectionTable
    (
        constitutiveLaw, 
        hyperElasticity, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hyperElasticity::hyperElasticity
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    constitutiveLaw(mesh, lawDict)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hyperElasticity::~hyperElasticity()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::hyperElasticity::correct
(
    volScalarField& sigmaHyd, 
    volSymmTensorField& sigmaDev, 
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
    const volScalarField& mu(mesh_.lookupObject<volScalarField>("mu"));
    const volScalarField& lambda(mesh_.lookupObject<volScalarField>("lambda"));
    const volTensorField& F = mesh_.lookupObject<volTensorField>("F");

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];

        // Calculate the right Cauchy–Green deformation tensor
        symmTensor C = symm(F[cellI].T() & F[cellI]);

        // Calculate the Green strain tensor
        symmTensor E = 0.5*(C - I);

        // Calculate the 2nd Piola Kirchhoff stress
        symmTensor S = 2.0*mu[cellI]*E + lambda[cellI]*tr(E)*I;

        // Calculate the Jacobian of the deformation gradient
        scalar J = det(F[cellI]);

        // Convert the 2nd Piola Kirchhoff stress to the Cauchy stress
        // sigma = (1.0/J)*symm(F() & S & F().T());
        symmTensor sigma = (1.0/J)*transform(F[cellI], S);

        // hyperElasticity laws does not split between hydrostatic and 
        // deviatoric component. 
        // NOTE: maybe sigmaHyd and sigmaDev could be obtained directly from
        // the 2nd Piola Kirchoff stress tensor
        sigmaHyd[cellI] = tr(sigma)/3.0;
        sigmaDev[cellI] = sigma - sigmaHyd[cellI]*I;

        const cell& c = mesh_.cells()[cellI]; 

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);   

            if (patchID > -1 and F.boundaryField()[patchID].size())
            {
                const label faceID = 
                mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                const scalarField& muF = mu.boundaryField()[patchID];
                const scalarField& lambdaF = lambda.boundaryField()[patchID];
                const tensorField& FF = F.boundaryField()[patchID];

                symmTensorField& sigmaDevF = 
                sigmaDev.boundaryFieldRef()[patchID];
                scalarField& sigmaHydF = 
                sigmaHyd.boundaryFieldRef()[patchID];

                symmTensor C = symm(FF[faceID].T() & FF[faceID]);
                symmTensor E = 0.5*(C - I);
                symmTensor S = 2.0*muF[faceID]*E + lambdaF[faceID]*tr(E)*I;
                scalar J = det(FF[faceID]);
                symmTensor sigma = (1.0/J)*transform(FF[faceID], S);

                sigmaHydF[faceID] = tr(sigma)/3.0;
                sigmaDevF[faceID] = sigma - sigmaHydF[faceID]*I;
            }
        }
    }
}

// ************************************************************************* //


