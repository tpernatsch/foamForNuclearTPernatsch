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

#include "hyperElasticMisesPlasticCreep.H"
#include "addToRunTimeSelectionTable.H"

#include "fvm.H"
#include "fvc.H"
#include "calculatedFvPatchFields.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hyperElasticMisesPlasticCreep, 0);
    addToRunTimeSelectionTable
    (
        misesPlasticity, 
        hyperElasticMisesPlasticCreep, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hyperElasticMisesPlasticCreep::hyperElasticMisesPlasticCreep
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    misesPlasticity(mesh, lawDict),
    creepModel_(creepModel::New(mesh, lawDict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hyperElasticMisesPlasticCreep::~hyperElasticMisesPlasticCreep()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::hyperElasticMisesPlasticCreep::correct
(
    volScalarField& sigmaHyd, 
    volSymmTensorField& sigmaDev, 
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
    // Calculate new epsilon creep increment and update elastic strain
    creepModel_->correctCreep(epsilonEl, addr);

    // Correct deviatoric and hydrostatic stresses according to plastic law
    misesPlasticity::correct(sigmaHyd, sigmaDev, epsilonEl, addr);

    // Scale and rotate stress --> transform S into sigma
    const volTensorField& F = mesh_.lookupObject<volTensorField>("F");
    
    forAll(addr, addrI)
    {
        // Cell adrress. Act only on element of addr
        const label cellI = addr[addrI];

        // Convert the 2nd Piola Kirchhoff stress to the Cauchy stress
        // sigma = (1.0/J)*symm(F() & S & F().T());
        symmTensor S = sigmaHyd[cellI]*I + sigmaDev[cellI];

        // Calculate the Jacobian of the deformation gradient
        scalar J = det(F[cellI]);

        symmTensor sigmaReal = (1.0/J)*transform(F[cellI], S);

        sigmaHyd[cellI] = tr(sigmaReal)/3.0;
        sigmaDev[cellI] = sigmaReal - sigmaHyd[cellI]*I;

        const cell& c = mesh_.cells()[cellI];  

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

            if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
            {
                const label faceID = 
                mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);          

                // Take references to the boundary patch fields for efficiency
                symmTensorField& sigmaDevP = 
                sigmaDev.boundaryFieldRef()[patchID];

                scalarField& sigmaHydP = 
                sigmaHyd.boundaryFieldRef()[patchID];

                // Convert the 2nd Piola Kirchhoff stress to the Cauchy stress
                // sigma = (1.0/J)*symm(F() & S & F().T());
                const tensorField& Fp = 
                mesh_.lookupObject<volTensorField>("F").boundaryField()[patchID];

                symmTensor S = sigmaHydP[faceID]*I + sigmaDevP[faceID];

                // Calculate the Jacobian of the deformation gradient
                scalar J = det(Fp[faceID]);
                
                symmTensor sigmaReal = (1.0/J)*transform(Fp[faceID], S);

                sigmaHydP[faceID] = tr(sigmaReal)/3.0;
                sigmaDevP[faceID] = sigmaReal - sigmaHydP[faceID]*I;
            }
        }
    }
}


void Foam::hyperElasticMisesPlasticCreep::correctEpsilonEl
(
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
    misesPlasticity::correctEpsilonEl(epsilonEl, addr);
};


void Foam::hyperElasticMisesPlasticCreep::correctAxialStrain
(
    volSymmTensorField& epsilon,
    const labelList& addr
)
{
    creepModel_->correctAxialStrain(epsilon, addr);
}


void Foam::hyperElasticMisesPlasticCreep::correctAdditionalStrain
(
    volSymmTensorField& additionalStrain,
    const labelList& addr
)
{
    forAll(addr, addrI)
    {
        const label cellI(addr[addrI]);

        additionalStrain[cellI] = creepModel_->additionalStrain()[cellI]
                                + epsilonP_.internalField()[cellI];
    }
}


Foam::scalar Foam::hyperElasticMisesPlasticCreep::nextDeltaT
(
    const labelList& addr
)
{
    return creepModel_->nextDeltaT(addr);
}
// ************************************************************************* //


