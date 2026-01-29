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

#include "ZircaloyPlasticInstabilityBison.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"
#include "offbeatTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ZircaloyPlasticInstabilityBison, 0);
    addToRunTimeSelectionTable(failureModel, ZircaloyPlasticInstabilityBison, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ZircaloyPlasticInstabilityBison::ZircaloyPlasticInstabilityBison
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict
)
:
    failureModel(mesh, dict, baseDict),
#ifdef OPENFOAMFOUNDATION
    patchNames_(failureModelDict_.lookup<wordList>("patchNames")),
#elif OPENFOAMESI
    patchNames_(failureModelDict_.get<wordList>("patchNames")),
#endif
    burstStrainRate_
    (
        failureModelDict_.lookupOrDefault<scalar>("burstStrainRate", 2.78e-2)
    )
{    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ZircaloyPlasticInstabilityBison::~ZircaloyPlasticInstabilityBison()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::ZircaloyPlasticInstabilityBison::isFailed
(
    const labelList& addr
)
{
    if
    ( 
        !mesh_.foundObject<volScalarField>("DepsilonCreepEq")
        or 
        !mesh_.foundObject<volScalarField>("DEpsilonPEq")
    )
    {
        FatalErrorIn("Foam::ZircaloyPlasticInstabilityBison::isFailed(const labelList& addr)")
        << "The selected failure criterion needs the plasticCreep model to be activated.\n"
        << exit(FatalError); 
    }

    // Reference to creep deformation increment
    const volSymmTensorField& DepsilonCreep = 
        mesh_.lookupObject<volSymmTensorField>("DepsilonCreep");
    
    // Reference to plastic deformation increment
    const volSymmTensorField& DepsilonP = 
        mesh_.lookupObject<volSymmTensorField>("DEpsilonP");

    // Compute plastic deformation rate
    const volSymmTensorField plasticDefRate = 
        (DepsilonCreep+DepsilonP)/mesh_.time().deltaTValue();

    // Compute effective plastic strain rate
    const volScalarField effPlasticStrainRate = sqrt(2.0/3.0*magSqr(plasticDefRate));

    // User run time 
    const offbeatTime& userRunTime(refCast<const offbeatTime>(mesh_.time()));
    
    // Set failure switch to false
    failed_ = false;
    
    forAll(patchNames_, i)
    {
        // Get name of the patch 
        const label patchID = mesh_.boundaryMesh().findPatchID(patchNames_[i]);

        if ( patchID == -1 )
        {
            FatalErrorIn("Foam::ZircaloyPlasticInstabilityBison::isFailed(const labelList& addr)")
            << "Patch name '" << patchNames_[i] << "' does not exist." 
            << exit(FatalError); 
        }

        // Get maximum value and its position
        forAll(effPlasticStrainRate.boundaryField()[patchID], j)
        {
            const scalar eRateP =  effPlasticStrainRate.boundaryField()[patchID][j];

            if ( eRateP >= burstStrainRate_ )
            {
                if( !mesh_.foundObject<volScalarField>("failedMaterial"))
                {
                    // Initialize failedMaterial_ (from parent class "failureModel") 
                    initializeFailedMaterialField();
                }

                failed_ = true;

                // Print out information about failure
                Info << "\nFailed at time " <<  userRunTime.userTime() << userRunTime.unit() 
                     << " with a plastic strain rate of " 
                     << eRateP << " 1/s on patch " <<  patchNames_[i]
                     << endl;

                // Set the field failedMaterial_ to 1 
                failedMaterial_().boundaryFieldRef()[patchID][j] = 1;

            }
        }
    }

    return failed_;}


// ************************************************************************* //
