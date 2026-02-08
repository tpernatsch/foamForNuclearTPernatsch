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

#include "ZircaloyOverstressBison.H"
#include "addToRunTimeSelectionTable.H"
#include "offbeatTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ZircaloyOverstressBison, 0);
    addToRunTimeSelectionTable(failureModel, ZircaloyOverstressBison, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ZircaloyOverstressBison::ZircaloyOverstressBison
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
    betaFrac_(nullptr)
{    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ZircaloyOverstressBison::~ZircaloyOverstressBison()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::ZircaloyOverstressBison::isFailed
(
    const labelList& addr
)
{

    if( betaFrac_ == nullptr )
    {
        if( !mesh_.foundObject<volScalarField>("betaFraction") )
        {
            FatalErrorInFunction()
            << "The selected failure criterion needs beta transition model to be activated.\n"
            << exit(FatalError); 
        }

        betaFrac_ = &mesh_.lookupObject<volScalarField>("betaFraction");
    }

    if( !mesh_.foundObject<volSymmTensorField>("sigmaCyl") )
    {
        FatalErrorInFunction()
        << "The selected failure criterion needs the field sigmaCyl to be activated.\n"
        << exit(FatalError); 
    }

    // Constant reference to the sigmaCyl field
    const volSymmTensorField& sigmaCyl(mesh_.lookupObject<volSymmTensorField>("sigmaCyl"));
    
    // Reference to user run time
    const offbeatTime& userRunTime(refCast<const offbeatTime>(mesh_.time()));

    // Set failure switch to false
    failed_ = false;

    forAll( patchNames_, i )
    {
        const label patchID = mesh_.boundaryMesh().findPatchID(patchNames_[i]);
        
        if ( patchID == -1 )
        {
            FatalErrorIn("Foam::ZircaloyPlasticInstabilityBison::isFailed(const labelList& addr)")
            << "Patch name '" << patchNames_[i] << "' does not exist." 
            << exit(FatalError); 
        }

        // Take hoop component of the stress
        const scalarField& sigmaHoop(sigmaCyl.component(symmTensor::YY)->boundaryField()[patchID]);

        // Constant reference to the T field
        const scalarField& T(mesh_.lookupObject<volScalarField>("T").boundaryField()[patchID]);

        // Constant reference to the beta field
        const scalarField& betaFrac(betaFrac_->boundaryField()[patchID]);

        // Initialize parameter for the correlation
        scalar a;
        scalar b;

        forAll(T, j)
        {
            const scalar bf(betaFrac[j]);

            if ( bf == 0)
            {
                a = 830e6;
                b = 1e-3;
            }
            else if ( bf < 0.5)
            {
                scalar aAlpha(830e6);
                scalar a5050(3000e6);
                scalar bAlpha(1e-3);
                scalar b5050(3e-3);

                a = exp( bf/0.5*(log(a5050)-log(aAlpha))+log(aAlpha) );
                b = bf/0.5*(b5050-bAlpha) + bAlpha;

            }
            else
            {
                scalar a5050(3000e6);
                scalar aBeta(2300e6);

                a = exp( (bf-0.5)/0.5*(log(aBeta)-log(a5050))+log(a5050) );
                b = 3e-3;
            }

            // Calculate the burst stress
            // TODO: OXIDATION contribution is not considered here yet !!!
            const scalar burstStress( a*exp(-b*T[j]) );

            if ( sigmaHoop[j] >= burstStress )
            {
                if( !mesh_.foundObject<volScalarField>("failedMaterial"))
                {
                    // Initialize failedMaterial field (from parent class "failureModel") 
                    initializeFailedMaterialField();
                }

                failed_ = true;

                // Print out information about failure
                Info << "\nFailed at time " <<  userRunTime.userTime() << userRunTime.unit() 
                     << " with local hoop stress = " << sigmaHoop[j]
                     << " on patch " << patchNames_[i]
                     << endl;

                // Set the field failedMaterial_ to 1 
                failedMaterial_().boundaryFieldRef()[patchID][j] = 1;
            }

        }

    }

    return failed_;
}


// ************************************************************************* //
