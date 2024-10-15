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

#include "ZyOverstrainBISON.H"
#include "addToRunTimeSelectionTable.H"
#include "offbeatTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ZyOverstrainBISON, 0);
    addToRunTimeSelectionTable(failureModel, ZyOverstrainBISON, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ZyOverstrainBISON::ZyOverstrainBISON
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    failureModel(mesh, dict),
#ifdef OPENFOAMFOUNDATION
    patchNames_(failureModelDict_.lookup<wordList>("patchNames")),
#elif OPENFOAMESI
    patchNames_(failureModelDict_.get<wordList>("patchNames")),
#endif
    hoopStrainLimit_
    (
        failureModelDict_.lookupOrDefault<scalar>("hoopStrainLimit", 0.4)
    )
{    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ZyOverstrainBISON::~ZyOverstrainBISON()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::ZyOverstrainBISON::isFailed
(
    const labelList& addr
)
{
    if( !mesh_.foundObject<volSymmTensorField>("epsilonCreep") )
    {
        FatalErrorInFunction()
        << "The selected failure criterion needs the creep model to be activated.\n"
        << exit(FatalError); 
    }

    // Constant reference to the epsilonCreep field
    const volSymmTensorField& epsilonCreep(mesh_.lookupObject<volSymmTensorField>("epsilonCreep"));

    // User run time
    const offbeatTime& userRunTime(refCast<const offbeatTime>(mesh_.time()));

    // Set failure switch to false
    failed_ = false;

    forAll( patchNames_, i )
    {
        const label patchID = mesh_.boundaryMesh().findPatchID(patchNames_[i]);
        
        if ( patchID == -1 )
        {
            FatalErrorIn("Foam::ZyPlasticInstabilityBISON::isFailed(const labelList& addr)")
            << "Patch name '" << patchNames_[i] << "' does not exist." 
            << exit(FatalError); 
        }

        // Take reference to eps creep patch field
        const symmTensorField& epsCreepP(epsilonCreep.boundaryField()[patchID]);
        
        // Initialize the epsilonCreep tensor in cylindrical coordinates
        symmTensorField epsCreepCylP(epsCreepP.size(), symmTensor::zero);

        forAll(epsCreepCylP, j)
        {
            // Transform in cyl coordinates
            tensor R = tensor::zero;
            const vector& Ci = mesh_.boundaryMesh()[patchID].faceCentres()[j];
            const scalar  radius = sqrt(Ci[0]*Ci[0] + Ci[1]*Ci[1]);
            R.xx() =  Ci[0]/radius;
            R.xy() =  Ci[1]/radius;
            R.yy() =  Ci[0]/radius;
            R.yx() = -Ci[1]/radius;
            R.zz() = 1;
            epsCreepCylP[j] = symm( R & epsCreepP[j] & R.T() );

            if ( epsCreepCylP[j].component(symmTensor::YY) >= hoopStrainLimit_ )
            {
                if( !mesh_.foundObject<volScalarField>("failedMaterial"))
                {
                    // Initialize failedMaterial field (from parent class "failureModel") 
                    initializeFailedMaterialField();
                }

                failed_ = true;

                // Print out information about failure
                Info << "Failed at time " <<  userRunTime.userTime() << userRunTime.unit() 
                     << ", with creep hoop strain = " << epsCreepCylP[j].component(symmTensor::YY)
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
