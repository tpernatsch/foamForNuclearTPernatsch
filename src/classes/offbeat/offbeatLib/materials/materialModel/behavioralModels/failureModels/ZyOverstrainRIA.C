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

#include "ZyOverstrainRIA.H"
#include "addToRunTimeSelectionTable.H"
#include "offbeatTime.H"

#include "fvm.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "globalFieldLists.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ZyOverstrainRIA, 0);
    addToRunTimeSelectionTable(failureModel, ZyOverstrainRIA, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ZyOverstrainRIA::ZyOverstrainRIA
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    failureModel(mesh, dict),
    hoopStrainLimit_
    (
        IOobject
        (
            "hoopStrainLimitRIA",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("hoopStrainLimitRIA", dimless, 0.053),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    epsDotYY_
    (
        IOobject
        (
            "epsDotYY",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("epsDotYY", dimless, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    K1_
    (   
        IOobject
        (
            "K1",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("K1", dimensionSet(1,-0.5,-2,0,0,0,0), 0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    CFP_
    (
        IOobject
        (
            "CFP",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("CFP", dimless, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    CH_Value_(readScalar(dict.lookup("CH"))),
    phiValue_(readScalar(dict.lookup("phi"))),
    crackLengthValue_(readScalar(dict.lookup("crackLength")))
{    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ZyOverstrainRIA::~ZyOverstrainRIA()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::ZyOverstrainRIA::isFailed
(
    const labelList& addr
)
{

    // Assign mapper object to pointer
    if(mapper_ == nullptr)
    {
        mapper_ = &mesh_.lookupObject<sliceMapper>("sliceMapper");
    }

    phiValue_ *= 1e25;
    
    // Constant reference to the epsilonP field
    const volSymmTensorField& epsilonP(mesh_.lookupObject<volSymmTensorField>("epsilonP"));
    const volSymmTensorField& epsilonEl(mesh_.lookupObject<volSymmTensorField>("epsEl"));
    const volSymmTensorField epsilonElOld = epsilonEl.oldTime();

    // Constant reference to the T field
    const volScalarField& Temp(mesh_.lookupObject<volScalarField>("T"));

    // Constant reference to DEpsilonP field
    const volSymmTensorField& DEpsP(mesh_.lookupObject<volSymmTensorField>("DEpsilonP"));
    const volSymmTensorField& DEpsEl = epsilonEl - epsilonElOld;

    // Sigma field
    const volSymmTensorField& sig(mesh_.lookupObject<volSymmTensorField>("sigma"));

    scalar currentDeltaT = mesh_.time().deltaT().value();

    forAll(addr, addrI){
        const label cellI = addr[addrI];
        scalar f1 = 1;
        scalar f2 = 1;
        scalar f3 = 1;
        scalar epsf0 = 2.82 + 1.22e-2 * Temp[cellI];

        epsDotYY_[cellI] = DEpsP[cellI].component(symmTensor::YY) * (1/currentDeltaT)
        + DEpsEl[cellI].component(symmTensor::YY) * (1/currentDeltaT) ;

        if(epsDotYY_[cellI] < SMALL){
            f1 = 0.046-0.031*log10(SMALL);
        }
        else if(epsDotYY_[cellI] <= 1){
            f1 = 0.046 - 0.31*log10(epsDotYY_[cellI]);
        }
        else{
            f1 = 0.046;
        }

        scalar gamma = 0;
        if(epsDotYY_[cellI] > SMALL){
            gamma = 6.52e-4 + 2.21e-3 * (6 + log10(epsDotYY_[cellI]))*(1 -  tanh((Temp[cellI] - 298)/8.5));
        }
        if(epsDotYY_[cellI] < SMALL){
            gamma = 0;
        }
        else{
            gamma = 6.52e-4 + 2.21e-3 * (6 + log10(epsDotYY_[cellI]))*(1 -  tanh((Temp[cellI] - 298)/8.5));
        }
        //- Terminal solid solubility of hydrogen
        scalar cTSS = 1.28e5 * exp(-35900 / 8.143 *Temp[cellI]);
        //- Hydrogen in excess (amount > to cTSS)
        scalar CHEX = max(0, CH_Value_ - cTSS);
        f2 = 0.01 + 0.99*exp(- gamma * CHEX);
        f3 = 0.05 + 0.95 *exp(-2.89e-24*phiValue_);
        hoopStrainLimit_[cellI] = 3 * f1*f2*f3 / (f1*f2 + f1*f3 + f2*f3) * epsf0 / 100;

        //CFP AND K1 FIELD
        CFP_[cellI] = (epsilonP[cellI].component(symmTensor::YY) + epsilonEl[cellI].component(symmTensor::YY)) 
        / hoopStrainLimit_[cellI];
        
        K1_[cellI] = sig[cellI].component(symmTensor::YY) 
        *  sqrt(constant::mathematical::pi * crackLengthValue_);
    }
    // Set failure switch to false
    failed_ = false;

    // Initialize epsilonCPCyl as a copy of epsilonP and * to zero
    volSymmTensorField epsilonPCyl("epsilonPCyl", epsilonP);
    epsilonPCyl *= scalar(0);

    // Transform the epsilonP tensor in cylindrical coordinates
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];

        tensor R = tensor::zero;

        const vector& Ci = mesh_.C()[cellI];
        const scalar  radius = sqrt(Ci[0]*Ci[0] + Ci[1]*Ci[1]);

        R.xx() =  Ci[0]/radius;
        R.xy() =  Ci[1]/radius;
        R.yy() =  Ci[0]/radius;
        R.yx() = -Ci[1]/radius;
        R.zz() = 1;

        epsilonPCyl[cellI] = symm( R & epsilonP[cellI] & R.T() );
    }

    // Calculate the slice-average values of epsilonCreep and hoop limit
    const symmTensorField& epsilonPSliceAvg(mapper_->sliceAverage(epsilonPCyl));
    const scalarField& hoopStrainLimitSliceAvg(mapper_->sliceAverage(hoopStrainLimit_));

    const offbeatTime& userRunTime(
        refCast<const offbeatTime>(mesh_.time()));

    // Check if the avg epsilonCreep per slice overcmes the hoop strain
    forAll(epsilonPSliceAvg, i)
    {
        const symmTensor e = epsilonPSliceAvg[i];
        const scalar criticalValue = hoopStrainLimitSliceAvg[i];

        if ( e.component(symmTensor::YY) >= criticalValue )
        {
            if( !mesh_.foundObject<volScalarField>("failedMaterial"))
            {
                // Initialize failedMaterial field (from parent class "failureModel") 
                initializeFailedMaterialField();
            }

            failed_ = true;

            // Print out information about failure
            Info << "Failed at time " <<  userRunTime.userTime() << userRunTime.unit() 
                 << ", in slice n. " << i 
                 << ", with slice-avg. creep hoop strain = " 
                 << e.component(symmTensor::YY)
                 << endl;

            // Set the field failedMaterial_ to 1 in correspondance of the slice
            const labelList addrFailed = mapper_->sliceAddrList()[i];
            forAll( addrFailed, addrI )
            {
                const label cellI = addrFailed[addrI];
                failedMaterial_()[cellI] = 1;
            }
        }
    }

    return failed_;

}


// ************************************************************************* //
