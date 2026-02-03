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

#include "ZircaloyRIAJernkvistModified.H"
#include "addToRunTimeSelectionTable.H"
#include "offbeatTime.H"

#include "fvm.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "globalFieldLists.H"
#include "fundamentalConstants.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ZircaloyRIAJernkvistModified, 0);
    addToRunTimeSelectionTable(failureModel, ZircaloyRIAJernkvistModified, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
Foam::scalar Foam::ZircaloyRIAJernkvistModified::computeS0
(
    const scalar& rimThickness
)
{
    return min(max(0.05, 1 - (1-0.05)/70 * (rimThickness-50)),1);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ZircaloyRIAJernkvistModified::ZircaloyRIAJernkvistModified
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict
)
:
    failureModel(mesh, dict, baseDict),
#ifdef OPENFOAMFOUNDATION    
    phiValue_(failureModelDict_.lookup<scalar>("phi")),
    CH_(failureModelDict_.lookup<scalar>("meanHydrogenConcentration")),
#elif OPENFOAMESI
    phiValue_(failureModelDict_.get<scalar>("phi")),
    CH_(failureModelDict_.get<scalar>("meanHydrogenConcentration")),
#endif
    useRimThicknessReductionFactor_(failureModelDict_.lookupOrDefault<bool>("rimThicknessReductionFactor", false)),
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
    )
{    
    if(useRimThicknessReductionFactor_)
    {
    #ifdef OPENFOAMFOUNDATION    
        rimThickness_ = failureModelDict_.lookup<scalar>("rimThickness");
        RCO_ = failureModelDict_.lookup<scalar>("cladOuterRadius");
        CHRim_ = failureModelDict_.lookup<scalar>("rimHydrogenConcentration");
    #elif OPENFOAMESI
        rimThickness_ = failureModelDict_.get<scalar>("rimThickness");
        RCO_ = failureModelDict_.get<scalar>("cladOuterRadius");
        CHRim_ = failureModelDict_.get<scalar>("rimHydrogenConcentration");
    #endif
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ZircaloyRIAJernkvistModified::~ZircaloyRIAJernkvistModified()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::ZircaloyRIAJernkvistModified::isFailed
(
    const labelList& addr
)
{

	Info << "RIA FAILURE CRITERION!!!!" << endl;
//    Info << addr << endl;
    // not working???
    // const scalar R = constant::physicoChemical::R.value();
    const scalar R = 8.31446261815324;

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

    // Constant reference to DEpsilonP field
    const volSymmTensorField& DEpsP(mesh_.lookupObject<volSymmTensorField>("DEpsilonP"));
    const volSymmTensorField& DEpsEl = epsilonEl - epsilonElOld;

    // Constant reference to the T field
    const volScalarField& T(mesh_.lookupObject<volScalarField>("T"));
    const volScalarField& TOld = T.oldTime();
    const volScalarField dT = T - TOld;

    // Sigma field
    const volSymmTensorField& sig(mesh_.lookupObject<volSymmTensorField>("sigma"));

    scalar currentDeltaT = mesh_.time().deltaT().value();
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];

        const scalar hoopStrainRate = (DEpsP[cellI].component(symmTensor::YY)
        +(epsilonEl[cellI] - epsilonElOld[cellI]).component(symmTensor::YY)) / currentDeltaT;

        //- Base value of the critical hoop strain
        scalar epsf0 = (2.82 + 1.22e-2 * T[cellI]) / 100;

        //- Ductility reduction factors
        //- f1: effect of strain rate, f2: effect of hydrogen, f3: effect of irradiation,
        scalar f1 = 1;
        scalar f2 = 1;
        scalar f3 = 1;
        //- f1 calculation
        if(mag(hoopStrainRate) <= 1)
        {
            f1 = 0.046 - 0.31*log10(mag(hoopStrainRate));
            f1=min(1, f1);
        }
        else
        {
            f1 = 0.046;
        }

        //- gamma coefficient
        scalar gamma = 6.52e-4 + 2.21e-3 * (6+log10(mag(hoopStrainRate))) * (1 - tanh((T[cellI] - 450 )/50));
        const scalar TSS = 99000 * exp(-34523/R*T[cellI]);
        const scalar Cpp = max(0, CH_ - TSS);

        f2 = min(1, 0.01 + 0.99 * exp(- gamma * Cpp));
        f3 = min(1, 0.05 + 0.95 * exp(-2.89e-24*phiValue_));
        //hoopStrainLimit_[cellI] = 3 * epsf0 * SRecovery * f1*f2*f3 / (f1*f2 + f1*f3 + f2*f3);
        hoopStrainLimit_[cellI] = 3 * epsf0 * f1*f2*f3 / (f1*f2 + f1*f3 + f2*f3);
    }

    //Info << CFP_ << endl;
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
    const scalarField dTSliceAvg(mapper_->sliceAverage(dT));

    // Info << epsilonPSliceAvg << endl;

    const offbeatTime& userRunTime(
        refCast<const offbeatTime>(mesh_.time()));

    // Check if the avg epsilonP per slice overcmes the hoop strain
    forAll(epsilonPSliceAvg, i)
    {
        const symmTensor e = epsilonPSliceAvg[i];
        scalar criticalValue = hoopStrainLimitSliceAvg[i];
        const scalar dTSlice = dTSliceAvg[i];

        //- List of cladd cells of current axialSlice
        const labelList addrSlice = mapper_->sliceAddrList()[1];

        if(useRimThicknessReductionFactor_)
        {   
            // get z of the current axialSlice
            const label cellZZ = addrSlice[0];
            const scalarField Cz = mesh_.C().component(vector::Z);
            const scalar zzAxialSlice = Cz[cellZZ];
            // rim depth
            const scalar xxRIM = RCO_ - rimThickness_ * 1e-6;
            //- get nearest cell center
            const label cellRIM = mesh_.findCell(point(xxRIM,0, zzAxialSlice));

            //- Get Temperature at blister's depth
            const scalar TRIM = T[cellRIM];
            //Info << "TRim --> " << TRIM << endl;
            
            // H concentration at depth of the rim
            const scalar CssRim = max(CHRim_, 99000 * exp(-34523/R*TRIM));

            const scalar S0 = computeS0(rimThickness_);
            const scalar SRecovery = S0 + CssRim/CHRim_*(1-S0);
            
            criticalValue *= SRecovery;
            forAll(addrSlice, addrI)
            {
                const label cellI = addr[addrI];
                if(dTSlice > 0)
                {
                    CFP_[cellI] = max(e.component(symmTensor::YY)/criticalValue ,0) > CFP_.oldTime().internalField()[cellI] ? 
                    max(e.component(symmTensor::YY)/criticalValue ,0) : CFP_.oldTime().internalField()[cellI];       
                }
            }
        }
        else
        {
            // get z of the current axialSlice
            forAll(addrSlice, addrI)
            {
                const label cellI = addr[addrI];
                if(dTSlice > 0)
                {
                    CFP_[cellI] = max(e.component(symmTensor::YY)/criticalValue ,0) > CFP_.oldTime().internalField()[cellI] ? 
                    max(e.component(symmTensor::YY)/criticalValue ,0) : CFP_.oldTime().internalField()[cellI];        
                }
            }
        }

        if ( e.component(symmTensor::YY) >= criticalValue && dTSlice > 0)
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
