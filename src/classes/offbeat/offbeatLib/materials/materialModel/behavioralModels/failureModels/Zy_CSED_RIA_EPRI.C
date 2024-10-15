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

#include "Zy_CSED_RIA_EPRI.H"
#include "addToRunTimeSelectionTable.H"
#include "offbeatTime.H"

#include "fvm.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "globalFieldLists.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Zy_CSED_RIA_EPRI, 0);
    addToRunTimeSelectionTable(failureModel, Zy_CSED_RIA_EPRI, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Zy_CSED_RIA_EPRI::Zy_CSED_RIA_EPRI
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    failureModel(mesh, dict),
    SED_
    (
        IOobject
        (
            "SED",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("SED", dimless, 0.0),
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
    CH_Value_(readScalar(dict.lookup("CH"))),
    crackLengthValue_(readScalar(dict.lookup("crackLength")))
{    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Zy_CSED_RIA_EPRI::~Zy_CSED_RIA_EPRI()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::Zy_CSED_RIA_EPRI::isFailed
(
    const labelList& addr
)
{

    // Assign mapper object to pointer
    if(mapper_ == nullptr)
    {
        mapper_ = &mesh_.lookupObject<sliceMapper>("sliceMapper");
    }

    // Constant reference to the epsilonP field
    const volSymmTensorField& DEpsP(mesh_.lookupObject<volSymmTensorField>("DEpsilonP"));
    
    // Epsilon Elastic field
    const volSymmTensorField& epsilonEl(mesh_.lookupObject<volSymmTensorField>("epsEl"));
    //Delta epsilonEl
    const volSymmTensorField& epsilonElOld = epsilonEl.oldTime();
    const volSymmTensorField DEpsEl = epsilonEl - epsilonElOld;

    // Mechanical strain (plastic + elastic)
    const volSymmTensorField DEpsilonM = DEpsP + DEpsEl;
    
    // Sigma field and DSigma
    const volSymmTensorField& sig(mesh_.lookupObject<volSymmTensorField>("sigma"));
    const volScalarField& sigYY(sig.component(symmTensor::YY));

    // SED critical value
    scalar criticalValue = (12.69 * exp(-0.002409 * CH_Value_) + 1.861) * 1e6 ;

    //CFP AND K1 FIELD
    forAll(addr, addrI){
        const label cellI = addr[addrI];
        SED_[cellI] += (DEpsilonM.internalField()[cellI] & sig.internalField()[cellI]) && I;
        CFP_[cellI] = SED_[cellI]/criticalValue;
        K1_[cellI] = sigYY[cellI] * sqrt(constant::mathematical::pi * crackLengthValue_);
    }

    // Set failure switch to false
    failed_ = false;
    // Calculate the slice-average values of SED
    const scalarField& SEDSliceAvg(mapper_->sliceAverage(SED_));

    const offbeatTime& userRunTime(
    refCast<const offbeatTime>(mesh_.time()));

    // Check if the avg SED per slice is greater than the critical value
    forAll(SEDSliceAvg, i)
    {
        const scalar avgCSED = SEDSliceAvg[i];

        if ( avgCSED >= criticalValue)
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
                 << ", with slice-avg. SED = " 
                 << avgCSED
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
