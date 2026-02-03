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
#include "isotropicCracking.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"
#include "fvm.H"
#include "fvc.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "Time.H"

namespace Foam
{
    defineTypeNameAndDebug(isotropicCracking, 0);
    addToRunTimeSelectionTable
    (
        damageModel,
        isotropicCracking,
        dictionary
    );
}

Foam::isotropicCracking::isotropicCracking
(
    const fvMesh& mesh,
    const dictionary& damageDict,
    const dictionary& baseDict
)
:   
    damageModel(mesh, damageDict, baseDict),
    mapper_(nullptr),
    nCracks_
    (
        createOrLookup<scalar>(mesh, "nCracks", dimless, 0.0,
                zeroGradientFvPatchField<scalar>::typeName)
    ),
    isotropicCrackingType_(damageDict.lookupOrDefault<word>("isotropicCrackingType", "Barani")),
    nCracksMax_(12),
    nCracksFixed_(1)
{
    if (isotropicCrackingType_ == "Barani")
    {
        if(damageDict.found("nCracksMax"))
        {
            damageDict.lookup("nCracksMax") >> nCracksMax_;
        }
    }
    else if (isotropicCrackingType_ == "JankusWeeks")
    {
        damageDict.lookup("nCracksMax") >> nCracksFixed_;
    }
    else
    {                
        FatalErrorInFunction()
            << "Isotropic cracking model " << isotropicCrackingType_ << " does not"
            << " exist. Only choices are:" << nl 
            << "- Barani (default)" << nl
            << "- JankusWeeks"
            << abort(FatalError);
    }
}

Foam::isotropicCracking::~isotropicCracking()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::isotropicCracking::soften
(
    scalarField& E,
    scalarField& nu,
    const labelList& addr
)
{
    // Calculate number of cracks
    if(isotropicCrackingType_ == "Barani")
    {
        if(mapper_ == nullptr)
        {
            mapper_ = &mesh_.lookupObject<sliceMapper>("sliceMapper");
        }

        scalar nSlices(mapper_->nSlices());
        if(nSlices<=0)
        {
            FatalErrorInFunction()
            << "Number of slices cannot be less than 1." << nl
            << "Check that a sliceMapper different than \"none\" is used "
            << "or check the slice definition."
            << abort(FatalError);
        }

        if(!mesh_.foundObject<volScalarField>("lhgr"))
        {
            FatalErrorInFunction()
                << "Linear heat generation rate field (lhgr) not found." << nl
                << "Check if the selected heat source model defines the lhgr "
                << "as this field is needed by the isotropic cracking model." << nl
                << "Alternatively, deactivate the isotropic cracking model."
                << exit(FatalError);
        }

        // Reference to lhgr (in W/m)
        const volScalarField& lhgr(
            mesh_.lookupObject<volScalarField>("lhgr"));

        // Calculate the slice-average values of lhgr (in W/m)
        const scalarField& lhgrSliceAvg(mapper_->sliceAverage(lhgr));
    
        //- Linear Heat Rate in kW/m for first crack
        scalar lhr0(5);

        //- Number of first cracks
        scalar n0(1);

        //- Asymptotic number of  cracks
        scalar nInf(12);

        //- Exponential period in kW/m
        scalar tau(21);

        forAll(addr, addrI)
        {
            const label& cellI = addr[addrI];
            const label& sliceID = mapper_->sliceID()[cellI];

            //- Linear Heat Rate in kW/m
            scalar lhrSlice(lhgrSliceAvg[sliceID]);

            scalar nCracks
            (
                n0 + (nInf - n0)*(1 - exp(- (lhrSlice/1000 - lhr0)/tau))
            );

            nCracks_.ref()[cellI]= 
            min
            (
                max(nCracks_.internalField()[cellI], nCracks),
                nCracksMax_
            );
        }
    }
    else if (isotropicCrackingType_ == "JankusWeeks")
    {
        forAll(addr, addrI)
        {
            const label& cellI = addr[addrI];
            
            nCracks_.ref()[cellI] = nCracksFixed_;
        }
    }
    else
    {                
        FatalErrorInFunction()
            << "Isotropic cracking model " << isotropicCrackingType_ << " does not"
            << " exist. Only choices are:" << nl 
            << "- Barani (default)" << nl
            << "- JankusWeeks"
            << abort(FatalError);
    }
    
    // Soften E and nu
    forAll(addr, i)
    {   
        const label cellI = addr[i];

        scalar& nui = nu[cellI];
        scalar& Ei  = E[cellI];
        const scalar& nCracksi = nCracks_.internalField()[cellI];

        if (isotropicCrackingType_ == "Barani")
        {
            const scalar fnu = (2/3.0)*(2 - nui)/(2 + nui)/(1 - nui);
            Ei =  pow(fnu, nCracksi)*Ei;
            nui = nui/(pow(2, nCracksi) + (pow(2, nCracksi) - 1)*nui);
        }
        else if (isotropicCrackingType_ == "JankusWeeks")
        {
            Ei =  pow(2.0/3.0, nCracksi);
            nui =  pow(0.5, nCracksi)*nui;            
        }    
        else
        {                
            FatalErrorInFunction()
                << "Isotropic cracking model " << isotropicCrackingType_ << " does not"
                << " exist. Only choices are:" << nl 
                << "- Barani (default)" << nl
                << "- JankusWeeks"
                << abort(FatalError);
        }
    }
}