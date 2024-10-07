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

#include "YoungModulusSckCenUPuO2.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(YoungModulusSckCenUPuO2, 0);
    addToRunTimeSelectionTable
    (
        YoungModulusModel, 
        YoungModulusSckCenUPuO2, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::YoungModulusSckCenUPuO2::YoungModulusSckCenUPuO2
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    YoungModulusModel(mesh, dict, defaultModel),       
    y_(2.0-dict.lookupOrDefault<scalar>("oxygenToMetalRatio",2.0)),
    ratioPuMetal_(readScalar(dict.lookup("PuWeightFraction"))),
    dEdy_(),
    E1_(219.12),
    E2_(-0.0154),
    E3_(-9.0e-6),
    EUO2_(218.74),
    EPuO2_(249.45),
    perturb(1.0),
    isotropicCracking_(dict.lookup("isotropicCracking")),
    nCracks_
    (
        createOrLookup<scalar>(mesh, "nCracks", dimless, 0.0,
                zeroGradientFvPatchField<scalar>::typeName)
    ),
    nCracksMax_(1),
    mapper_(nullptr)
{
    if(dict.found("YoungModulus"))
    {
        const dictionary& YoungModulusDict = dict.subDict("YoungModulus");

        E1_ = YoungModulusDict.lookupOrDefault<scalar>("E1", 219.12);
        E2_ = YoungModulusDict.lookupOrDefault<scalar>("E2", -0.0154);
        E3_ = YoungModulusDict.lookupOrDefault<scalar>("E3", -9.0e-6);
        EUO2_  = YoungModulusDict.lookupOrDefault<scalar>("EUO2", 218.74);
        EPuO2_ = YoungModulusDict.lookupOrDefault<scalar>("EPuO2", 249.45);

        perturb = YoungModulusDict.lookupOrDefault<scalar>("perturb", 1.0);
    }

    if(isotropicCracking_)
    {
        dict.lookup("nCracksMax") >> nCracksMax_;
    }

    if( y_ <= 0.037 )
    {
        dEdy_ = -586;
    }
    else
    {
        dEdy_ = -126.59;
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::YoungModulusSckCenUPuO2::~YoungModulusSckCenUPuO2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::YoungModulusSckCenUPuO2::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{    
    // const scalarField& nu
    // (
    //     mesh_.lookupObject<volScalarField>("nu").internalField()
    // );
        
    if(isotropicCracking_)
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

    forAll(addr, i)
    {   
        const label cellI = addr[i];

        const scalar Ti = T[cellI];

        const scalar Emox0 = (1-ratioPuMetal_)*EUO2_ + ratioPuMetal_*EPuO2_;
        
        const scalar Emox = Emox0 + y_*dEdy_;

        const scalar fT   = E1_ + E2_*Ti     + E3_*pow(Ti,2);
        const scalar f273 = E1_ + E2_*273.00 + E3_*pow(273.00,2);

        const scalar nominalValue = Emox*fT/f273*1e9;

        const scalar nui = 0.276;
        const scalar fnu = (2/3.0)*(2 - nui)/(2 + nui)/(1 - nui);
        
        if(isotropicCracking_)
        {
            sf[cellI] =  pow(fnu, nCracks_.internalField()[cellI])*nominalValue*perturb;
        }
        else
        {
            sf[cellI] =  nominalValue*perturb;
        }
    }
}

// ************************************************************************* //
