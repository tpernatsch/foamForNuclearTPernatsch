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

#include "YoungModulusMatproUO2.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(YoungModulusMatproUO2, 0);
    addToRunTimeSelectionTable
    (
        YoungModulusModel, 
        YoungModulusMatproUO2, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::YoungModulusMatproUO2::YoungModulusMatproUO2
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    YoungModulusModel(mesh, dict, defaultModel),       
    densityFrac_(dict.lookupOrDefault("densityFraction", 0.95)),
    par1(2.334e11),
    par2(2.752),
    par3(1.0915e-4),
    perturb(1.0),
    isotropicCracking_(dict.lookup("isotropicCracking")),
    isotropicCrackingType_("Barani"),
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

        par1 = YoungModulusDict.lookupOrDefault<scalar>("par1", 2.334e11);
        par2 = YoungModulusDict.lookupOrDefault<scalar>("par2", 2.752);
        par3 = YoungModulusDict.lookupOrDefault<scalar>("par3", 1.0915e-4);

        perturb = YoungModulusDict.lookupOrDefault<scalar>("perturb", 1.0);
    }

    if(isotropicCracking_)
    {
        dict.lookup("nCracksMax") >> nCracksMax_;
        isotropicCrackingType_ = dict.lookupOrDefault<word>("isotropicCrackingType", "Barani");

        if(isotropicCrackingType_ != "Barani"
            and isotropicCrackingType_ != "JankusWeeks")
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

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::YoungModulusMatproUO2::~YoungModulusMatproUO2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::YoungModulusMatproUO2::correct
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
        else
        {
            forAll(addr, addrI)
            {
                const label& cellI = addr[addrI];
                
                nCracks_.ref()[cellI] = nCracksMax_;
            }
        }
    }



    forAll(addr, i)
    {   
        const label cellI = addr[i];

        const scalar Ti = T[cellI];
        const scalar nui = 0.316;

        const scalar nominalValue = par1*(1 - par2*(1-densityFrac_))*(1 - par3*Ti);
        
        const scalar fnu = (2/3.0)*(2 - nui)/(2 + nui)/(1 - nui);

        if(isotropicCracking_)
        {
            if(isotropicCrackingType_ == "Barani")
            {
                sf[cellI] =  pow(fnu, nCracks_.internalField()[cellI])*nominalValue*perturb;
            }
            else
            {
                sf[cellI] =  pow(2.0/3.0, nCracks_.internalField()[cellI])*nominalValue*perturb;                
            }            
        }
        else
        {
            sf[cellI] =  nominalValue*perturb;
        }
    }
}

// ************************************************************************* //
