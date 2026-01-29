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

#include "YoungModulusMatproUPuO2.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(YoungModulusMatproUPuO2, 0);
    addToRunTimeSelectionTable
    (
        YoungModulusModel, 
        YoungModulusMatproUPuO2, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::YoungModulusMatproUPuO2::YoungModulusMatproUPuO2
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    YoungModulusModel(mesh, dict, baseDict, defaultModel),       
    densityFrac_
    (
        dict.lookupOrDefault<scalar>
        (
            "densityFraction",
            baseDict.lookupOrDefault<scalar>("densityFraction", 0.945)
        )
    ),
    N_U_(mesh_.lookupObject<volScalarField>("N_U")),
    N_Pu_(mesh_.lookupObject<volScalarField>("N_Pu")),
    N_Am_(nullptr),
    par1(2.334e11),
    par2(2.752),
    par3(1.0915e-4),
    B(0),    
    OM_
    (
        dict.lookupOrDefault<scalar>
        (
            "oxygenMetalRatio",
            baseDict.lookupOrDefault<scalar>("oxygenMetalRatio", 2.0)
        )
    ),
    x_(2 - OM_),
    isotropicCracking_(dict.found("damageModel")? 
        false : dict.lookupOrDefault<bool>("isotropicCracking", false)),
    isotropicCrackingType_("Barani"),
    nCracks_
    (
        createOrLookup<scalar>(mesh, "nCracks", dimless, 0.0,
                zeroGradientFvPatchField<scalar>::typeName)
    ),
    nCracksMax_(1),
    mapper_(nullptr)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'YoungModulus' dictionary
    if((dict.found("YoungModulus")) || (dict.dictName() == "YoungModulus"))
    {
        const dictionary& YoungModulusDict = dict.found("YoungModulus")?
        dict.subDict("YoungModulus"):
        dict;

        par1 = YoungModulusDict.lookupOrDefault<scalar>("par1", 2.334e11);
        par2 = YoungModulusDict.lookupOrDefault<scalar>("par2", 2.752);
        par3 = YoungModulusDict.lookupOrDefault<scalar>("par3", 1.0915e-4);
    }

    if(isotropicCracking_)
    {        
        WarningIn("YoungModulusMatproUPuO2(const fvMesh&, const dictionary&, const word)")
            << "The use of 'isotropicCracking' model for this class is still supported "
            << "but is deprecated." << nl 
            << "Please use the 'damageModel' class instead." << nl << endl;
            
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

    if( x_ < 0 )
    {
        // Hypostechiometric fuel
        B = 1.75;
    }
    else
    {
        // Hyperstechiometric fuel
        B = 1.35;
    }

    // Check if minor actinides have been initialized by materials
    if( mesh_.foundObject<volScalarField>("N_Am") )
    {
        N_Am_ = &mesh_.lookupObjectRef<volScalarField>("N_Am");
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::YoungModulusMatproUPuO2::~YoungModulusMatproUPuO2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::YoungModulusMatproUPuO2::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{ 
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

        // Number densities for heavy metals in fuel
        const scalar N_U = N_U_[cellI];
        const scalar N_Pu = N_Pu_[cellI];
        const scalar N_Am = (N_Am_ == nullptr) ? 0.0 : N_Am_->internalField()[cellI];

        // Compute concentrations atoms/atomsHM
        const scalar conc_at_Pu = N_Pu / (N_Pu + N_Am + N_U);    

        // Conversion factors from atoms/atomsHM -> mass/mass_fuel
        // Approximate conversion factors are used, computed in this way:
        // 
        // For Pu:
        // at.Pu/at.HM = k_Pu * m_Pu/m_HM
        // Assuming: 
        // MM_Pu ~ 239 g/mol, OM_ratio = 2, MM_fuel ~ 238.5 g/mol, MM_O = 16:
        // => k_Pu = (MM_HM + OM_ratio * MM_O) / MM_Pu
        const scalar k_Pu = 1.13;

        // Assign concentrations in mass/mass_fuel
        const scalar c_w_Pu = conc_at_Pu / k_Pu;

        const scalar nui = 0.276;

        const scalar ES = par1*(1 - par2*(1-densityFrac_))*(1 - par3*Ti);

        const scalar nominalValue = ES*exp(-B*x_)*(1+0.05*c_w_Pu);
        
        const scalar fnu = (2/3.0)*(2 - nui)/(2 + nui)/(1 - nui);

        if(isotropicCracking_)
        {
            if(isotropicCrackingType_ == "Barani")
            {
                sf[cellI] =  pow(fnu, nCracks_.internalField()[cellI])*nominalValue;
            }
            else
            {
                sf[cellI] =  pow(2.0/3.0, nCracks_.internalField()[cellI])*nominalValue;                
            }            
        }
        else
        {
            sf[cellI] =  nominalValue;
        }
    }
}

// ************************************************************************* //
