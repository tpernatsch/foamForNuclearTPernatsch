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

#include "densificationFrapcon.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(densificationFrapcon, 0);
    addToRunTimeSelectionTable(densificationModel, densificationFrapcon, dictionary);

    const char* densificationFrapcon::group_ = 
        ("behavioralModels::densification::UO2Frapcon");

    defineParameter(densificationFrapcon, F_epsilonDensification_, 
        "F_epsilonDensification", (dimless), 1.0);
        
    defineParameter(densificationFrapcon, delta_epsilonDensification_, 
        "delta_epsilonDensification", (dimless), 0.0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::densificationFrapcon::densificationFrapcon
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict
)
:
    densificationModel(mesh, dict, baseDict),
    burnupName_(dict.lookupOrDefault<word>("burnupName", "Bu" )),
    Bu_(nullptr),
    densityFrac_
    (
        dict.lookupOrDefault<scalar>
        (
            "densityFraction",
            baseDict.lookupOrDefault<scalar>("densityFraction", 0.95)
        )
    ),
    resintDensChange_(readScalar(dict.lookup("resinteringDensityChange"))),
    Tsinter_(dict.lookupOrDefault<scalar>("Tsintering", 1800)),
    par1(22.2),
    par2(1453),
    par3(66.6),
    par4(2.0),
    par5(35),
    par6(3.0)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'densification' dictionary
    if((dict.found("densification")) || (dict.dictName() == "densification"))
    {
        const dictionary& densificationDict = dict.found("densification")?
        dict.subDict("densification"):
        dict;

        par1 = densificationDict.lookupOrDefault<scalar>("par1", 22.2);
        par2 = densificationDict.lookupOrDefault<scalar>("par2", 1453);
        par3 = densificationDict.lookupOrDefault<scalar>("par3", 66.6);
        par4 = densificationDict.lookupOrDefault<scalar>("par4", 2.0);
        par5 = densificationDict.lookupOrDefault<scalar>("par5", 35.0);
        par6 = densificationDict.lookupOrDefault<scalar>("par6", 3.0);
    }  
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::densificationFrapcon::~densificationFrapcon()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::densificationFrapcon::correct
(
    const scalarField& T, 
    const labelList& addr
)
{
    // Check which F and delta should be used
    const dimensionedScalar& F = 
    isInUserParameters(this->group(), "F_epsilonDensification") ?
    this->F_epsilonDensification() : 
    densificationModel::F_epsilonDensification();

    const dimensionedScalar& delta = 
    isInUserParameters(this->group(), "delta_epsilonDensification") ?
    this->delta_epsilonDensification() : 
    densificationModel::delta_epsilonDensification();

    if(Bu_ == nullptr)
    {        
        Bu_ = &mesh_.lookupObject<volScalarField>(burnupName_);
    }

    const scalarField oldTimeDensification_ = 
    epsilonDensification_.oldTime().internalField();
    
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        
        double dLOverL_max(0.0);
        
        // Convert burnup to MWd/kgU
        const scalar Bui = Bu_->internalField()[cellI]/1000/0.881;
        
        if(resintDensChange_ == 0)
        {
            if(T[cellI] < 950)
            {
                dLOverL_max = 
                -par1*(1 - densityFrac_)/(Tsinter_ - par2);
            }
            // Interpolation region to smooth transition (arbitrary 100K wide range)
            else if(T[cellI] < 1050)
            {
                scalar dLOverL_950 = 
                -par1*(1 - densityFrac_)/(Tsinter_ - par2);

                scalar dLOverL_1050 = 
                -par1*(1 - densityFrac_)/(Tsinter_ - par2);

                dLOverL_max = dLOverL_950 + 
                (T[cellI] - 950)/(1050 - 950)*(dLOverL_1050 - dLOverL_950);  
            }
            else
            {
                dLOverL_max = 
                -par3*(1 - densityFrac_)/(Tsinter_ - par2);            
            }
        }
        else
        {
            dLOverL_max = 
            -resintDensChange_/(densityFrac_*100 + resintDensChange_)/3;

            if(T[cellI] < 950)
            {
                dLOverL_max = -(0.0015*resintDensChange_*10960/100)/100;
            }
            // Interpolation region to smooth transition (arbitrary 100K wide range)
            else if(T[cellI] < 1050)
            {
                scalar dLOverL_950 = 
                -(0.0015*resintDensChange_*10960/100)/100;;

                scalar dLOverL_1050 = 
                -(0.00285*resintDensChange_*10960/100)/100;

                dLOverL_max = dLOverL_950 + 
                (T[cellI] - 950)/(1050 - 950)*(dLOverL_1050 - dLOverL_950);  
            }
            else
            {
                dLOverL_max = -(0.00285*resintDensChange_*10960/100)/100;
            }
        }
        
        // Iterate to find B such that dL/L=0 at Bui=0
        double B0, B = 1;
        int nIters = 0;
        B = 1;
        do
        {
            B0 = B;
            B = -log(-par4*exp(-par5*B) - dLOverL_max)/par6;
            
            nIters++;
            
            if (nIters > 10)
            {
            FatalErrorIn("Foam::densification(...)")
                    << "Solution failed to converge."
                    << Foam::abort(FatalError);
            }
        } while (abs(B-B0) > 1e-6);


        scalar nominalValue = 
        dLOverL_max + exp(-par6*(Bui+B)) + par4*exp(-par5*(Bui+B));  

        // Perturb epsilon densification for sensitivity analysis
        nominalValue = nominalValue * F.value() + delta.value();                          
                                
        //- Check if new densification value is lower than previous one           
        //- Keep the smallest value
        scalarField& densificationI = epsilonDensification_.ref();

        densificationI[cellI] = 
        min(nominalValue, oldTimeDensification_[cellI] );
    }


}


// ************************************************************************* //
