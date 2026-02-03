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

#include "conductivityMatproUO2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conductivityMatproUO2, 0);
    addToRunTimeSelectionTable
    (
        conductivityModel, 
        conductivityMatproUO2, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityMatproUO2::conductivityMatproUO2
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    conductivityModel(mesh, dict, baseDict, defaultModel),       
    burnupName_(dict.lookupOrDefault<word>("burnupName", "Bu" )),
    Bu_(nullptr),
    densityFrac_(mesh.nCells(), dict.lookupOrDefault<scalar>
    (
        "densityFraction",
        readScalar(baseDict.lookup("densityFraction"))
    )),
    GdContent_
    (
        dict.lookupOrDefault<scalar>
        (
            "GdContent",
            baseDict.lookupOrDefault<scalar>("GdContent", 0.0)
        )
    ),
    par1(0.0452),
    par2(0.000246),
    par3(0.00187),
    par4(1.1599),
    par5(0.9),
    par6(0.04),
    par7(0.038),
    par8(0.28),
    par9(396.0),
    par10(6380.0),
    par11(3.5e9),
    par12(2.0), 
    par13(16360.0)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'conductivity' dictionary
    if((dict.found("conductivity")) || (dict.dictName() == "conductivity"))
    {
        const dictionary& conductivityDict = dict.found("conductivity")?
        dict.subDict("conductivity"):
        dict;

        par1 = conductivityDict.lookupOrDefault<scalar>("par1", 0.0452);
        par2 = conductivityDict.lookupOrDefault<scalar>("par2", 0.000246);
        par3 = conductivityDict.lookupOrDefault<scalar>("par3", 0.00187);
        par4 = conductivityDict.lookupOrDefault<scalar>("par4", 1.1599);
        par5 = conductivityDict.lookupOrDefault<scalar>("par5", 0.9);
        par6 = conductivityDict.lookupOrDefault<scalar>("par6", 0.04);
        par7 = conductivityDict.lookupOrDefault<scalar>("par7", 0.038);
        par8 = conductivityDict.lookupOrDefault<scalar>("par8", 0.28);
        par9 = conductivityDict.lookupOrDefault<scalar>("par9", 396.0);
        par10 = conductivityDict.lookupOrDefault<scalar>("par10", 6380.0);
        par11 = conductivityDict.lookupOrDefault<scalar>("par11", 3.5e9);
        par12 = conductivityDict.lookupOrDefault<scalar>("par12", 2.0); 
        par13 = conductivityDict.lookupOrDefault<scalar>("par13", 16360.0);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductivityMatproUO2::~conductivityMatproUO2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::conductivityMatproUO2::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{   
    if(Bu_ == nullptr)
    {        
        Bu_ = &mesh_.lookupObject<volScalarField>(burnupName_);
    }
 
    const scalarField& Bu = Bu_->internalField();

    //- To include porosity evolution, the new cell volume is considered to be
    //- equal to original volume times (1+3epsilonSwelling + 3epsilonDensification)
    scalarField volumeCorr(addr.size(), 1);

    if(mesh_.foundObject<volScalarField>("Intragranular_gas_swelling"))
    {
        const scalarField& Intragranular_gas_swelling = 
        mesh_.lookupObject<volScalarField>("Intragranular_gas_swelling").internalField();
        const scalarField& Intergranular_gas_swelling = 
        mesh_.lookupObject<volScalarField>("Intergranular_gas_swelling").internalField();

        forAll(addr, i)
        {  
            const label cellI = addr[i];  
            volumeCorr[i] += 3*Intragranular_gas_swelling[cellI]/3;
            volumeCorr[i] += 3*Intergranular_gas_swelling[cellI]/3;
        }
    }

    if(mesh_.foundObject<volScalarField>("epsilonDensification"))
    {
        const scalarField& epsilonDensification = 
        mesh_.lookupObject<volScalarField>("epsilonDensification").internalField();

        forAll(addr, i)
        {  
            const label cellI = addr[i];  
            volumeCorr[i] += 3*epsilonDensification[cellI];
        }
    }

    // Add correction for porosity evolution due to porosity migration
    if(mesh_.foundObject<volScalarField>("porosity"))
    {
        const scalarField& p = 
        mesh_.lookupObject<volScalarField>("porosity").internalField();

        forAll(addr, i)
        {  
            const label cellI = addr[i];
            
            // Set cutoff at 95% porosity
            densityFrac_[i] = max(1-p[cellI], 0.05);
        }
    }

    //- Correct conductivity for every cell in addr
    forAll(addr, i)
    {   
        const label cellI = addr[i];              

        // Convert burnup to MWd/kgU
        const scalar Bui = Bu[cellI]/1000/0.881;

        const scalar Ti = T[cellI];
            
        double term0 = par1 ;
        double term1 = par2*Ti;
        double term2 = par3*Bui; 
        double term3 = par4*GdContent_; 
        double term4 = (1-par5*exp(-par6*Bui))*par7*pow(Bui,par8);    
        double term5 = 1/(1+par9*exp(-par10/Ti));
        double term6 = (par11/pow(Ti,par12))*exp(-par13/Ti);
        
        double k95 = 1/(term0 + term1 + term2 + term3 + term4*term5) + term6;    
        
        // Update initial density fraction with correction due to porosity 
        // evolution
        const scalar denFracNew = densityFrac_[i]/volumeCorr[i];
        const scalar nominalValue = k95*1.0789*denFracNew/(1+0.5*(1-denFracNew));

        sf[cellI] = nominalValue;
    }
}

// ************************************************************************* //
