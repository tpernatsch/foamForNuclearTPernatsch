
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

#include "conductivityUPuO2Philipponneau.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conductivityUPuO2Philipponneau, 0);
    addToRunTimeSelectionTable
    (
        conductivityModel, 
        conductivityUPuO2Philipponneau, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityUPuO2Philipponneau::conductivityUPuO2Philipponneau
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
    B_(2.885e-4),
    C_(76.38e-12),
    OM_(mesh_.lookupObject<volScalarField>("oxygenMetalRatio"))
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'conductivity' dictionary
    if((dict.found("conductivity")) || (dict.dictName() == "conductivity"))
    {
        const dictionary& conductivityDict = dict.found("conductivity")?
        dict.subDict("conductivity"):
        dict;

        B_ = conductivityDict.lookupOrDefault<scalar>("B", 2.885e-4);
        C_ = conductivityDict.lookupOrDefault<scalar>("C", 76.38e-12);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductivityUPuO2Philipponneau::~conductivityUPuO2Philipponneau()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::conductivityUPuO2Philipponneau::correct
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
    
    // Reference to Bu field in MWd/Mton 
    const scalarField& Bu = Bu_->internalField();

    // Convert Bu in GWd/ton and then divide by 9.5 to get it in %FIMA 
    const scalarField BuFIMA = Bu * 1e-3 / 9.5 / 100;

    // Reference to porosity field
    const volScalarField& porosity(
        mesh_.lookupObject<volScalarField>("porosity"));

    // Correct conductivity for every cell in addr
    forAll(addr, i)
    {   
        const label cellI = addr[i];              

        // Fractional burnup of cellI
        const scalar tau = BuFIMA[cellI];

        // Porosity of this cell - with cutoff
        const scalar& poroI = min(porosity[cellI], 0.95);

        // Copmute deviation from stoichiometry
        const scalar x = 2-OM_[cellI];

        // Compute A
        const scalar A = 1.528 * Foam::sqrt(x+0.00931) - 0.1055 + 0.44*tau;

        // Compute alpha
        const scalar alpha = 1/0.864 * (1-poroI)/(1 + 2*poroI);

        // Cell temperature
        const scalar Ti = T[cellI];

        // Computing thermal conductivity
        const scalar nominalValue = ( 1/(A + B_*Ti) + C_*pow(Ti,3.0) ) * alpha;

        sf[cellI] = nominalValue;
    }
}

// ************************************************************************* //
