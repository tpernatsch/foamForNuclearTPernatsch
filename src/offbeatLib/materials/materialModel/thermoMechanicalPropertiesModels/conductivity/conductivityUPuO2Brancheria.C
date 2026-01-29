
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

#include "conductivityUPuO2Brancheria.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conductivityUPuO2Brancheria, 0);
    addToRunTimeSelectionTable
    (
        conductivityModel, 
        conductivityUPuO2Brancheria, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityUPuO2Brancheria::conductivityUPuO2Brancheria
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    conductivityModel(mesh, dict, baseDict, defaultModel),
    densityFrac_
    (
        dict.lookupOrDefault<scalar>
        (
            "densityFraction",
            baseDict.lookupOrDefault<scalar>("densityFraction", 0.945)
        )
    ),
    p_(nullptr),
    perturb(1)  
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'conductivity' dictionary
    if((dict.found("conductivity")) || (dict.dictName() == "conductivity"))
    {
        const dictionary& conductivityDict = dict.found("conductivity")?
        dict.subDict("conductivity"):
        dict;
        
        perturb = conductivityDict.lookupOrDefault<scalar>("perturb", 1.0);
    }

}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductivityUPuO2Brancheria::~conductivityUPuO2Brancheria()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::conductivityUPuO2Brancheria::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{
    // init scalar porosity
    scalar pi(1-densityFrac_);

    // Assign porosity pointer
    if (p_ == nullptr)
    {
        if ( mesh_.foundObject<volScalarField>("porosity") )
        {
            p_ = &mesh_.lookupObject<volScalarField>("porosity");
        }
    }

    // Correct conductivity for every cell in addr
    forAll(addr, i)
    {   
        const label cellI = addr[i];

        // Temperature
        const scalar Ti = T[cellI];

        // Porosity
        if (p_ != nullptr) 
        {
            pi = p_->internalField()[cellI];
        }

        // Compute D1 coeff depending on porosity of this cell (cutoff @ 15%)
        scalar D1 = (pi <= 0.05)
        ? 1 - 1.5*pi
        : 1 - min(pi, 0.15) - 10 * pow(min(pi, 0.15), 2.0);

        // Compute conductivity in W/m K
        scalar nominalValue = 
            1e2 * (D1 * (pow(2.88+0.0252*Ti, -1) + 5.83e-13*pow(Ti,3.0)));
        
        // Assign k to scalar field
        sf[cellI] = nominalValue*perturb;
    }
}

// ************************************************************************* //
