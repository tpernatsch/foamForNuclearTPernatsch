
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

#include "conductivityUPuO2LanningBeyer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conductivityUPuO2LanningBeyer, 0);
    addToRunTimeSelectionTable
    (
        conductivityModel, 
        conductivityUPuO2LanningBeyer, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityUPuO2LanningBeyer::conductivityUPuO2LanningBeyer
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
    OM_(mesh_.lookupObject<volScalarField>("oxygenMetalRatio")),
    densityFrac_
    (
        dict.lookupOrDefault<scalar>
        (
            "densityFraction",
            baseDict.lookupOrDefault<scalar>("densityFraction", 0.945)
        )
    ),
    porosity_(1 - densityFrac_)
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductivityUPuO2LanningBeyer::~conductivityUPuO2LanningBeyer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::conductivityUPuO2LanningBeyer::correct
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

    //- Correct conductivity for every cell in addr
    forAll(addr, i)
    {   
        const label cellI = addr[i];

        const scalar Bui = Bu[cellI]*1e-3;

        const scalar Ti = T[cellI];

        const scalar x = 2 - OM_[cellI];

        // Computing term for correlation
        const scalar A = 0.035 + 2.8*x;
        const scalar B = (2.86-7.15*x)*1e-4;
        const scalar f = 1.87e-3*Bui;
        const scalar g = 0.038*pow(Bui,0.28);
        const scalar h = pow(1+396*exp(-6380/Ti),-1);
        const scalar C = 1.5e9;
        const scalar D = 13520;

        // Computing k95 (= thermal conductivity at 95% TD )            
        const scalar k95 = 
        pow((A+B*Ti+f+(1-0.9*exp(-0.04*Bui))*g*h),-1.0)+C*pow(Ti,-2.0)*exp(-D/Ti);

        //- Computing k
        const scalar nominalValue = k95 * 1.0789 * (1-porosity_) / (1+0.5*porosity_);

        sf[cellI] = nominalValue;
    }
}

// ************************************************************************* //
