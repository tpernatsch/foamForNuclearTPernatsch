
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

#include "conductivityUPuO2Kato.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conductivityUPuO2Kato, 0);
    addToRunTimeSelectionTable
    (
        conductivityModel, 
        conductivityUPuO2Kato, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityUPuO2Kato::conductivityUPuO2Kato
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    conductivityModel(mesh, dict, defaultModel),
    densityFrac_(readScalar(dict.lookup("densityFraction"))),
    p_(nullptr),
    Am_(readScalar(dict.lookup("ratioAmMetal"))),
    Np_(readScalar(dict.lookup("ratioNpMetal"))),
    x_(2-dict.lookupOrDefault<scalar>("oxygenMetalRatio",2)),
    perturb(1)  
{
    if(dict.found("conductivity"))
    {
        const dictionary& conductivityDict = dict.subDict("conductivity");
        perturb = conductivityDict.lookupOrDefault<scalar>("perturb", 1.0);
    }

}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductivityUPuO2Kato::~conductivityUPuO2Kato()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::conductivityUPuO2Kato::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{

    scalar PiLim(1);
    scalar Pi(1);

    if (p_ == nullptr)
    {
        if ( mesh_.foundObject<volScalarField>("porosity") )
        {
            p_ = &mesh_.lookupObject<volScalarField>("porosity");
        }
    }

    //- Correct conductivity for every cell in addr
    forAll(addr, i)
    {   
        const label cellI = addr[i];

        // Temperature
        const scalar Ti = T[cellI];

        // Porosity
        if (p_ != nullptr) 
        {
            Pi = p_->internalField()[cellI];
        }
        
        // Maximum limit for porosity in this correlation, otherwise k->0 !
        if( !mesh_.foundObject<fvMesh>("referenceMesh") )
        {
            PiLim = min(Pi, 0.95); 
        }

        // Compute correction factor f(p) - Noascone
        const scalar fp = (1-PiLim)/(1+0.5*PiLim);

        // Compute denominator
        const scalar den = 1.595e-2 + 2.713*x_ + 3.583e-1*Am_ 
                         + 6.317e-2 *Np_ + (-2.625*x_ + 2.493)*1e-4*Ti;

        // Compute K(T)
        const scalar k = (1/(den)+1.541e11/pow(Ti,2.5)*exp(-1.522e4/Ti));

        // Correct k - Noascone method
        const scalar nominalValue = fp * k;

        // Correct k - Barani method
        // const scalar kappaHe(0.69);
        // const scalar nominalValue = k * (kappaHe+2*k-2*Pi*(k-kappaHe))
        //                               / (kappaHe+2*k + Pi*(k-kappaHe));
        
        // Assign k to scalar field
        sf[cellI] = nominalValue*perturb;
    }
}

// ************************************************************************* //
