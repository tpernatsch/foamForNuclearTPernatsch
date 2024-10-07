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

#include "thermalExpansionLemehovUPuO2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalExpansionLemehovUPuO2, 0);
    addToRunTimeSelectionTable
    (
        thermalExpansionModel, 
        thermalExpansionLemehovUPuO2, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalExpansionLemehovUPuO2::thermalExpansionLemehovUPuO2
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    thermalExpansionModel(mesh, dict, defaultModel),
    b0_(-0.3080),
    b1_(3.4303),
    b2_(-1.9157),
    b3_(3.4636),
    by_(3.98),
    TmUO2_(3127),
    TmPuO2_(2804),
    OM_(dict.lookupOrDefault<scalar>("oxygenMetalRatio",2.0)),
    ratioPuMetal_(readScalar(dict.lookup("ratioPuMetal"))),
    perturb(1.0)
{
    if(dict.found("thermalExpansion"))
    {
        const dictionary& thermalExpansionDict = dict.subDict("thermalExpansion");

        b0_ = thermalExpansionDict.lookupOrDefault<scalar>("b0", -0.3080);
        b1_ = thermalExpansionDict.lookupOrDefault<scalar>("b1",  3.4303);
        b2_ = thermalExpansionDict.lookupOrDefault<scalar>("b2", -1.9157);
        b3_ = thermalExpansionDict.lookupOrDefault<scalar>("b3",  3.4636);
        by_ = thermalExpansionDict.lookupOrDefault<scalar>("by",  3.98);
        TmUO2_  = thermalExpansionDict.lookupOrDefault<scalar>("TmUO2",  3127);
        TmPuO2_ = thermalExpansionDict.lookupOrDefault<scalar>("TmPuO2", 2804);

        perturb = thermalExpansionDict.lookupOrDefault<scalar>("perturb", 1.0);   
        
    }
}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalExpansionLemehovUPuO2::~thermalExpansionLemehovUPuO2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::thermalExpansionLemehovUPuO2::correct
(
    symmTensorField& sf, 
    const scalarField& T, 
    const labelList& addr
)
const
{  
    forAll(addr, i)
    {   
        const label cellI = addr[i];    

        const scalar Ti = T[cellI];

        // Compute melting temperature
        const scalar Tm = TmUO2_*(1-ratioPuMetal_) + ratioPuMetal_*TmPuO2_;
        
        // Compute dL/L0
        const scalar dLOverL0 = b0_
                              + b1_*Ti/Tm
                              + b2_*pow(Ti/Tm,2.0)
                              + b3_*pow(Ti/Tm,3.0);

        // Notice that this value is multiplied by 0.01 because from
        // the reference it looks like the value is given in % units
        // (see reference for more details)
        const symmTensor nominalValue = 
            0.01*(dLOverL0 * (1+by_*(2-OM_)))*I;

        sf[cellI] = nominalValue*perturb;

    }
}   
 
// ************************************************************************* //
