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

#include "thermalExpansionMAMOX.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalExpansionMAMOX, 0);
    addToRunTimeSelectionTable
    (
        thermalExpansionModel, 
        thermalExpansionMAMOX, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::scalar Foam::thermalExpansionMAMOX::a0( scalar cPu, scalar x) const
{
    const scalar b0 = -2.8809;
    const scalar b1 = 0.0301;
    const scalar b2 = -4.3954;
    const scalar b3 = 0.0156;
    const scalar b4 = -15.1759;
    const scalar b5 = 2.5642;
    return 1e-3*(b0 + b1*cPu + b2*x + b3*pow(cPu,2.0) + b4*pow(x,2.0) + b5*cPu*x);
}
const Foam::scalar Foam::thermalExpansionMAMOX::a1( scalar cPu, scalar x) const
{
    const scalar b0 = 9.5024;
    const scalar b1 = -0.1864;
    const scalar b2 = 15.8173;
    const scalar b3 = -0.0229;
    const scalar b4 = 7.6258;
    const scalar b5 = -7.5789;
    return 1e-6*(b0 + b1*cPu + b2*x + b3*pow(cPu,2.0) + b4*pow(x,2.0) + b5*cPu*x);
}
const Foam::scalar Foam::thermalExpansionMAMOX::a2( scalar cPu, scalar x) const
{
    const scalar b0 = 2.0894;
    const scalar b1 = 2.9483;
    const scalar b2 = -19.9227;
    const scalar b3 = -1.0355;
    const scalar b4 = 73.8931;
    const scalar b5 = 11.6442;
    return 1e-10*(b0 + b1*cPu + b2*x + b3*pow(cPu,2.0) + b4*pow(x,2.0) + b5*cPu*x);
}
const Foam::scalar Foam::thermalExpansionMAMOX::a3( scalar cPu, scalar x) const
{
    const scalar b0 = 4.4096;
    const scalar b1 = -1.4263;
    const scalar b2 = 23.5638;
    const scalar b3 = 0.0251;
    const scalar b4 = -54.751;
    const scalar b5 = -14.418;
    return 1e-13*(b0 + b1*cPu + b2*x + b3*pow(cPu,2.0) + b4*pow(x,2.0) + b5*cPu*x);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalExpansionMAMOX::thermalExpansionMAMOX
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    thermalExpansionModel(mesh, dict, baseDict, defaultModel),
    cPu_(readScalar(dict.lookup("ratioPuMetal"))),    
    OM_
    (
        dict.lookupOrDefault<scalar>
        (
            "oxygenMetalRatio",
            baseDict.lookupOrDefault<scalar>("oxygenMetalRatio", 2.0)
        )
    ),
    x_(2.0 - OM_),
    perturb(1.0)
{
}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalExpansionMAMOX::~thermalExpansionMAMOX()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::thermalExpansionMAMOX::correct
(
    symmTensorField& sf, 
    const scalarField& T, 
    const labelList& addr
)
const
{  
    // Get reference temperature
    scalar Tref(Tref_);

    scalar dL_L0_ref = a0(cPu_, x_) 
                    + a1(cPu_, x_) * Tref 
                    + a2(cPu_, x_) * pow(Tref, 2.0) 
                    + a3(cPu_, x_) * pow(Tref, 3.0);

    forAll(addr, i)
    {   
        const label cellI = addr[i];    

        const scalar Ti = T[cellI];

        scalar dL_L0 = a0(cPu_, x_) 
                    + a1(cPu_, x_) * Ti 
                    + a2(cPu_, x_) * pow(Ti, 2.0) 
                    + a3(cPu_, x_) * pow(Ti, 3.0);
    
        
        // Compute thermal expansion tensor
        symmTensor nominalValue = (dL_L0-dL_L0_ref)*I;

        // NOTE: avoids instabilities when trying to simulate a material at 
        // constant temperature equal to Tref
        if
        (
            mag(nominalValue.xx()) < 1e-7 &&
            mag(nominalValue.yy()) < 1e-7 &&
            mag(nominalValue.zz()) < 1e-7
        )
        {
            nominalValue *= 0;
        }

        sf[cellI] = nominalValue*perturb; 

    }
}   
 
// ************************************************************************* //
