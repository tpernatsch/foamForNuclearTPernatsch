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

#include "heatCapacityMatproUO2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatCapacityMatproUO2, 0);
    addToRunTimeSelectionTable
    (
        heatCapacityModel, 
        heatCapacityMatproUO2, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatCapacityMatproUO2::heatCapacityMatproUO2
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    heatCapacityModel(mesh, dict, defaultModel),
    OM(2.0),
    K1(296.7),
    K2(2.43e-2),
    K3(8.745e7),
    phi(535.285),
    ED(1.577e5),
    R(Foam::constant::physicoChemical::R.value()),
    par1(2.0),
    par2(2.0),
    par3(2.0),
    perturb(1.0)
{
    if(dict.found("heatCapacity"))
    {
        const dictionary& heatCapacityDict = dict.subDict("heatCapacity");

        OM = heatCapacityDict.lookupOrDefault<scalar>("OM", 2.0);
        K1 = heatCapacityDict.lookupOrDefault<scalar>("K1", 296.7);
        K2 = heatCapacityDict.lookupOrDefault<scalar>("K2", 2.43e-2);
        K3 = heatCapacityDict.lookupOrDefault<scalar>("K3", 8.745e7);
        phi = heatCapacityDict.lookupOrDefault<scalar>("phi", 535.285);
        ED = heatCapacityDict.lookupOrDefault<scalar>("ED", 1.577e5);
        par1 = heatCapacityDict.lookupOrDefault<scalar>("par1", 2.0);
        par2 = heatCapacityDict.lookupOrDefault<scalar>("par2", 2.0);
        par3 = heatCapacityDict.lookupOrDefault<scalar>("par3", 2.0);

        perturb = heatCapacityDict.lookupOrDefault<scalar>("perturb", 1.0);
        
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatCapacityMatproUO2::~heatCapacityMatproUO2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::heatCapacityMatproUO2::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
const
{  
    forAll(addr, i)
    {   
        const label cellI = addr[i];      

        const scalar Ti = T[cellI];
        
        const scalar nominalValue =  K1*pow(phi, par1)*exp(phi/Ti)/pow(Ti*(exp(phi/Ti)-1), par2)
         + K2*Ti
         + (OM/2)*K3*ED/(R*pow(Ti, par3))*exp(-ED/(R*Ti));

        sf[cellI] = nominalValue*perturb;
    }
}



    
 
// ************************************************************************* //
