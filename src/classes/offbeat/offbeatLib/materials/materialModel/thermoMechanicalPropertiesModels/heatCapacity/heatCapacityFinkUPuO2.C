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

#include "heatCapacityFinkUPuO2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatCapacityFinkUPuO2, 0);
    addToRunTimeSelectionTable
    (
        heatCapacityModel, 
        heatCapacityFinkUPuO2, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatCapacityFinkUPuO2::heatCapacityFinkUPuO2
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    heatCapacityModel(mesh, dict, defaultModel),
    C1_(322.49),
    C2_(1.4679e-2),
    C3_(0),
    theta_(587.41),
    Ea_(18531.7),
    perturb(1.0)
{
    if(dict.found("heatCapacity"))
    {
        const dictionary& heatCapacityDict = dict.subDict("heatCapacity");

        C1_ = heatCapacityDict.lookupOrDefault<scalar>("C1", 322.49);
        C2_ = heatCapacityDict.lookupOrDefault<scalar>("C2", 1.4679e-2);
        C3_ = heatCapacityDict.lookupOrDefault<scalar>("C3", 0);
        theta_ = heatCapacityDict.lookupOrDefault<scalar>("theta", 587.41);
        Ea_ = heatCapacityDict.lookupOrDefault<scalar>("Ea", 18531.7);

        perturb = heatCapacityDict.lookupOrDefault<scalar>("perturb", 1.0);
        
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatCapacityFinkUPuO2::~heatCapacityFinkUPuO2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::heatCapacityFinkUPuO2::correct
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

        const scalar nominalValue = 
            C1_*pow(theta_/Ti,2.0)*exp(theta_/Ti)/pow(exp(theta_/Ti)-1.0,2.0)
          + 2*C2_*Ti 
          + C3_*Ea_*exp(-Ea_/Ti)/pow(Ti,2.0); 

        sf[cellI] = nominalValue*perturb;
    }
}

// ************************************************************************* //
