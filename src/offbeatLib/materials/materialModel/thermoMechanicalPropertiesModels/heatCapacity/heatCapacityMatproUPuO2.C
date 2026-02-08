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

#include "heatCapacityMatproUPuO2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatCapacityMatproUPuO2, 0);
    addToRunTimeSelectionTable
    (
        heatCapacityModel, 
        heatCapacityMatproUPuO2, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatCapacityMatproUPuO2::heatCapacityMatproUPuO2
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    heatCapacityModel(mesh, dict, baseDict, defaultModel),
    OM_
    (
        dict.lookupOrDefault<scalar>
        (
            "oxygenMetalRatio",
            baseDict.lookupOrDefault<scalar>("oxygenMetalRatio", 2.0)
        )
    ),
    K1_(347.4),
    K2_(3.95e-4),
    K3_(3.86e7),
    theta_(571),
    ED_(1.967e5),
    R_(Foam::constant::physicoChemical::R.value()),
    wPuO2_(readScalar(dict.subDict("isotopes").subDict("Pu").lookup("ratioOverMetal"))),
    wUO2_(readScalar(dict.subDict("isotopes").subDict("U").lookup("ratioOverMetal")))
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'heatCapacity' dictionary
    if((dict.found("heatCapacity")) || (dict.dictName() == "heatCapacity"))
    {
        const dictionary& heatCapacityDict = dict.found("heatCapacity")?
        dict.subDict("heatCapacity"):
        dict;

        OM_ = heatCapacityDict.lookupOrDefault<scalar>("OM", 2.0);
        K1_ = heatCapacityDict.lookupOrDefault<scalar>("K1", 347.4);
        K2_ = heatCapacityDict.lookupOrDefault<scalar>("K2", 3.95e4);
        K3_ = heatCapacityDict.lookupOrDefault<scalar>("K3", 3.86e7);
        theta_ = heatCapacityDict.lookupOrDefault<scalar>("theta", 571);
        ED_ = heatCapacityDict.lookupOrDefault<scalar>("ED", 1.967e5);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatCapacityMatproUPuO2::~heatCapacityMatproUPuO2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::heatCapacityMatproUPuO2::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
const
{   
    //- Heat capacity in case of pure UO2 fuel
    scalar hc_UO2;

    //- Heat capacity in case of pure PuO2 fuel
    scalar hc_PuO2;
    
    forAll(addr, i)
    {   
        const label cellI = addr[i];      

        const scalar Ti = T[cellI];

        hc_UO2 =  296.7*pow(535.285,2)*exp(535.285/Ti)
                  /pow(Ti*(exp(535.285/Ti)-1),2)
                  + 2.43e-2*Ti
                  + (OM_/2)*8.745e7*1.577e5/(R_*pow(Ti,2))*exp(-1.577e5/(R_*Ti));
        
        hc_PuO2 =  K1_*pow(theta_,2)*exp(theta_/Ti)
                    /pow(Ti*(exp(theta_/Ti)-1),2)
                    + K2_*Ti
                    + (OM_/2)*K3_*ED_/(R_*pow(Ti,2))*exp(-ED_/(R_*Ti));

        //- Weighted sum of the hc of UO2 and PuO2
        const scalar nominalValue = wUO2_*hc_UO2 + wPuO2_*hc_PuO2;

        sf[cellI] = nominalValue;
    }
}

// ************************************************************************* //
