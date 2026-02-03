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

#include "heatCapacityMatproZircaloy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatCapacityMatproZircaloy, 0);
    addToRunTimeSelectionTable
    (
        heatCapacityModel, 
        heatCapacityMatproZircaloy, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatCapacityMatproZircaloy::heatCapacityMatproZircaloy
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    heatCapacityModel(mesh, dict, baseDict, defaultModel)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatCapacityMatproZircaloy::~heatCapacityMatproZircaloy()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::heatCapacityMatproZircaloy::correct
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

        // NOTE: lower bound from literature (300) extended to 273 for allowing 
        // simulation below 290 K
        if (Ti < 272.9 || Ti > 2099)
        {
            WarningInFunction
                << "Supplied temperature, " << Ti << ", is outside of range 273 < T < 2099 K";
        }

        //- Initialize values
        scalar nominalValue(0.0);
        scalar Tup(0);
        scalar Tlow(0);
        scalar CpUp(0);
        scalar CpLow(0);

        //- For each Temperature range two interpolation boundaries for Cp
        //- are defined:
        if (Ti <= 400)
        {
            Tlow = 300;
            Tup  = 400;
            CpLow = 281;
            CpUp  = 302;
        }
        else if (Ti <= 640)
        {
            Tlow = 400;
            Tup  = 640;
            CpLow = 302;
            CpUp  = 331;
        }
        else if (Ti <= 1090)
        {
            Tlow = 640;
            Tup  = 1090;
            CpLow = 331;
            CpUp  = 375;
        }
        else if (Ti <= 1093)
        {
            Tlow = 1090;
            Tup  = 1093;
            CpLow = 375;
            CpUp  = 502;
        }
        else if (Ti <= 1113)
        {
            Tlow = 1090;
            Tup  = 1113;
            CpLow = 502;
            CpUp  = 590;
        }
        else if (Ti <= 1133)
        {
            Tlow = 1113;
            Tup  = 1133;
            CpLow = 590;
            CpUp  = 615;
        }
        else if (Ti <= 1153)
        {
            Tlow = 1133;
            Tup  = 1153;
            CpLow = 615;
            CpUp  = 719;
        }
        else if (Ti <= 1173)
        {
            Tlow = 1153;
            Tup  = 1173;
            CpLow = 719;
            CpUp  = 816;
        }
        else if (Ti <= 1193)
        {
            Tlow = 1173;
            Tup  = 1193;
            CpLow = 816;
            CpUp  = 770;
        }
        else if (Ti <= 1213)
        {
            Tlow = 1193;
            Tup  = 1213;
            CpLow = 770;
            CpUp  = 619;
        }
        else if (Ti <= 1233)
        {
            Tlow = 1213;
            Tup  = 1233;
            CpLow = 619;
            CpUp  = 469;
        }
        else if (Ti <= 1248)
        {
            Tlow = 1233;
            Tup  = 1248;
            CpLow = 469;
            CpUp  = 356;
        }
        else
        {
            Tlow = 1248;
            Tup  = 2099;
            CpLow = 356;
            CpUp  = 356;            
        }

        //- Linear interpolation to find the Cp value (J/KgK)
        nominalValue = CpLow + (CpUp - CpLow) * (Ti - Tlow) / (Tup - Tlow);

        sf[cellI] = nominalValue;

    }
}




    
 
// ************************************************************************* //
