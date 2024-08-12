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

#include "conductivityRelapZy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conductivityRelapZy, 0);
    addToRunTimeSelectionTable
    (
        conductivityModel, 
        conductivityRelapZy, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityRelapZy::conductivityRelapZy
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    conductivityModel(mesh, dict, defaultModel),     
    par1(7.51),
    par2(2.09e-2),
    par3(1.45e-5),
    par4(7.67e-9),
    perturb(1)
{
    if(dict.found("conductivity"))
    {
        const dictionary& conductivityDict = dict.subDict("conductivity");

        par1 = conductivityDict.lookupOrDefault<scalar>("par1", 7.51);
        par2 = conductivityDict.lookupOrDefault<scalar>("par2", 2.09e-2);
        par3 = conductivityDict.lookupOrDefault<scalar>("par3", 1.45e-5);
        par4 = conductivityDict.lookupOrDefault<scalar>("par4", 7.67e-9);

        perturb = conductivityDict.lookupOrDefault<scalar>("perturb", 1.0);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductivityRelapZy::~conductivityRelapZy()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::conductivityRelapZy::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{    

    forAll(addr, i)
    {   
        const label cellI = addr[i];              

        const scalar Ti = T[cellI];

        //- The validity range is changed by one degree on each side to allow 
        //- for rounding errors
        if (Ti < 290 || Ti > 1801)
        {
            WarningInFunction
                << "Supplied temperature, " << Ti << ", out of range 290 < T < 1800 K";
        }
       
        const scalar nominalValue =  par1 + par2*Ti - par3*Ti*Ti + par4*Ti*Ti*Ti;

        sf[cellI] = nominalValue*perturb;
    }
}


    


// ************************************************************************* //
