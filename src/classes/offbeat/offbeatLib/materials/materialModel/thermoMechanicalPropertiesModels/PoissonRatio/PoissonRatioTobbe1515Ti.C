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

#include "PoissonRatioTobbe1515Ti.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PoissonRatioTobbe1515Ti, 0);
    addToRunTimeSelectionTable
    (
        PoissonRatioModel, 
        PoissonRatioTobbe1515Ti, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PoissonRatioTobbe1515Ti::PoissonRatioTobbe1515Ti
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    PoissonRatioModel(mesh, dict, defaultModel),
    par1(0.277),
    par2(6e-5),
    perturb(1.0)
{
    if(dict.found("PoissonRatio"))
    {
        const dictionary& PoissonRatioDict = dict.subDict("PoissonRatio");

        par1 = PoissonRatioDict.lookupOrDefault<scalar>("par1", 0.277);        
        par2 = PoissonRatioDict.lookupOrDefault<scalar>("par2", 6e-5);

        perturb = PoissonRatioDict.lookupOrDefault<scalar>("perturb", 1.0);
    }

}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PoissonRatioTobbe1515Ti::~PoissonRatioTobbe1515Ti()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::PoissonRatioTobbe1515Ti::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{     
    forAll(addr, i)
    {   
        const label cellI = addr[i];        

        const scalar Ti = T[cellI] - 273.15;  

        const scalar nominalValue = par1 + par2*Ti;

        sf[cellI] = nominalValue*perturb;
    }
}

        
        


// ************************************************************************* //
