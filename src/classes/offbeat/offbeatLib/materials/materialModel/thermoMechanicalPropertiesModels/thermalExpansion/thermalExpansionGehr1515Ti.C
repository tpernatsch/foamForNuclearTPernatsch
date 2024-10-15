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

#include "thermalExpansionGehr1515Ti.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalExpansionGehr1515Ti, 0);
    addToRunTimeSelectionTable
    (
        thermalExpansionModel, 
        thermalExpansionGehr1515Ti, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalExpansionGehr1515Ti::thermalExpansionGehr1515Ti
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    thermalExpansionModel(mesh, dict, defaultModel),
    par1(-3.101e-4),
    par2(1.545e-5),
    par3(2.75e-9),
    perturb(1.0)
{
    if(dict.found("thermalExpansion"))
    {
        const dictionary& thermalExpansionDict = dict.subDict("thermalExpansion");

        par1 = thermalExpansionDict.lookupOrDefault<scalar>("par1", -3.101e-4);
        par2 = thermalExpansionDict.lookupOrDefault<scalar>("par2", 1.545e-5);
        par3 = thermalExpansionDict.lookupOrDefault<scalar>("par3", 2.75e-9);

        perturb = thermalExpansionDict.lookupOrDefault<scalar>("perturb", 1.0);    
    }

}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalExpansionGehr1515Ti::~thermalExpansionGehr1515Ti()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::thermalExpansionGehr1515Ti::correct
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

        // Temperature must be given in °C
        const scalar Ti = T[cellI] - 273.15;

        Foam::symmTensor nominalValue = (par1 + par2 * Ti + par3*pow(Ti,2))*I;

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

        if( T[cellI] <= 293 )
        {
            nominalValue *= 0;
        }

        sf[cellI] = nominalValue*perturb;
    }
}

// ************************************************************************* //
