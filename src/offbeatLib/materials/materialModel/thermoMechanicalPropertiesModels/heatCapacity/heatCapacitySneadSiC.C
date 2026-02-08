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

#include "heatCapacitySneadSiC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatCapacitySneadSiC, 0);
    addToRunTimeSelectionTable
    (
        heatCapacityModel,
        heatCapacitySneadSiC,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatCapacitySneadSiC::heatCapacitySneadSiC
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    heatCapacityModel(mesh, dict, baseDict, defaultModel),
    par1(925.65),
    par2(0.3772),
    par3(-7.9259e-5),
    par4(-3.1946e7)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'heatCapacity' dictionary
    if((dict.found("heatCapacity")) || (dict.dictName() == "heatCapacity"))
    {
        const dictionary& heatCapacityDict = dict.found("heatCapacity")?
        dict.subDict("heatCapacity"):
        dict;

        par1 = heatCapacityDict.lookupOrDefault<scalar>("par1", 925.65);
        par2 = heatCapacityDict.lookupOrDefault<scalar>("par2", 0.3772);
        par3 = heatCapacityDict.lookupOrDefault<scalar>("par3", -7.9259e-5);
        par4 = heatCapacityDict.lookupOrDefault<scalar>("par4", -3.19446e7);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatCapacitySneadSiC::~heatCapacitySneadSiC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::heatCapacitySneadSiC::correct
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

        if (Ti < 200.0|| Ti > 2400.0)
        {
            WarningInFunction
                << "Foam::heatCapacitySneadSiC::correct() :\n"
                << "Supplied temperature, " << Ti
                << ", is outside of range 200 < T < 2400 K " << endl;
        }

        const scalar nominalValue = par1 + par2 * Ti + par3 * pow(Ti,2) + par4 * pow(Ti,-2);

        sf[cellI] = nominalValue;

    }
}






// ************************************************************************* //
