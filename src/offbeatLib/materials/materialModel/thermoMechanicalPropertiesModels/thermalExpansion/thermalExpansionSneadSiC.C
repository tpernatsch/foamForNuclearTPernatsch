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

#include "thermalExpansionSneadSiC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalExpansionSneadSiC, 0);
    addToRunTimeSelectionTable
    (
        thermalExpansionModel,
        thermalExpansionSneadSiC,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalExpansionSneadSiC::thermalExpansionSneadSiC
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    thermalExpansionModel(mesh, dict, baseDict, defaultModel),
    par1(-1.8276),
    par2(0.0178),
    par3(-1.5544e-5),
    par4(4.5246e-9),
    par5(5e-6),
    Treference(298.15)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'thermalExpansion' dictionary
    if((dict.found("thermalExpansion")) || (dict.dictName() == "thermalExpansion"))
    {
        const dictionary& thermalExpansionDict = dict.found("thermalExpansion")?
        dict.subDict("thermalExpansion"):
        dict;

        par1 = thermalExpansionDict.lookupOrDefault<scalar>("par1", -1.8276);
        par2 = thermalExpansionDict.lookupOrDefault<scalar>("par2", 0.0178);
        par3 = thermalExpansionDict.lookupOrDefault<scalar>("par3", -1.5544e-5);
        par4 = thermalExpansionDict.lookupOrDefault<scalar>("par4", 4.5246e-9);
        par5 = thermalExpansionDict.lookupOrDefault<scalar>("par5", 5e-6);
        Treference = thermalExpansionDict.lookupOrDefault<scalar>("Treference", 298.15);
    }
}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalExpansionSneadSiC::~thermalExpansionSneadSiC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::thermalExpansionSneadSiC::correct
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

        // Temperature
        const scalar Ti = T[cellI];
        //To distinguish the stress free temperature(which is the reference temperature
        //in OFFBEAT) and the reference temperature (which is used to determine the mean
        //thermal expansion coef), here we use Tsf as stree free temperature.
        const scalar Tsf = Tref_;

        //-mean thermal expansion coefficient
        scalar alpha;

        //-mean thermal expansion coefficient at stress free temperature
        scalar alpha_sf;

        if (Ti <1273.15)
        {
          alpha = 1e-6*(par1+par2*Ti+par3*pow(Ti,2)+par4*pow(Ti,3));
        }
        else
        {
          alpha = par5;
        }

        alpha_sf = 1e-6*(par1+par2*Tsf+par3*pow(Tsf,2)+par4*pow(Tsf,3));

        // When Treference = Tsf, the equation is the same as the instantaneous one
        symmTensor nominalValue = (alpha*(Ti-Treference)-alpha_sf*(Tsf-Treference))
                                   / (1+alpha_sf*(Tsf-Treference))*I;


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

        sf[cellI] = nominalValue;
    }
}

// ************************************************************************* //
