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

#include "conductivityPARFUMESiC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conductivityPARFUMESiC, 0);
    addToRunTimeSelectionTable
    (
        conductivityModel,
        conductivityPARFUMESiC,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityPARFUMESiC::conductivityPARFUMESiC
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel
)
:
    conductivityModel(mesh, dict, defaultModel),
    par1(17885.0),
    par2(2.0),
    perturb(1)
{
  if(dict.found("conductivity"))
  {
      const dictionary& conductivityDict = dict.subDict("conductivity");

      par1 = conductivityDict.lookupOrDefault<scalar>("par1", 17885.0);
      par2 = conductivityDict.lookupOrDefault<scalar>("par2", 2.0);

      perturb = conductivityDict.lookupOrDefault<scalar>("perturb", 1.0);
  }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductivityPARFUMESiC::~conductivityPARFUMESiC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::conductivityPARFUMESiC::correct
(
    scalarField& sf,
    const scalarField& T,
    const labelList& addr
)
{

    forAll(addr, i)
    {
        const label cellI = addr[i];

        const scalar Ti = T[cellI]; //-temperature in Kelvin

        scalar nominalValue(0.0);

        nominalValue = par1 / Ti + par2;

        sf[cellI] = nominalValue*perturb;
    }
}





// ************************************************************************* //
