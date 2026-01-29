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

#include "conductivityParfumeBuffer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conductivityParfumeBuffer, 0);
    addToRunTimeSelectionTable
    (
        conductivityModel,
        conductivityParfumeBuffer,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityParfumeBuffer::conductivityParfumeBuffer
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    conductivityModel(mesh, dict, baseDict, defaultModel),
    densityName_(dict.lookupOrDefault<word>("densityName", "rho" )),
    rho_(nullptr),
    rho_init_(dict.lookupOrDefault("initialDensity", 1000.0)),
    rho_theo_
    (
        dict.lookupOrDefault<scalar>
        (
            "theoreticalDensity",
            baseDict.lookupOrDefault<scalar>("theoreticalDensity", 2250.0)
        )
    ),
    k_init_(dict.lookupOrDefault("initialConductivity", 0.5)),
    k_theo_(dict.lookupOrDefault("theoreticalConductivity", 4.0))
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductivityParfumeBuffer::~conductivityParfumeBuffer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::conductivityParfumeBuffer::correct
(
    scalarField& sf,
    const scalarField& T,
    const labelList& addr
)
{

    forAll(addr, i)
    {
        const label cellI = addr[i];

        if (rho_ == nullptr)
        {
            rho_ = &mesh_.lookupObject<volScalarField>(densityName_);
        }

        const volScalarField& rhoI = *rho_;

        const scalar nominalValue =  k_init_*k_theo_*rho_theo_*(rho_theo_-rho_init_) /
                                     (k_theo_*rho_theo_*(rho_theo_-rhoI[cellI])+k_init_*rhoI[cellI]*(rhoI[cellI]-rho_init_));

        sf[cellI] = nominalValue;
    }
}





// ************************************************************************* //
