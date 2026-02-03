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

#include "Inconel600.H"
#include "addToRunTimeSelectionTable.H"
#include "userParameters.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Inconel600, 0);
    addToRunTimeSelectionTable(materialModel, Inconel600, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Inconel600::Inconel600
(
    const fvMesh& mesh, 
    const dictionary& materialModelDict,
    const labelList& addr
)
:
    materialModel(mesh, materialModelDict, addr),
    rho_("rho", dimensionSet(1, -3, 0, 0, 0), materialModelDict)
{    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Inconel600::~Inconel600()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::Inconel600::rho
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{
    forAll(addr, i)
    {
        sf[addr[i]] = rho_.value();
    }
}

void Foam::Inconel600::Cp
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{
    forAll(addr, i)
    {
        const label cellI = addr[i];
        sf[cellI] = -2.80E-07*T[cellI]*T[cellI]*T[cellI] + 
                    6.65E-04*T[cellI]*T[cellI] - 
                    2.81E-01*T[cellI] +
                    494.663;
    }
}

void Foam::Inconel600::k
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{
    forAll(addr, i)
    {
        const label cellI = addr[i];
        sf[cellI] = 1.66E-02*T[cellI] + 
                    9.48;
    }
}

void Foam::Inconel600::emissivity
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{
    forAll(addr, i)
    {
        const label cellI = addr[i];
        sf[cellI] = 1.98E-04*T[cellI] +
                    5.73E-01;
    }
}


void Foam::Inconel600::E
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{
    forAll(addr, i)
    {
        const label cellI = addr[i];
        sf[cellI] = -7.07E+07*T[cellI] +
                    2.39E+11;
    }
}

void Foam::Inconel600::nu
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{
    forAll(addr, i)
    {
        const label cellI = addr[i];
        sf[addr[i]] = -3.66E-13*T[cellI]*T[cellI]*T[cellI]*T[cellI] +
                    1.18E-09*T[cellI]*T[cellI]*T[cellI] - 
                    1.21E-06*T[cellI]*T[cellI] + 
                    4.32E-04*T[cellI] +
                    2.74E-01;
    }
}


void Foam::Inconel600::alphaT
(
    symmTensorField& sf, 
    const scalarField& T, 
    const labelList& addr
) const
{
    forAll(addr, i)
    {
        const label cellI = addr[i];
        const scalar Ti = T[cellI];
        sf[cellI] = (1.19E-05- + 3.88E-09*Ti)*(Ti - 293.15)*I;
       // sf[cellI] = 0.0060373-0.000011175*Ti+0.000000013557*Ti*Ti-0.0000000000046688*Ti*Ti*Ti+6.9032E-016*Ti*Ti*Ti*Ti;
    }
}

// ************************************************************************* //
