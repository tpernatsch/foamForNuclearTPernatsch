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

#include "conductivityUO2IFA601.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conductivityUO2IFA601, 0);
    addToRunTimeSelectionTable
    (
        conductivityModel, 
        conductivityUO2IFA601, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityUO2IFA601::conductivityUO2IFA601
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    conductivityModel(mesh, dict, defaultModel),       
    burnupName_(dict.lookupOrDefault<word>("burnupName", "Bu" )),
    Bu_(nullptr),
    par1(40),
    par2(5),
    par3(0.24),
    par4(0.457),
    par5(0.00433),
    par6(0.11),
    par7(0.0000132),
    par8(0.00188),
    perturb(1)
{
    if(dict.found("conductivity"))
    {
        const dictionary& conductivityDict = dict.subDict("conductivity");

        par1 = conductivityDict.lookupOrDefault<scalar>("par1", 40);
        par2 = conductivityDict.lookupOrDefault<scalar>("par2", 5);
        par3 = conductivityDict.lookupOrDefault<scalar>("par3", 0.24);
        par4 = conductivityDict.lookupOrDefault<scalar>("par4", 0.457);
        par5 = conductivityDict.lookupOrDefault<scalar>("par5", 0.00433);
        par6 = conductivityDict.lookupOrDefault<scalar>("par6", 0.11);
        par7 = conductivityDict.lookupOrDefault<scalar>("par7", 0.0000132);
        par8 = conductivityDict.lookupOrDefault<scalar>("par8", 0.00188);

        perturb = conductivityDict.lookupOrDefault<scalar>("perturb", 1.0);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductivityUO2IFA601::~conductivityUO2IFA601()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::conductivityUO2IFA601::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{   
    if(Bu_ == nullptr)
    {        
        Bu_ = &mesh_.lookupObject<volScalarField>(burnupName_);
    }
 
    const scalarField& Bu = Bu_->internalField();

    forAll(addr, i)
    {   
        const label cellI = addr[i];              

        // Convert burnup to MWd/kgU
        const scalar Bui = Bu[cellI]/1000;

        const scalar Ti = T[cellI];
            
        double term1 = par1 + par2*Bui;
        double term2 = par3; 
        if(Bui > 50)
        {
            if(Bui <= 80)
            {
                term2 = par4 - par5*Bui;
            }
            else
            {
                term2 = par6;
            }
        }

        const scalar nominalValue = 1/(term1 + term2*Ti) + par7*exp(par8*Ti);

        sf[cellI] = nominalValue*perturb*1000;
    }
}

// ************************************************************************* //
