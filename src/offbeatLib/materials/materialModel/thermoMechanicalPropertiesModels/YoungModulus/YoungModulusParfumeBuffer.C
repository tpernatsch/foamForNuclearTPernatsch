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

#include "YoungModulusParfumeBuffer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(YoungModulusParfumeBuffer, 0);
    addToRunTimeSelectionTable
    (
        YoungModulusModel,
        YoungModulusParfumeBuffer,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::YoungModulusParfumeBuffer::YoungModulusParfumeBuffer
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    YoungModulusModel(mesh, dict, baseDict, defaultModel),
    fastFluenceName_(dict.lookupOrDefault<word>("fastFluenceName", "fastFluence")),
    phi_(nullptr),
    densityName_(dict.lookupOrDefault<word>("densityName", "rho")),
    rho_(nullptr),
    par1(25.5),
    par2(0.384),
    par3(0.324e-3),
    par4(0.23),
    par5(1.5e-4),
    par6(20.0)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'YoungModulus' dictionary
    if((dict.found("YoungModulus")) || (dict.dictName() == "YoungModulus"))
    {
        const dictionary& YoungModulusDict = dict.found("YoungModulus")?
        dict.subDict("YoungModulus"):
        dict;

        par1 = YoungModulusDict.lookupOrDefault<scalar>("par1", 25.5);
        par2 = YoungModulusDict.lookupOrDefault<scalar>("par2", 0.384);
        par3 = YoungModulusDict.lookupOrDefault<scalar>("par3", 0.324e-3);
        par4 = YoungModulusDict.lookupOrDefault<scalar>("par4", 0.23);
        par5 = YoungModulusDict.lookupOrDefault<scalar>("par5", 1.5e-4);
        par6 = YoungModulusDict.lookupOrDefault<scalar>("par6", 20.0);
    }

}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::YoungModulusParfumeBuffer::~YoungModulusParfumeBuffer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::YoungModulusParfumeBuffer::correct
(
    scalarField& sf,
    const scalarField& T,
    const labelList& addr
)
{
    if(phi_ == nullptr)
    {
        if(mesh_.foundObject<volScalarField>(fastFluenceName_))
        {
            phi_ = &mesh_.lookupObject<volScalarField>(fastFluenceName_);
        }
    }

    if(rho_ == nullptr)
    {
        rho_ = &mesh_.lookupObject<volScalarField>(densityName_);
    }

    forAll(addr, i)
    {
        const label cellI = addr[i];

        const scalar Ti = T[cellI] - 273.15;  //-temperature in Celsius
        
        // Get neutron fluence and convert to 1e25n/m2
        scalar phi = 0.0;
        if(phi_ != nullptr)
        {
            phi  = phi_->internalField()[cellI] * 1e4 / 1e25;
        }

        const scalar rho = rho_->internalField()[cellI];

        //- The functional relationship is calculated with fast fluence value of 3.96x10^25 n/m2
        //- for fast fluences (E>0.18 MeV) greater than 3.96x10^25 n/m22.
        phi = std::min(phi, 3.96);

        scalar nominalValue(0.0);

        nominalValue = (par1*(par2+par3*rho)*(1.0+par4*phi)*(1.0+par5*(Ti-par6)))*1e9;

        sf[cellI] = nominalValue;
    }
}




// ************************************************************************* //
