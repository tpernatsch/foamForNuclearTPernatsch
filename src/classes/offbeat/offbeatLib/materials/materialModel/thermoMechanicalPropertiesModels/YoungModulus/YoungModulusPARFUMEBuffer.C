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

#include "YoungModulusPARFUMEBuffer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(YoungModulusPARFUMEBuffer, 0);
    addToRunTimeSelectionTable
    (
        YoungModulusModel,
        YoungModulusPARFUMEBuffer,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::YoungModulusPARFUMEBuffer::YoungModulusPARFUMEBuffer
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel
)
:
    YoungModulusModel(mesh, dict, defaultModel),
    fastFluenceName_(dict.lookupOrDefault<word>("fastFluenceName", "fastFluence")),
    phi_(nullptr),
    densityName_(dict.lookupOrDefault<word>("densityName", "rho")),
    rho_(nullptr),
    par1(25.5),
    par2(0.384),
    par3(0.324e-3),
    par4(0.23),
    par5(1.5e-4),
    par6(20.0),
    perturb(1)
{
    if(dict.found("YoungModulus"))
    {
        const dictionary& YoungModulusDict = dict.subDict("YoungModulus");

        par1 = YoungModulusDict.lookupOrDefault<scalar>("par1", 25.5);
        par2 = YoungModulusDict.lookupOrDefault<scalar>("par2", 0.384);
        par3 = YoungModulusDict.lookupOrDefault<scalar>("par3", 0.324e-3);
        par4 = YoungModulusDict.lookupOrDefault<scalar>("par4", 0.23);
        par5 = YoungModulusDict.lookupOrDefault<scalar>("par5", 1.5e-4);
        par6 = YoungModulusDict.lookupOrDefault<scalar>("par6", 20.0);

        perturb = YoungModulusDict.lookupOrDefault<scalar>("perturb", 1.0);
    }

}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::YoungModulusPARFUMEBuffer::~YoungModulusPARFUMEBuffer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::YoungModulusPARFUMEBuffer::correct
(
    scalarField& sf,
    const scalarField& T,
    const labelList& addr
)
{
    if(phi_ == nullptr)
    {
        phi_ = &mesh_.lookupObject<volScalarField>(fastFluenceName_);
    }
    if(rho_ == nullptr)
    {
        rho_ = &mesh_.lookupObject<volScalarField>(densityName_);
    }

    forAll(addr, i)
    {
        const label cellI = addr[i];

        const scalar Ti = T[cellI] - 273.15;  //-temperature in Celsius
        scalar phi = phi_->internalField()[cellI] * 1e4 / 1e25; //-fast neutron fluence in 1e25n/m2;
        const scalar rho = rho_->internalField()[cellI];

        //- The functional relationship is calculated with fast fluence value of 3.96x10^25 n/m2
        //- for fast fluences (E>0.18 MeV) greater than 3.96x10^25 n/m22.
        phi = std::min(phi, 3.96);

        scalar nominalValue(0.0);

        nominalValue = (par1*(par2+par3*rho)*(1.0+par4*phi)*(1.0+par5*(Ti-par6)))*1e9;

        sf[cellI] = nominalValue*perturb;
    }
}




// ************************************************************************* //
