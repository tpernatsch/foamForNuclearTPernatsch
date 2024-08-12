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

#include "YoungModulusPARFUMEPyC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(YoungModulusPARFUMEPyC, 0);
    addToRunTimeSelectionTable
    (
        YoungModulusModel,
        YoungModulusPARFUMEPyC,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::YoungModulusPARFUMEPyC::YoungModulusPARFUMEPyC
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
    BAF_(dict.lookupOrDefault("asFabricatedAnisotropy", 1.0)),
    Lc_(dict.lookupOrDefault("crystalliteDiameter", 30.0)),
    par1(25.5),
    par2(0.384),
    par3(0.324e-3),
    par4(1.463),
    par5(-0.463),
    par6(2.985),
    par7(-0.0662),
    par8(0.23),
    par9(0.00015),
    par10(20.0),
    par11(0.481),
    par12(0.519),
    perturb(1.0)
{

  if (BAF_ < 1.0)
  {
    WarningInFunction
        << "The Initial Bacon Anisotropy Factor (BAF) must be greater than 1.0";
  }
    if(dict.found("YoungModulus"))
    {
        const dictionary& YoungModulusDict = dict.subDict("YoungModulus");

        par1 = YoungModulusDict.lookupOrDefault<scalar>("par1", 25.5);
        par2 = YoungModulusDict.lookupOrDefault<scalar>("par2", 0.384);
        par3 = YoungModulusDict.lookupOrDefault<scalar>("par3", 0.324e-3);
        par4 = YoungModulusDict.lookupOrDefault<scalar>("par4", 1.463);
        par5 = YoungModulusDict.lookupOrDefault<scalar>("par5", -0.463);
        par6 = YoungModulusDict.lookupOrDefault<scalar>("par6", 2.985);
        par7 = YoungModulusDict.lookupOrDefault<scalar>("par7", -0.0662);
        par8 = YoungModulusDict.lookupOrDefault<scalar>("par8", 0.23);
        par9 = YoungModulusDict.lookupOrDefault<scalar>("par9", 0.00015);
        par10 = YoungModulusDict.lookupOrDefault<scalar>("par10", 20.0);
        par11 = YoungModulusDict.lookupOrDefault<scalar>("par11", 0.481);
        par12 = YoungModulusDict.lookupOrDefault<scalar>("par12", 0.519);

        perturb = YoungModulusDict.lookupOrDefault<scalar>("perturb", 1.0);
    }

}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::YoungModulusPARFUMEPyC::~YoungModulusPARFUMEPyC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::YoungModulusPARFUMEPyC::correct
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
        scalar phi = phi_->internalField()[cellI]*1e4 / 1e25; //-fast neutron fluence in 10e25n/m2;
        const scalar rho = rho_->internalField()[cellI];

        //- The functional relationship is calculated with fast fluence value of 3.96x102525 n/m2
        //- for fast fluences (E>0.18 MeV) greater than 3.96x102525 n/m22.
        phi = std::min(phi, 3.96);

        if (rho < 1800.0 || rho > 2000.0)
        {
          WarningInFunction
              << "Supplied density, " << rho << ", is out of range 1800 < rho < 2000 kg/m3";
        }
        //- radial component
        const scalar Er = par1*(par2+par3*rho)*(par4+par5*BAF_)*(par6+par7*Lc_)
                          *(1.0+par8*phi)*(1.0+par9*(Ti-par10));
        //- tangential component
        const scalar Et = par1*(par2+par3*rho)*(par11+par12*BAF_)*(par6+par7*Lc_)
                          *(1.0+par8*phi)*(1.0+par9*(Ti-par10));

        scalar nominalValue(0.0);

        //- TO DO: Normally, the The PyC material should be transverse isotropic, with the
        //- different value for radial and tangential Young's Modulus. However, in the case
        //- the BAF = 1.0 and Lc_ = 30 (which usually are the default values for TRISO simulation), Er = Et.
        //- Hence, here we use isotropic value instead temporarily.

        nominalValue = ((Er + 2.0*Et) / 3.0)*1e9;

        sf[cellI] = nominalValue*perturb;
    }
}




// ************************************************************************* //
