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

#include "UPuO2.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "gapGasModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(UPuO2, 0);
    addToRunTimeSelectionTable(fuelMaterial, UPuO2, dictionary);

    const char* UPuO2::group_ = ("materials::UPuO2");

    defineParameter(UPuO2, F_rho_, "F_rho", (dimless), 1.0);
    defineParameter(UPuO2, F_Cp_, "F_Cp", (dimless), 1.0);
    defineParameter(UPuO2, F_k_, "F_k", (dimless), 1.0);
    defineParameter(UPuO2, F_emissivity_, "F_emissivity", (dimless), 1.0);
    defineParameter(UPuO2, F_E_, "F_E", (dimless), 1.0);
    defineParameter(UPuO2, F_nu_, "F_nu", (dimless), 1.0);
    defineParameter(UPuO2, F_alphaT_, "F_alphaT", (dimless), 1.0);

    defineParameter(UPuO2, delta_rho_, "delta_rho", (dimMass/dimVolume), 0.0);
    defineParameter(UPuO2, delta_Cp_, "delta_Cp", (0, 2, -2, -1, 0), 0.0);
    defineParameter(UPuO2, delta_k_, "delta_k", (1, 1, -3, -1, 0), 0.0);
    defineParameter(UPuO2, delta_emissivity_, "delta_emissivity", (dimless), 0.0);
    defineParameter(UPuO2, delta_E_, "delta_E", (dimPressure), 0.0);
    defineParameter(UPuO2, delta_nu_, "delta_nu", (dimless), 0.0);
    defineParameter(UPuO2, delta_alphaT_, "delta_alphaT", (dimless), 0.0);  
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UPuO2::UPuO2
(
    const fvMesh& mesh, 
    const dictionary& materialModelDict,
    const labelList& addr
)
:
    fuelMaterial(mesh, materialModelDict, addr),
    PuMetalRatio_(readScalar(isotopesDict_.subDict("Pu").lookup("ratioOverMetal"))),
    UMetalRatio_(readScalar(isotopesDict_.subDict("U").lookup("ratioOverMetal")))
{
    density_ =
        densityModel::New(mesh, materialModelDict, "UPuO2Constant");

    heatCapacity_ =
        heatCapacityModel::New(mesh, materialModelDict, "UPuO2Matpro");

    conductivity_ =
        conductivityModel::New(mesh, materialModelDict, "MaUPuO2Magni");

    emissivity_ =
        emissivityModel::New(mesh, materialModelDict, "UO2Relap");

    YoungModulus_ =
        YoungModulusModel::New(mesh, materialModelDict, "UPuO2Matpro");

    PoissonRatio_ =
        PoissonRatioModel::New(mesh, materialModelDict, "UPuO2Constant");

    thermalExpansion_ =
        thermalExpansionModel::New(mesh, materialModelDict, "UPuO2Matpro");

    densification_ = 
        densificationModel::New(mesh, materialModelDict, "UO2Frapcon");

    swelling_ = 
        swellingModel::New(mesh, materialModelDict, "UO2Frapcon");

    relocation_ = 
        relocationModel::New(mesh, materialModelDict, "UO2Frapcon");

    failure_ = 
        failureModel::New(mesh, materialModelDict, "none");

    poreVelocity_ =
    poreVelocityModel::New
    (
        mesh,
        materialModelDict, 
        "none"
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::UPuO2::~UPuO2()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::UPuO2::correctBehavioralModels
(
    const scalarField& T
)
{   
    // Update fields using behavioral models
    densification_->correct(T, addr_);
    swelling_->correct(T, addr_);
    relocation_->correct(T, addr_);
    poreVelocity_->correct(addr_);
}

void Foam::UPuO2::checkFailure()
{   
    // Check material failure
    failure_->checkFailure(addr_);
}


void Foam::UPuO2::additionalStrains
(
    symmTensorField& sf
)
const
{    
    // Reference to internal fields of dependencies
    const volScalarField& epsDen = densification_->epsilonDensification();
    const volSymmTensorField& epsSwl = swelling_->epsilonSwelling();
    const volScalarField& epsRel = relocation_->epsilonRelocation();

    forAll(addr_, i)
    {
        const label cellI = addr_[i];   

        sf[cellI] = (epsDen[cellI]*I + epsSwl[cellI]);
        sf[cellI].xx() += epsRel[cellI];
        sf[cellI].yy() += epsRel[cellI];
    }
}
// ************************************************************************* //
