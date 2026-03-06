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

#include "Steel1515Ti.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Steel1515Ti, 0);
    addToRunTimeSelectionTable(materialModel, Steel1515Ti, dictionary);

    const char* Steel1515Ti::group_ = ("materials::Steel1515Ti");

    defineParameter(Steel1515Ti, F_rho_, "F_rho", (dimless), 1.0);
    defineParameter(Steel1515Ti, F_Cp_, "F_Cp", (dimless), 1.0);
    defineParameter(Steel1515Ti, F_k_, "F_k", (dimless), 1.0);
    defineParameter(Steel1515Ti, F_emissivity_, "F_emissivity", (dimless), 1.0);
    defineParameter(Steel1515Ti, F_E_, "F_E", (dimless), 1.0);
    defineParameter(Steel1515Ti, F_nu_, "F_nu", (dimless), 1.0);
    defineParameter(Steel1515Ti, F_alphaT_, "F_alphaT", (dimless), 1.0);

    defineParameter(Steel1515Ti, delta_rho_, "delta_rho", (dimMass/dimVolume), 0.0);
    defineParameter(Steel1515Ti, delta_Cp_, "delta_Cp", (0, 2, -2, -1, 0), 0.0);
    defineParameter(Steel1515Ti, delta_k_, "delta_k", (1, 1, -3, -1, 0), 0.0);
    defineParameter(Steel1515Ti, delta_emissivity_, "delta_emissivity", (dimless), 0.0);
    defineParameter(Steel1515Ti, delta_E_, "delta_E", (dimPressure), 0.0);
    defineParameter(Steel1515Ti, delta_nu_, "delta_nu", (dimless), 0.0);
    defineParameter(Steel1515Ti, delta_alphaT_, "delta_alphaT", (dimless), 0.0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Steel1515Ti::Steel1515Ti
(
    const fvMesh& mesh, 
    const dictionary& materialModelDict,
    const labelList& addr
)
:
    materialModel(mesh, materialModelDict, addr),
    swelling_(swellingModel::New(mesh, materialModelDict, "Steel1515TiGeneralized"))
{
    density_ =
        densityModel::New(mesh, materialModelDict, "Steel1515TiSchumann");
    
    heatCapacity_ =
        heatCapacityModel::New(mesh, materialModelDict, "Steel1515TiBanerjee");
    
    conductivity_ =
        conductivityModel::New(mesh, materialModelDict, "Steel1515TiTobbe");
    
    emissivity_ =
        emissivityModel::New(mesh, materialModelDict, "ZircaloyConstant");
    
    YoungModulus_ =
        YoungModulusModel::New(mesh, materialModelDict, "Steel1515TiTobbe");
    
    PoissonRatio_ =
        PoissonRatioModel::New(mesh, materialModelDict, "Steel1515TiTobbe");
    
    thermalExpansion_ =
        thermalExpansionModel::New(mesh, materialModelDict, "Steel1515TiGehr");
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Steel1515Ti::~Steel1515Ti()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::Steel1515Ti::correctBehavioralModels
(
    const scalarField& T
)
{   
    // //- Update fields using behavioral models
    swelling_->correct(T, addr_);
}


void Foam::Steel1515Ti::additionalStrains
(
    symmTensorField& sf
)
const
{    
    //- Reference to internal fields of dependencies
    const symmTensorField& epsSwl = 
    swelling_->epsilonSwelling();

    forAll(addr_, i)
    {
        const label cellI = addr_[i];

        sf[cellI] = (epsSwl[cellI]);
    }
}
