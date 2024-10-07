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

#include "materialModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(materialModel, 0);
    defineRunTimeSelectionTable(materialModel, dictionary);  

    const char* materialModel::group_ = ("materials");

    defineParameter(materialModel, F_rho_, "F_rho", (dimless), 1.0);
    defineParameter(materialModel, F_Cp_, "F_Cp", (dimless), 1.0);
    defineParameter(materialModel, F_k_, "F_k", (dimless), 1.0);
    defineParameter(materialModel, F_emissivity_, "F_emissivity", (dimless), 1.0);
    defineParameter(materialModel, F_E_, "F_E", (dimless), 1.0);
    defineParameter(materialModel, F_nu_, "F_nu", (dimless), 1.0);
    defineParameter(materialModel, F_alphaT_, "F_alphaT", (dimless), 1.0);

    defineParameter(materialModel, delta_rho_, "delta_rho", (dimMass/dimVolume), 0.0);
    defineParameter(materialModel, delta_Cp_, "delta_Cp", (0, 2, -2, -1, 0), 0.0);
    defineParameter(materialModel, delta_k_, "delta_k", (1, 1, -3, -1, 0), 0.0);
    defineParameter(materialModel, delta_emissivity_, "delta_emissivity", (dimless), 0.0);
    defineParameter(materialModel, delta_E_, "delta_E", (dimPressure), 0.0);
    defineParameter(materialModel, delta_nu_, "delta_nu", (dimless), 0.0);  
    defineParameter(materialModel, delta_alphaT_, "delta_alphaT", (dimless), 0.0);  
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::materialModel::materialModel
(
    const fvMesh& mesh, 
    const dictionary& materialModelDict,
    const labelList& addr
)
:
    mesh_(mesh),
    materialModelDict_(materialModelDict), 
    addr_(addr), 
    density_(),
    heatCapacity_(),
    conductivity_(),
    emissivity_(),
    YoungModulus_(),
    PoissonRatio_(),
    thermalExpansion_()
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::materialModel>
Foam::materialModel::New
(
    const fvMesh& mesh, 
    const dictionary& materialModelDict,
    const labelList& addr
)
{
    word materialName;
    materialModelDict.lookup("material") >> materialName;

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(materialName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("materialModel::New(const fvMesh&, const dictionary&)")
            << "Unknown materialModel type "
            << materialName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }
    
    Info<< "     Selecting material model --> "
    << materialName << endl;

    if (debug)
    {
        Info<< "Selecting materialModel type "
            << materialName << endl;
    }

    autoPtr<Foam::materialModel> 
    materialModelPtr(cstrIter()(mesh, materialModelDict, addr));

    return materialModelPtr;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::materialModel::~materialModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::materialModel::rho
(
    scalarField& sf, 
    const scalarField& T
) 
{
    density_->correct(sf, T, addr_);

    const dimensionedScalar& F = 
    isInUserParameters(this->group(), "F_rho") ?
    this->F_rho() : materialModel::F_rho();

    const dimensionedScalar& delta = 
    isInUserParameters(this->group(), "delta_rho") ?
    this->delta_rho() : materialModel::delta_rho();

    applyParameters(sf, addr_, F, delta);
}

void Foam::materialModel::Cp
(
    scalarField& sf, 
    const scalarField& T
) 
{
    heatCapacity_->correct(sf, T, addr_);

    const dimensionedScalar& F = 
    isInUserParameters(this->group(), "F_Cp") ?
    this->F_Cp() : materialModel::F_Cp();

    const dimensionedScalar& delta = 
    isInUserParameters(this->group(), "delta_Cp") ?
    this->delta_Cp() : materialModel::delta_Cp();

    applyParameters(sf, addr_, F, delta);
}

void Foam::materialModel::k
(
    scalarField& sf, 
    const scalarField& T
) 
{
    conductivity_->correct(sf, T, addr_);

    const dimensionedScalar& F = 
    isInUserParameters(this->group(), "F_k") ?
    this->F_k() : materialModel::F_k();

    const dimensionedScalar& delta = 
    isInUserParameters(this->group(), "delta_k") ?
    this->delta_k() : materialModel::delta_k();

    applyParameters(sf, addr_, F, delta);
}

void Foam::materialModel::emissivity
(
    scalarField& sf, 
    const scalarField& T
) 
{
    emissivity_->correct(sf, T, addr_);

    const dimensionedScalar& F = 
    isInUserParameters(this->group(), "F_emissivity") ?
    this->F_emissivity() : materialModel::F_emissivity();

    const dimensionedScalar& delta = 
    isInUserParameters(this->group(), "delta_emissivity") ?
    this->delta_emissivity() : materialModel::delta_emissivity();

    applyParameters(sf, addr_, F, delta);
}


void Foam::materialModel::E
(
    scalarField& sf, 
    const scalarField& T
) 
{
    YoungModulus_->correct(sf, T, addr_);

    const dimensionedScalar& F = 
    isInUserParameters(this->group(), "F_E") ?
    this->F_E() : materialModel::F_E();

    const dimensionedScalar& delta = 
    isInUserParameters(this->group(), "delta_E") ?
    this->delta_E() : materialModel::delta_E();

    applyParameters(sf, addr_, F, delta);
}

void Foam::materialModel::nu
(
    scalarField& sf, 
    const scalarField& T
) 
{
    PoissonRatio_->correct(sf, T, addr_);

    const dimensionedScalar& F = 
    isInUserParameters(this->group(), "F_nu") ?
    this->F_nu() : materialModel::F_nu();

    const dimensionedScalar& delta = 
    isInUserParameters(this->group(), "delta_nu") ?
    this->delta_nu() : materialModel::delta_nu();

    applyParameters(sf, addr_, F, delta);
}


void Foam::materialModel::alphaT
(
    symmTensorField& sf, 
    const scalarField& T,
    const labelList& addr
) const
{    
    thermalExpansion_->correct(sf, T, addr);

    const dimensionedScalar& F = 
    isInUserParameters(this->group(), "F_alphaT") ?
    this->F_alphaT() : materialModel::F_alphaT();

    const dimensionedScalar& delta = 
    isInUserParameters(this->group(), "delta_alphaT") ?
    this->delta_alphaT() : materialModel::delta_alphaT();

    applyParameters(sf, addr, F, delta);
}


// ************************************************************************* //
