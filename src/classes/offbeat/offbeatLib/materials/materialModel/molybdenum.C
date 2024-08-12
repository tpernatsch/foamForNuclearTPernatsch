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

#include "molybdenum.H"
#include "addToRunTimeSelectionTable.H"
#include "userParameters.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(molybdenum, 0);
    addToRunTimeSelectionTable(materialModel, molybdenum, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::molybdenum::molybdenum
(
    const fvMesh& mesh, 
    const dictionary& materialModelDict,
    const labelList& addr
)
:
    materialModel(mesh, materialModelDict, addr),
    swelling_(swellingModel::New(mesh, materialModelDict, "FeCrAl"))
{  
    density_ =
        densityModel::New(mesh, materialModelDict, "MoConstant");
    
    heatCapacity_ =
        heatCapacityModel::New(mesh, materialModelDict, "Mo");
    
    conductivity_ =
        conductivityModel::New(mesh, materialModelDict, "Mo");
    
    emissivity_ =
        emissivityModel::New(mesh, materialModelDict, "MoConstant");
    
    YoungModulus_ =
        YoungModulusModel::New(mesh, materialModelDict, "Mo");
    
    PoissonRatio_ =
        PoissonRatioModel::New(mesh, materialModelDict, "MoConstant");
    
    thermalExpansion_ =
        thermalExpansionModel::New(mesh, materialModelDict, "Mo");  
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::molybdenum::~molybdenum()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::molybdenum::correctBehavioralModels
(
    const scalarField& T
)
{   
    //- Update fields using behavioral models
    swelling_->correct(T, addr_);
}


void Foam::molybdenum::additionalStrains
(
    symmTensorField& sf
)
const
{    
    //- Reference to internal fields of dependencies
    const symmTensorField& epsSwl = swelling_->epsilonSwelling();

    forAll(addr_, i)
    {
        const label cellI = addr_[i];

        sf[cellI] = (epsSwl[cellI]);
    }
}
// ************************************************************************* //
