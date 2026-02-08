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

#include "constantMaterial.H"
#include "addToRunTimeSelectionTable.H"
#include "userParameters.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantMaterial, 0);
    addToRunTimeSelectionTable(materialModel, constantMaterial, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantMaterial::constantMaterial
(
    const fvMesh& mesh, 
    const dictionary& materialModelDict,
    const labelList& addr
)
:
    materialModel(mesh, materialModelDict, addr)
{   
    density_ =
        densityModel::New(mesh, materialModelDict, "constant");

    heatCapacity_ =
        heatCapacityModel::New(mesh, materialModelDict, "constant");

    conductivity_ =
        conductivityModel::New(mesh, materialModelDict, "constant");

    emissivity_ =
        emissivityModel::New(mesh, materialModelDict, "constant");

    YoungModulus_ =
        YoungModulusModel::New(mesh, materialModelDict, "constant");

    PoissonRatio_ =
        PoissonRatioModel::New(mesh, materialModelDict, "constant");

    thermalExpansion_ =
        thermalExpansionModel::New(mesh, materialModelDict, "constant");
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantMaterial::~constantMaterial()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
