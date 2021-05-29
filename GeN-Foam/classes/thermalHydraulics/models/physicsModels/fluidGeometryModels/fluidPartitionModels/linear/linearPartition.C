/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "linearPartition.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidPartitionModels
{
    defineTypeNameAndDebug(linearPartition, 0);
    addToRunTimeSelectionTable
    (
        fluidPartitionModel, 
        linearPartition, 
        fluidPartitionModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidPartitionModels::linearPartition::linearPartition
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const fluid& dispersed,
    const fluid& continuous,
    const structureModel& structure
)
:
    fluidPartitionModel
    (
        objReg,
        dict,
        dispersed,
        continuous,
        structure
    ),
    fracD_
    (
        IOobject
        (
            IOobject::groupName("frac."+dispersed.name(), typeName),
            mesh_.time().timeName(),
            objReg
        ),
        mesh_,
        dimensionedScalar(0),
        zeroGradientFvPatchScalarField::typeName
    ),
    fracC_
    (
        IOobject
        (
            IOobject::groupName("frac."+continuous.name(), typeName),
            mesh_.time().timeName(),
            objReg
        ),
        mesh_,
        dimensionedScalar(0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    /*
    if (dict.found(dispersed.name()))
    {
        dimensionedScalar valueD(dispersed.name(), dimless, dict);
        forAll(mesh_.cells(), celli)
        {
            fracD_[celli] = valueD.value();
        }
        if (dict.found(continuous.name()))
        {
            dimensionedScalar valueC(continuous.name(), dimless, dict);
            forAll(mesh_.cells(), celli)
            {
                fracC_[celli] = valueC.value();
            }
        }
        else
        {
            fracC_ = 1.0 - fracD_;
        }
    }
    else if (dict.found(continuous.name()))
    {
        dimensionedScalar valueC(continuous.name(), dimless, dict);
        forAll(mesh_.cells(), celli)
        {
            fracC_[celli] = valueC.value();
        }
        fracD_ = 1.0 - fracC_;
    }
    else
    {
        FatalErrorInFunction
            << "No values provided for interfacial area fractions!"
            << exit(FatalError);
    }
    
    fracD_.correctBoundaryConditions();
    fracC_.correctBoundaryConditions();
    volScalarField fracSum(fracD_+fracC_);
    
    if 
    ( 
        (max(fracSum).value() != 1) 
        or 
        (min(fracSum).value() != 1)
        or 
        (max(fracD_).value() > 1)
        or 
        (min(fracD_).value() < 0)
        or 
        (max(fracC_).value() > 1)
        or 
        (min(fracC_).value() < 0)
    )
    {
        FatalErrorInFunction
            << "Provided interfacial area fractions "
            << "do not sum up to 1 or are outside bounds [0,1]!"
            << exit(FatalError);
    }
    */
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidPartitionModels::linearPartition::~linearPartition()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::fluidPartitionModels::linearPartition::fracD() const
{
    return dispersed_.normalized();
}


Foam::tmp<Foam::volScalarField> 
Foam::fluidPartitionModels::linearPartition::fracC() const
{
    return continuous_.normalized();
}


// ************************************************************************* //
