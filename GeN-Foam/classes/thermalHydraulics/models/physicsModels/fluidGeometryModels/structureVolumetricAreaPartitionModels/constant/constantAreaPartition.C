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

#include "constantAreaPartition.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace structureVolumetricAreaPartitionModels
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable
    (
        structureVolumetricAreaPartitionModel, 
        constant, 
        structureVolumetricAreaPartitionModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::structureVolumetricAreaPartitionModels::constant::constant
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const fluid& fluid1,
    const fluid& fluid2,
    const structureModel& structure
)
:
    structureVolumetricAreaPartitionModel
    (
        objReg,
        dict,
        fluid1,
        fluid2,
        structure
    ),
    frac1_
    (
        IOobject
        (
            IOobject::groupName("frac."+fluid1.name(), typeName),
            mesh_.time().timeName(),
            objReg
        ),
        mesh_,
        dimensionedScalar(0),
        zeroGradientFvPatchScalarField::typeName
    ),
    frac2_
    (
        IOobject
        (
            IOobject::groupName("frac."+fluid2.name(), typeName),
            mesh_.time().timeName(),
            objReg
        ),
        mesh_,
        dimensionedScalar(0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    if (dict.found("value1"))
    {
        dimensionedScalar value1("value1", dimless, dict);
        forAll(mesh_.cells(), celli)
        {
            frac1_[celli] = value1.value();
        }
        if (dict.found("value2"))
        {
            dimensionedScalar value2("value2", dimless, dict);
            forAll(mesh_.cells(), celli)
            {
                frac2_[celli] = value2.value();
            }
        }
        else
        {
            frac2_ = 1.0 - frac1_;
        }
    }
    else if (dict.found("value2"))
    {
        dimensionedScalar value2("value2", dimless, dict);
        forAll(mesh_.cells(), celli)
        {
            frac2_[celli] = value2.value();
        }
        frac1_ = 1.0 - frac2_;
    }
    else
    {
        FatalErrorInFunction
            << "No values provided for interfacial area fractions!"
            << exit(FatalError);
    }
    
    frac1_.correctBoundaryConditions();
    frac2_.correctBoundaryConditions();
    volScalarField fracSum(frac1_+frac2_);
    
    if 
    ( 
        (max(fracSum).value() != 1) 
        or 
        (min(fracSum).value() != 1)
        or 
        (max(frac1_).value() > 1)
        or 
        (min(frac1_).value() < 0)
        or 
        (max(frac2_).value() > 1)
        or 
        (min(frac2_).value() < 0)
    )
    {
        FatalErrorInFunction
            << "Provided interfacial area fractions "
            << "do not sum up to 1 or are outside bounds [0,1]!"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::structureVolumetricAreaPartitionModels::constant::~constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::structureVolumetricAreaPartitionModels::constant::frac1() const
{
    tmp<volScalarField> tfrac1
    (
        new volScalarField
        (
            frac1_
        )
    );
    
    return tfrac1;
}


Foam::tmp<Foam::volScalarField> 
Foam::structureVolumetricAreaPartitionModels::constant::frac2() const
{
    tmp<volScalarField> tfrac2
    (
        new volScalarField
        (
            frac2_
        )
    );
    
    return tfrac2;
}


// ************************************************************************* //
