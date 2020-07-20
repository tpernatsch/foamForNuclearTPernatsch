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

#include "AzadShirani.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace structureVolumetricAreaPartitionModels
{
    defineTypeNameAndDebug(AzadShirani, 0);
    addToRunTimeSelectionTable
    (
        structureVolumetricAreaPartitionModel, 
        AzadShirani, 
        structureVolumetricAreaPartitionModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::structureVolumetricAreaPartitionModels::AzadShirani::AzadShirani
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
    alpha_
    (
        (fluid1.isGas()) ? 
        fluid1.normalized() : fluid2.normalized()
    ),
    x_
    (
        (fluid1.isGas()) ? 
        fluid1.flowQuality() : fluid2.flowQuality()
    ),
    alphaCrit_
    (
        dict.lookupOrDefault<scalar>("alphaCrit", 0.957)
    ),
    xCrit_
    (
        dict.lookupOrDefault<scalar>("xCrit", 0.3)
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
    if 
    (
        !(fluid1.isLiquid() and fluid2.isGas()) and
        !(fluid2.isLiquid() and fluid1.isGas())
    )
    {
        FatalErrorInFunction
            << "The AzadShirani model only works for liquid-gas systems. Set "
            << "the stateOfMatter entry in "
            << "phaseProperties." << fluid1.name() << "Properties and/or "
            << "phaseProperties." << fluid2.name() << "Properties) to "
            << "distinguish between gas and liquid"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::structureVolumetricAreaPartitionModels::AzadShirani::~AzadShirani()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::structureVolumetricAreaPartitionModels::AzadShirani::frac1() const
{
    tmp<volScalarField> tfrac1
    (
        new volScalarField
        (
            IOobject
            (
                "",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("", dimless, 0.0)
        )
    );
    volScalarField& frac1 = tfrac1.ref();
    frac1 = 
        max
        (
            pow((1.0-alpha_)/(1.0-alphaCrit_),0.5), 1.0
        )*
        max
        (
            pow((1-x_)/(1-xCrit_), 1.5), 1.0
        );
    
    return tfrac1;
}


Foam::tmp<Foam::volScalarField> 
Foam::structureVolumetricAreaPartitionModels::AzadShirani::frac2() const
{    
    return 1.0-frac1();
}


// ************************************************************************* //
