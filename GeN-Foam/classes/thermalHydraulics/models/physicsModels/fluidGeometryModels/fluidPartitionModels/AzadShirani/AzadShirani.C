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
namespace fluidPartitionModels
{
    defineTypeNameAndDebug(AzadShirani, 0);
    addToRunTimeSelectionTable
    (
        fluidPartitionModel, 
        AzadShirani, 
        fluidPartitionModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidPartitionModels::AzadShirani::AzadShirani
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
    alpha_
    (
        (dispersed.isGas()) ? 
        dispersed.normalized() : continuous.normalized()
    ),
    x_
    (
        (dispersed.isGas()) ? 
        dispersed.flowQuality() : continuous.flowQuality()
    ),
    alphaCrit_
    (
        dict.lookupOrDefault<scalar>("alphaCrit", 0.957)
    ),
    xCrit_
    (
        dict.lookupOrDefault<scalar>("xCrit", 0.3)
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
    if 
    (
        !(dispersed.isLiquid() and continuous.isGas()) and
        !(continuous.isLiquid() and dispersed.isGas())
    )
    {
        FatalErrorInFunction
            << "The AzadShirani model only works for liquid-gas systems. Set "
            << "the stateOfMatter entry in "
            << "phaseProperties." << dispersed.name() << "Properties and/or "
            << "phaseProperties." << continuous.name() << "Properties) to "
            << "distinguish between gas and liquid"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidPartitionModels::AzadShirani::~AzadShirani()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::fluidPartitionModels::AzadShirani::fracC() const
{
    tmp<volScalarField> tfracC
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
    volScalarField& fracC = tfracC.ref();
    fracC = 
        min
        (
            pow((1.0-alpha_)/(1.0-alphaCrit_),0.5), 1.0
        )*
        min
        (
            pow((1-x_)/(1-xCrit_), 1.5), 1.0
        );
    
    return tfracC;
}


Foam::tmp<Foam::volScalarField> 
Foam::fluidPartitionModels::AzadShirani::fracD() const
{    
    return 1.0-fracC();
}


// ************************************************************************* //
