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

#include "constantPower.H"
#include "addToRunTimeSelectionTable.H"
#include "structureModel.H"
#include "fvMatrix.H"
#include "fvmDdt.H"
#include "fvmSup.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerModels
{
    defineTypeNameAndDebug(constantPower, 0);
    addToRunTimeSelectionTable
    (
        powerModel, 
        constantPower, 
        powerModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModels::constantPower::constantPower
(
    const structureModel& structure,
    const dictionary& dict,
    const word name
)
:
    powerModel
    (
        structure,
        dict,
        name
    ),
    powerDensity_
    (
        IOobject
        (
            "powerDensity."+typeName+"."+name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("powerDensity", dimPower/dimVol, dict),
        zeroGradientFvPatchScalarField::typeName

    ),
    T_
    (
        IOobject
        (
            "T."+typeName+"."+name,
            structure.mesh().time().timeName(),
            structure.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("T", dimTemperature, dict),
        zeroGradientFvPatchScalarField::typeName
    ),
    Cp_("Cp", dimEnergy/dimMass/dimTemperature, dict),
    rho_("rho", dimMass/dimVol, dict)
{
    powerDensity_.primitiveFieldRef() *= cellField_;
    powerDensity_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerModels::constantPower::~constantPower()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::powerModels::constantPower::T() const
{
    tmp<volScalarField> tT
    (
        new volScalarField(T_)
    );
    volScalarField& T = tT.ref();
    T.primitiveFieldRef() *= cellField_;
    T.correctBoundaryConditions();
    return tT;
}


void Foam::powerModels::constantPower::correct
(
    const volScalarField& HTSum,  // == SUM_j [htc_j*T_j*frac_j]
    const volScalarField& HSum    // == SUM_j [htc_j*frac_j]
)
{
    volScalarField alpha
    (
        Foam::max
        (
            structure_.alphaFields()[name_],
            dimensionedScalar("", dimless, 1e-69)
        )
    );

    fvScalarMatrix EEqn
    (
            fvm::ddt(alpha*rho_*Cp_, T_)
        ==
            iA_*HTSum
        -   fvm::Sp(iA_*HSum, T_)
        +   alpha*powerDensity_
    );
    EEqn.solve();
    T_.correctBoundaryConditions();
}


// ************************************************************************* //
