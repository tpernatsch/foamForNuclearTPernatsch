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

#include "NusseltHeatTransfer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(Nusselt, 0);
    addToRunTimeSelectionTable
    (
        heatTransferModel, 
        Nusselt, 
        FFHeatTransferModels
    );
    addToRunTimeSelectionTable
    (
        heatTransferModel, 
        Nusselt, 
        FSHeatTransferModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::Nusselt::Nusselt
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FFPair& FFPair,
    const word nameBulk,
    const word nameInterface
)
:
    heatTransferModel
    (
        objReg,
        dict,
        FFPair,
        nameBulk,
        nameInterface
    ),
    fluid_
    (
        (
            FFPair.fluid1().name() == nameBulk
        ) ?
        FFPair.fluid1() :
        FFPair.fluid2()
    ),
    Re_
    (
        FFPair.Re()
    ),
    Dh_
    (
        FFPair.DhDispersed()
    ),
    A_
    (
        dict.get<scalar>("const")
    ),
    B_
    (
        dict.get<scalar>("coeff")
    ),
    C_
    (
        dict.get<scalar>("expRe")
    ),
    D_
    (
        dict.get<scalar>("expPr")
    ),
    //- I need to init a ref to something, so I init it to
    //  the dispersion field which is otherwise not needed
    //  for an FF heatTransferModel
    cellField_(fluid_.dispersion())
{
    if (FFPair.fluid2().name() != nameBulk)
    {
        FatalErrorInFunction
            << "Phase " << nameBulk << " not found!"
            << exit(FatalError);
    }
}


Foam::heatTransferModels::Nusselt::Nusselt
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FSPair& FSPair,
    const word region
)
:
    heatTransferModel
    (
        objReg,
        dict,
        FSPair,
        region
    ),
    fluid_
    (
        FSPair.fluidRef()
    ),
    Re_
    (
        FSPair.Re()
    ),
    Dh_
    (
        FSPair.fluidRef().Dh()
    ),
    A_
    (
        dict.get<scalar>("const")
    ),
    B_
    (
        dict.get<scalar>("coeff")
    ),
    C_
    (
        dict.get<scalar>("expRe")
    ),
    D_
    (
        dict.get<scalar>("expPr")
    ),
    cellField_
    (
        FSPair.structure().cellFields()[region]
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::heatTransferModels::Nusselt::htc() const
{
    tmp<volScalarField> thtc
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
            dimensionedScalar("", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& htc = thtc.ref();

    //- Prandtl of bulk fluid computed on spot
    htc = 
            A_ 
        +   B_
        *   pow(Re_, C_)
        *   pow
            (
                    fluid_.thermo().Cp()
                *   fluid_.thermo().mu()
                /   fluid_.thermo().kappa(),
                D_
            );

    htc *= fluid_.thermo().kappa()/Dh_;

    if (withStructure_)
    {
        htc.primitiveFieldRef() *= cellField_;
        htc.correctBoundaryConditions();
    }

    return thtc;
}


// ************************************************************************* //
