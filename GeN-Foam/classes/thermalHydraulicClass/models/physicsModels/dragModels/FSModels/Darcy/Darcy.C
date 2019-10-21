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

#include "Darcy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(Darcy, 0);
    addToRunTimeSelectionTable(dragModel, Darcy, FSDragModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::Darcy::Darcy
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FSPair& FSPair,
    const word region
)
:
    dragModel
    (
        objReg,
        dict,
        FSPair,
        region
    ),
    FSPair_(FSPair),
    coeffs_("coeff", dimless, dict),
    exps_("exp", dimless, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volTensorField> Foam::dragModels::Darcy::Kd() const
{    
    tmp<volTensorField> tKd
    (
        new volTensorField
        (
            IOobject
            (
                "",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedTensor
            (
                "", 
                dimK, 
                tensor
                (
                    0, 0, 0,
                    0, 0, 0,
                    0, 0, 0
                )
            ),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volTensorField& Kd = tKd.ref();

    const volVectorField& lDh(FSPair_.structure().lDh());
    const volScalarField& magU(FSPair_.fluidRef().magU());
    tmp<volScalarField> lReX(mag(FSPair_.lRe().component(0)));
    lReX = max(lReX, dimensionedScalar("", dimless, 10));
    tmp<volScalarField> lReY(mag(FSPair_.lRe().component(1)));
    lReY = max(lReY, dimensionedScalar("", dimless, 10));
    tmp<volScalarField> lReZ(mag(FSPair_.lRe().component(2)));
    lReZ = max(lReZ, dimensionedScalar("", dimless, 10));

    volScalarField KdDarcyX
    (
        0.5*FSPair_.fluidRef().thermo().rho()*magU*
        (
            coeffs_[0]*pow(lReX, exps_[0])
        )
        /
        lDh.component(0)
    );
    Kd.replace(0, KdDarcyX);

    volScalarField KdDarcyY
    (
        0.5*FSPair_.fluidRef().thermo().rho()*magU*
        (
            coeffs_[1]*pow(lReY, exps_[1])
        )
        /
        lDh.component(1)
    );
    Kd.replace(4, KdDarcyY);

    volScalarField KdDarcyZ
    (
        0.5*FSPair_.fluidRef().thermo().rho()*magU*
        (
            coeffs_[2]*pow(lReZ, exps_[2])
        )
        /
        lDh.component(2)
    );
    Kd.replace(8, KdDarcyZ);
    
    Kd.correctBoundaryConditions();

    return tKd;
}


// ************************************************************************* //
