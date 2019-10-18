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

#include "ChurchillGunterShaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(ChurchillGunterShaw, 0);
    addToRunTimeSelectionTable(dragModel, ChurchillGunterShaw, FSDragModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::ChurchillGunterShaw::ChurchillGunterShaw
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
    surfaceRoughness_
    (
        "surfaceRoughness", dimless, dict
    )
{
    word principalAxisName(dict.get<word>("principalAxis"));
    if (principalAxisName == "localZ")
    {
        principalAxis_ = 2;
        transverseAxis1_ = 0;
        transverseAxis2_ = 1;
    }
    else if (principalAxisName == "localY")
    {
        principalAxis_ = 1;
        transverseAxis1_ = 0;
        transverseAxis2_ = 2;
    }
    else if (principalAxisName == "localX")
    {
        principalAxis_ = 0;
        transverseAxis1_ = 1;
        transverseAxis2_ = 2;
    }
    else
    {
        FatalErrorInFunction 
            << "Valid keywords for principalAxis are: localZ, localY, localX"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volTensorField> Foam::dragModels::ChurchillGunterShaw::Kd() const
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
    tmp<volScalarField> pRe(mag(FSPair_.lRe().component(principalAxis_)));
    pRe = max(pRe, dimensionedScalar("", dimless, 10));
    tmp<volScalarField> tRe1(mag(FSPair_.lRe().component(transverseAxis1_)));
    tRe1 = max(tRe1, dimensionedScalar("", dimless, 10));
    tmp<volScalarField> tRe2(mag(FSPair_.lRe().component(transverseAxis2_)));
    tRe2 = max(tRe2, dimensionedScalar("", dimless, 10));

    volScalarField KdChurchill
    (
        FSPair_.fluidRef().thermo().rho()*
        magU*
        4*
        pow
        (
            pow(8/pRe(), 12) +
            scalar(1) / 
            (
                pow
                (
                    (
                        pow
                        (
                            -2.457*
                            Foam::log
                            (
                                pow(7/pRe(), 0.9)+
                                0.27*surfaceRoughness_
                            ), 
                            16
                        ) +
                        pow
                        (
                            37530/pRe(),
                            16
                        )
                    ), 
                    1.5
                )
            ),
            scalar(1)/12
        )/
        lDh.component(principalAxis_)
    );
    Kd.replace(principalAxis_*4, KdChurchill);

    volScalarField KdGunterShaw1
    (
        FSPair_.fluidRef().thermo().rho()*
        magU*
        0.96*pow(tRe1(),-0.145)
        /lDh.component(transverseAxis1_)
    );
    Kd.replace(transverseAxis1_*4, KdGunterShaw1);

    volScalarField KdGunterShaw2
    (
        FSPair_.fluidRef().thermo().rho()*
        magU*
        0.96*pow(tRe2(),-0.145)
        /lDh.component(transverseAxis2_)
    );
    Kd.replace(transverseAxis2_*4, KdGunterShaw2);
    Kd.correctBoundaryConditions();

    return tKd;
}


// ************************************************************************* //
