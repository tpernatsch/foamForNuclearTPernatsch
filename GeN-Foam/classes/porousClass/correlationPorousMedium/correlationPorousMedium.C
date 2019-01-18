/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "correlationPorousMedium.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(correlationPorousMedium, 0);

    addToRunTimeSelectionTable
    (
        porousMedium,
        correlationPorousMedium,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::correlationPorousMedium::correlationPorousMedium
(
    const volScalarField& rho,
    const volVectorField& U,
    const rhoThermo& thermo
)
:
    porousMedium(rho, U, thermo),//correlationPorousMedium is derived from porousMedium -> use a porousMedium constructor to initialize that part
    reynoldsTurb_(
        IOobject(
            "porousMedium::reynoldsTurb",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("reynoldsTurb")),
        zeroGradientFvPatchScalarField::typeName
    ),
    reynoldsLam_(
        IOobject(
            "porousMedium::reynoldsLam",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("reynoldsLam")),
        zeroGradientFvPatchScalarField::typeName
    ),
    darcyConstTurb_(
        IOobject(
            "porousMedium::darcyConstTurb",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("darcyConstTurb")),
        zeroGradientFvPatchScalarField::typeName
    ),
    darcyConstLam_(
        IOobject(
            "porousMedium::darcyConstLam",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("darcyConstLam")),
        zeroGradientFvPatchScalarField::typeName
    ),
    darcyExpTurb_(
        IOobject(
            "porousMedium::darcyExpTurb",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("darcyExpTurb")),
        zeroGradientFvPatchScalarField::typeName
    ),
    darcyExpLam_(
        IOobject(
            "porousMedium::darcyExpLam",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("darcyExpLam")),
        zeroGradientFvPatchScalarField::typeName
    ),
    nusseltConstTurb1_(
        IOobject(
            "porousMedium::nusseltConstTurb1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("nusseltConstTurb1")),
        zeroGradientFvPatchScalarField::typeName
    ),
    nusseltConstLam1_(
        IOobject(
            "porousMedium::nusseltConstLam1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("nusseltConstLam1")),
        zeroGradientFvPatchScalarField::typeName
    ),
    nusseltConstTurb2_(
        IOobject(
            "porousMedium::nusseltConstTurb2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("nusseltConstTurb2")),
        zeroGradientFvPatchScalarField::typeName
    ),
    nusseltConstLam2_(
        IOobject(
            "porousMedium::nusseltConstLam2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("nusseltConstLam2")),
        zeroGradientFvPatchScalarField::typeName
    ),
    nusseltReExpTurb_(
        IOobject(
            "porousMedium::nusseltReExpTurb",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("nusseltReExpTurb")),
        zeroGradientFvPatchScalarField::typeName
    ),
    nusseltReExpLam_(
        IOobject(
            "porousMedium::nusseltReExpLam",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("nusseltReExpLam")),
        zeroGradientFvPatchScalarField::typeName
    ),
    nusseltPrExpTurb_(
        IOobject(
            "porousMedium::nusseltPrExpTurb",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("nusseltPrExpTurb")),
        zeroGradientFvPatchScalarField::typeName
    ),
    nusseltPrExpLam_(
        IOobject(
            "porousMedium::nusseltPrExpLam",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("nusseltPrExpLam")),
        zeroGradientFvPatchScalarField::typeName
    ),
    hydraulicDiameter_(
        IOobject(
            "porousMedium::hydraulicDiameter",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(IOdictionary::lookup("hydraulicDiameter")),
        zeroGradientFvPatchScalarField::typeName
    ),
    hydraulicDiameterStructure_(
        IOobject(
            "porousMedium::hydraulicDiameterStructure",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(IOdictionary::lookup("hydraulicDiameterStructure")),
        zeroGradientFvPatchScalarField::typeName
    ),
    pumpMomentumSource_(
        IOobject(
            "porousMedium::pumpMomentumSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("pumpMomentumSource")),
        zeroGradientFvPatchScalarField::typeName
    ),
    turbulenceIntensityConst_(
        IOobject(
            "porousMedium::turbulenceIntensityConst",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("turbulenceIntensityConst")),
        zeroGradientFvPatchScalarField::typeName
    ),
    turbulenceIntensityExp_(
        IOobject(
            "porousMedium::turbulenceIntensityExp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("turbulenceIntensityExp")),
        zeroGradientFvPatchScalarField::typeName
    ),
    turbulenceLengthScaleConst_(
        IOobject(
            "porousMedium::turbulenceLengthScaleConst",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("turbulenceLengthScaleConst")),
        zeroGradientFvPatchScalarField::typeName
    ),
    kepsilonConvergenceRate_(
        IOobject(
            "porousMedium::kepsilonConvergenceRate",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(IOdictionary::lookup("kepsilonConvergenceRate")),
        zeroGradientFvPatchScalarField::typeName
    ),
    turbulenceIntensity_(
        IOobject(
            "porousMedium::turbulenceIntensity",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(word(), dimless, vector::zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    fRe_(
        IOobject(
            "porousMedium::fRe",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(word(), dimless, scalar(0)),
        zeroGradientFvPatchScalarField::typeName
    )
{
    porousMedium::gamma_ = dimensionedScalar(IOdictionary::lookup("voidFraction"));
    porousMedium::volumetricAreaSS_ = dimensionedScalar(IOdictionary::lookup("volumetricAreaSS"));
    porousMedium::volumetricAreaFuel_ = dimensionedScalar(IOdictionary::lookup("volumetricAreaFuel"));
    porousMedium::volumetricAreaHX_ = dimensionedScalar(IOdictionary::lookup("volumetricAreaHX"));
    porousMedium::rotate_ = dimensionedTensor(coordinateSystem("",vector::zero,IOdictionary::lookup("localZaxis"),IOdictionary::lookup("localXaxis")).R().R().T());
    porousMedium::externalHeatTransferCoefficient_ = dimensionedScalar(IOdictionary::lookup("externalHeatTransferCoefficient"));
    porousMedium::externalT_ = dimensionedScalar(IOdictionary::lookup("externalT"));
    porousMedium::externalVolHeatSource_ = dimensionedScalar(IOdictionary::lookup("externalVolHeatSource"));
    porousMedium::externalRhoCp_ = dimensionedScalar(IOdictionary::lookup("externalRhoCp"));
}



Foam::correlationPorousMedium::correlationPorousMedium
(
    const volScalarField& rho,
    const volVectorField& U,
    const rhoThermo& thermo,
    const dimensionedScalar& gamma,
    const dimensionedScalar& volumetricAreaSS,
    const dimensionedScalar& volumetricAreaFuel,
    const dimensionedScalar& volumetricAreaHX,
    const dimensionedTensor& rotate,
    const dimensionedScalar& externalHeatTransferCoefficient,
    const dimensionedScalar& externalT,
    const dimensionedScalar& externalVolHeatSource,
    const dimensionedScalar& externalRhoCp,
    const dimensionedVector& reynoldsTurb,
    const dimensionedVector& reynoldsLam,
    const dimensionedVector& darcyConstTurb,
    const dimensionedVector& darcyConstLam,
    const dimensionedVector& darcyExpTurb,
    const dimensionedVector& darcyExpLam,
    const dimensionedVector& nusseltConstTurb1,
    const dimensionedVector& nusseltConstLam1,
    const dimensionedVector& nusseltConstTurb2,
    const dimensionedVector& nusseltConstLam2,
    const dimensionedVector& nusseltReExpTurb,
    const dimensionedVector& nusseltReExpLam,
    const dimensionedVector& nusseltPrExpTurb,
    const dimensionedVector& nusseltPrExpLam,
    const dimensionedScalar& hydraulicDiameter,
    const dimensionedScalar& hydraulicDiameterStructure,
    const dimensionedVector& pumpMomentumSource,
    const dimensionedVector& turbulenceIntensityConst,
    const dimensionedVector& turbulenceIntensityExp,
    const dimensionedVector& turbulenceLengthScaleConst,
    const dimensionedVector& kepsilonConvergenceRate
)
:
    porousMedium(rho, U, thermo),//correlationPorousMedium is derived from porousMedium -> use a porousMedium constructor to initialize that part
    reynoldsTurb_(
        IOobject(
            "porousMedium::reynoldsTurb",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        reynoldsTurb,
        zeroGradientFvPatchScalarField::typeName
    ),
    reynoldsLam_(
        IOobject(
            "porousMedium::reynoldsLam",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        reynoldsLam,
        zeroGradientFvPatchScalarField::typeName
    ),
    darcyConstTurb_(
        IOobject(
            "porousMedium::darcyConstTurb",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        darcyConstTurb,
        zeroGradientFvPatchScalarField::typeName
    ),
    darcyConstLam_(
        IOobject(
            "porousMedium::darcyConstLam",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        darcyConstLam,
        zeroGradientFvPatchScalarField::typeName
    ),
    darcyExpTurb_(
        IOobject(
            "porousMedium::darcyExpTurb",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        darcyExpTurb,
        zeroGradientFvPatchScalarField::typeName
    ),
    darcyExpLam_(
        IOobject(
            "porousMedium::darcyExpLam",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        darcyExpLam,
        zeroGradientFvPatchScalarField::typeName
    ),
    nusseltConstTurb1_(
        IOobject(
            "porousMedium::nusseltConstTurb1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        nusseltConstTurb1,
        zeroGradientFvPatchScalarField::typeName
    ),
    nusseltConstLam1_(
        IOobject(
            "porousMedium::nusseltConstLam1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        nusseltConstLam1,
        zeroGradientFvPatchScalarField::typeName
    ),
    nusseltConstTurb2_(
        IOobject(
            "porousMedium::nusseltConstTurb2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        nusseltConstTurb2,
        zeroGradientFvPatchScalarField::typeName
    ),
    nusseltConstLam2_(
        IOobject(
            "porousMedium::nusseltConstLam2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        nusseltConstLam2,
        zeroGradientFvPatchScalarField::typeName
    ),
    nusseltReExpTurb_(
        IOobject(
            "porousMedium::nusseltReExpTurb",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        nusseltReExpTurb,
        zeroGradientFvPatchScalarField::typeName
    ),
    nusseltReExpLam_(
        IOobject(
            "porousMedium::nusseltReExpLam",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        nusseltReExpLam,
        zeroGradientFvPatchScalarField::typeName
    ),
    nusseltPrExpTurb_(
        IOobject(
            "porousMedium::nusseltPrExpTurb",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        nusseltPrExpTurb,
        zeroGradientFvPatchScalarField::typeName
    ),
    nusseltPrExpLam_(
        IOobject(
            "porousMedium::nusseltPrExpLam",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        nusseltPrExpLam,
        zeroGradientFvPatchScalarField::typeName
    ),
    hydraulicDiameter_(
        IOobject(
            "porousMedium::hydraulicDiameter",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        hydraulicDiameter,
        zeroGradientFvPatchScalarField::typeName
    ),
    hydraulicDiameterStructure_(
        IOobject(
            "porousMedium::hydraulicDiameterStructure",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        hydraulicDiameterStructure,
        zeroGradientFvPatchScalarField::typeName
    ),
    pumpMomentumSource_(
        IOobject(
            "porousMedium::pumpMomentumSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        pumpMomentumSource,
        zeroGradientFvPatchScalarField::typeName
    ),
    turbulenceIntensityConst_(
        IOobject(
            "porousMedium::turbulenceIntensityConst",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        turbulenceIntensityConst,
        zeroGradientFvPatchScalarField::typeName
    ),
    turbulenceIntensityExp_(
        IOobject(
            "porousMedium::turbulenceIntensityExp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        turbulenceIntensityExp,
        zeroGradientFvPatchScalarField::typeName
    ),
    turbulenceLengthScaleConst_(
        IOobject(
            "porousMedium::turbulenceLengthScaleConst",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        turbulenceLengthScaleConst,
        zeroGradientFvPatchScalarField::typeName
    ),
    kepsilonConvergenceRate_(
        IOobject(
            "porousMedium::kepsilonConvergenceRate",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        kepsilonConvergenceRate,
        zeroGradientFvPatchScalarField::typeName
    ),
    turbulenceIntensity_(
        IOobject(
            "porousMedium::turbulenceIntensity",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(word(), dimless, vector::zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    fRe_(
        IOobject(
            "porousMedium::fRe",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(word(), dimless, scalar(0)),
        zeroGradientFvPatchScalarField::typeName
    )
{
    porousMedium::gamma_ = gamma;
    porousMedium::volumetricAreaSS_ = volumetricAreaSS;
    porousMedium::volumetricAreaFuel_ = volumetricAreaFuel;
    porousMedium::volumetricAreaHX_ = volumetricAreaHX;
    porousMedium::rotate_ = rotate;
    porousMedium::externalHeatTransferCoefficient_ = externalHeatTransferCoefficient;
    porousMedium::externalT_ = externalT;
    porousMedium::externalVolHeatSource_ = externalVolHeatSource;
    porousMedium::externalRhoCp_ = externalRhoCp;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::correlationPorousMedium::~correlationPorousMedium()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp< Foam::volScalarField >
Foam::correlationPorousMedium::PrandtlNumber() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "correlationPorousMedium::PrandtlNumber",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mu_*thermo_.Cp()/thermo_.kappa()
        )
    );
}

Foam::tmp< Foam::volScalarField >
Foam::correlationPorousMedium::ReynoldsNumber() const //Re number
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "correlationPorousMedium::ReynoldsNumber",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            max((rho_*(hydraulicDiameter_ * mag(U_))/mu_),SMALL)
        )
    );
}

Foam::tmp< Foam::volVectorField >
Foam::correlationPorousMedium::darcyFrictionFactor() const
{
    tmp<volVectorField> tDarcyFrictionFactor
    (
        new volVectorField
        (
            IOobject
            (
                "correlationPorousMedium::darcyFrictionFactor",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("", dimless, vector::zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    tmp<volScalarField> tRe(new volScalarField(ReynoldsNumber()));
    tmp<volVectorField> tUR(new volVectorField(rotate_ & U_));
    tmp<volScalarField> tmagU(new volScalarField(mag(U_)));

    volVectorField& darcyFrictionFactor = tDarcyFrictionFactor.ref();
    volScalarField& Re = tRe.ref();
    volVectorField& UR = tUR.ref();
    volScalarField& magU = tmagU.ref();

    //here below I multiply and divide by U and magU so that the inner product between heatTransferCoefficient and UR gives a weighted average
    tmp<volVectorField> tReynoldsTurbVect(new volVectorField(reynoldsTurb_));   
    volVectorField& reynoldsTurbVect = tReynoldsTurbVect.ref();
    reynoldsTurbVect.replace(0,reynoldsTurbVect.component(0)*UR.component(0)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(0))));
    reynoldsTurbVect.replace(1,reynoldsTurbVect.component(1)*UR.component(1)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(1))));
    reynoldsTurbVect.replace(2,reynoldsTurbVect.component(2)*UR.component(2)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(2))));
    tmp<volScalarField> treynoldsTurb = max(mag((reynoldsTurbVect & UR) / max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),magU)),SMALL);
    volScalarField& reynoldsTurb = treynoldsTurb.ref();

    tmp<volVectorField>  tReynoldsLamVect(new volVectorField(reynoldsLam_));
    volVectorField& reynoldsLamVect = tReynoldsLamVect.ref();
    reynoldsLamVect.replace(0,reynoldsLamVect.component(0)*UR.component(0)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(0))));
    reynoldsLamVect.replace(1,reynoldsLamVect.component(1)*UR.component(1)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(1))));
    reynoldsLamVect.replace(2,reynoldsLamVect.component(2)*UR.component(2)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(2))));
    tmp<volScalarField> treynoldsLam = max(mag((reynoldsLamVect & UR) / max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),magU)),SMALL);
    volScalarField& reynoldsLam = treynoldsLam.ref();

    tmp<volScalarField> darcyFrictionFactorTurbI(new volScalarField(darcyConstTurb_.component(0)*pow(max(reynoldsTurb,Re),darcyExpTurb_.component(0))));
    tmp<volScalarField> darcyFrictionFactorTurbJ(new volScalarField(darcyConstTurb_.component(1)*pow(max(reynoldsTurb,Re),darcyExpTurb_.component(1))));
    tmp<volScalarField> darcyFrictionFactorTurbK(new volScalarField(darcyConstTurb_.component(2)*pow(max(reynoldsTurb,Re),darcyExpTurb_.component(2))));

    tmp<volScalarField> darcyFrictionFactorLamI(new volScalarField(darcyConstLam_.component(0)*pow(min(reynoldsLam,Re),darcyExpLam_.component(0))));
    tmp<volScalarField> darcyFrictionFactorLamJ(new volScalarField(darcyConstLam_.component(1)*pow(min(reynoldsLam,Re),darcyExpLam_.component(1))));
    tmp<volScalarField> darcyFrictionFactorLamK(new volScalarField(darcyConstLam_.component(2)*pow(min(reynoldsLam,Re),darcyExpLam_.component(2))));

    tmp<volScalarField> tReCoeff(new volScalarField((Re - reynoldsLam) / max(SMALL,(reynoldsTurb - reynoldsLam))));
    volScalarField& ReCoeff = tReCoeff.ref();

    tmp<volScalarField> darcyFrictionFactorI(new volScalarField(darcyFrictionFactorLamI() +  max(0.0,min(1.0,ReCoeff)) * (darcyFrictionFactorTurbI - darcyFrictionFactorLamI())));
    tmp<volScalarField> darcyFrictionFactorJ(new volScalarField(darcyFrictionFactorLamJ() +  max(0.0,min(1.0,ReCoeff)) * (darcyFrictionFactorTurbJ - darcyFrictionFactorLamJ())));
    tmp<volScalarField> darcyFrictionFactorK(new volScalarField(darcyFrictionFactorLamK() +  max(0.0,min(1.0,ReCoeff)) * (darcyFrictionFactorTurbK - darcyFrictionFactorLamK())));

    darcyFrictionFactor.replace(0,darcyFrictionFactorI);
    darcyFrictionFactor.replace(1,darcyFrictionFactorJ);
    darcyFrictionFactor.replace(2,darcyFrictionFactorK);

    return tDarcyFrictionFactor;
}

Foam::tmp< Foam::volTensorField >
Foam::correlationPorousMedium::dragCoeff() const
{
    tmp<volTensorField> tDragCoeff
    (
        new volTensorField
        (
            IOobject
            (
                "correlationPorousMedium::dragCoeff",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor("", dimensionSet(1,-3,-1,0,0,0,0), tensor::zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    tmp<volVectorField> tUR(new volVectorField(rotate_ & U_));
    tmp<volVectorField> tdf(new volVectorField(darcyFrictionFactor()));

    volTensorField& dragCoeff = tDragCoeff.ref();
    volVectorField& UR = tUR.ref();
    volVectorField& df = tdf.ref();

    dragCoeff.replace(0,0.5* rho_ * mag(UR)*df.component(0)/hydraulicDiameter_);
    dragCoeff.replace(4,0.5* rho_ * mag(UR)*df.component(1)/hydraulicDiameter_);
    dragCoeff.replace(8,0.5* rho_ * mag(UR)*df.component(2)/hydraulicDiameter_);

    tmp<volTensorField> rotateBack(new volTensorField(rotate_.T()));

    dragCoeff = rotateBack & dragCoeff & rotate_;

    return tDragCoeff;
}

Foam::tmp< Foam::volVectorField >
Foam::correlationPorousMedium::pump() const
{
    return tmp<volVectorField>
    (
        new volVectorField
        (
            IOobject
            (
                "correlationPorousMedium::pump",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pumpMomentumSource()
        )
    );
}

Foam::tmp< Foam::volScalarField >
Foam::correlationPorousMedium::equilibriumK() const
{
    tmp<volVectorField> tEquilibriumK
    (
        new volVectorField
        (
            IOobject
            (
                "correlationPorousMedium::equilibriumK",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("", dimensionSet(0,2,-2,0,0,0,0), vector::zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    tmp<volScalarField> tRe(new volScalarField(ReynoldsNumber()));
    tmp<volVectorField> tUR(new volVectorField(rotate_ & U_));
    tmp<volScalarField> tmagU(new volScalarField(mag(U_)));

    volVectorField& equilibriumK = tEquilibriumK.ref();
    volScalarField& Re = tRe.ref();
    volVectorField& UR = tUR.ref();
    volScalarField& magU = tmagU.ref();

    // SR - Compute turbulence intensity separately and limit it by the laminar-turbulence factor fRe
    updateFRe();
    //Info << "porousMedium::fRe min/max/avg from k " << gMin(fRe_) << "/" << gAverage(fRe_) << "/" << gMax(fRe_) << endl;
    turbulenceIntensity_.replace(0, fRe_*turbulenceIntensityConst_.component(0)*pow(max(SMALL,Re),turbulenceIntensityExp_.component(0)));
    turbulenceIntensity_.replace(1, fRe_*turbulenceIntensityConst_.component(1)*pow(max(SMALL,Re),turbulenceIntensityExp_.component(1)));
    turbulenceIntensity_.replace(2, fRe_*turbulenceIntensityConst_.component(2)*pow(max(SMALL,Re),turbulenceIntensityExp_.component(2)));

    equilibriumK.replace(0,1.5* pow(magU*turbulenceIntensity_.component(0), 2.0));
    equilibriumK.replace(1,1.5* pow(magU*turbulenceIntensity_.component(1), 2.0));
    equilibriumK.replace(2,1.5* pow(magU*turbulenceIntensity_.component(2), 2.0));

    equilibriumK.replace(0,equilibriumK.component(0)*UR.component(0)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(0))));
    equilibriumK.replace(1,equilibriumK.component(1)*UR.component(1)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(1))));
    equilibriumK.replace(2,equilibriumK.component(2)*UR.component(2)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(2))));

    tmp<volScalarField> scalarEquilibriumK
    (
        new volScalarField
        (
            max(dimensionedScalar("", dimensionSet(0,3,-3,0,0,0,0),SMALL),(equilibriumK & UR))
            / 
            max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),magU)
        )
    );

    return scalarEquilibriumK;
}

Foam::tmp< Foam::volScalarField >
Foam::correlationPorousMedium::equilibriumEpsilon() const
{
    tmp<volVectorField> tEquilibriumK
    (
        new volVectorField
        (
            IOobject
            (
                "correlationPorousMedium::equilibriumK",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("", dimensionSet(0,2,-2,0,0,0,0), vector::zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tmp<volVectorField> tEquilibriumEpsilon
    (
        new volVectorField
        (
            IOobject
            (
                "correlationPorousMedium::equilibriumEpsilon",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("", dimensionSet(0,2,-3,0,0,0,0), vector::zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tmp<volScalarField> tRe(new volScalarField(ReynoldsNumber()));
    tmp<volVectorField> tUR(new volVectorField(rotate_ & U_));
    tmp<volScalarField> tmagU(new volScalarField(mag(U_)));

    volVectorField& equilibriumK = tEquilibriumK.ref();
    volVectorField& equilibriumEpsilon = tEquilibriumEpsilon.ref();
    volScalarField& Re = tRe.ref();
    volVectorField& UR = tUR.ref();
    volScalarField& magU = tmagU.ref();

    // SR - Compute turbulence intensity separately and limit it by the laminar-turbulence factor fRe
    updateFRe();
    //Info << "porousMedium::fRe min/max/avg from epsilon " << gMin(fRe_) << "/" << gAverage(fRe_) << "/" << gMax(fRe_) << endl;
    turbulenceIntensity_.replace(0, fRe_*turbulenceIntensityConst_.component(0)*pow(max(SMALL,Re),turbulenceIntensityExp_.component(0)));
    turbulenceIntensity_.replace(1, fRe_*turbulenceIntensityConst_.component(1)*pow(max(SMALL,Re),turbulenceIntensityExp_.component(1)));
    turbulenceIntensity_.replace(2, fRe_*turbulenceIntensityConst_.component(2)*pow(max(SMALL,Re),turbulenceIntensityExp_.component(2)));

    equilibriumK.replace(0,1.5* pow(magU*turbulenceIntensity_.component(0), 2.0));
    equilibriumK.replace(1,1.5* pow(magU*turbulenceIntensity_.component(1), 2.0));
    equilibriumK.replace(2,1.5* pow(magU*turbulenceIntensity_.component(2), 2.0));

    equilibriumEpsilon.replace(0,pow(0.09,3.0/4.0)*pow(mag(equilibriumK.component(0)),3.0/2.0)/((turbulenceLengthScaleConst_.component(0)+SMALL)*hydraulicDiameter_));
    equilibriumEpsilon.replace(1,pow(0.09,3.0/4.0)*pow(mag(equilibriumK.component(1)),3.0/2.0)/((turbulenceLengthScaleConst_.component(1)+SMALL)*hydraulicDiameter_));
    equilibriumEpsilon.replace(2,pow(0.09,3.0/4.0)*pow(mag(equilibriumK.component(2)),3.0/2.0)/((turbulenceLengthScaleConst_.component(2)+SMALL)*hydraulicDiameter_));

    equilibriumEpsilon.replace(0,equilibriumEpsilon.component(0)*UR.component(0)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(0))));
    equilibriumEpsilon.replace(1,equilibriumEpsilon.component(1)*UR.component(1)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(1))));
    equilibriumEpsilon.replace(2,equilibriumEpsilon.component(2)*UR.component(2)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(2))));

    tmp<volScalarField> scalarEquilibriumEpsilon(new volScalarField(
        max(dimensionedScalar("", dimensionSet(0,3,-4,0,0,0,0),SMALL),(equilibriumEpsilon & UR))
        / max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),magU)));


    return scalarEquilibriumEpsilon;
}

// SR
void Foam::correlationPorousMedium::updateFRe() const
{
    volScalarField Re(ReynoldsNumber());
    vector highReVector(max(reynoldsTurb()).value());
    vector lowReVector(max(reynoldsLam()).value());
    scalar highRe(max(max(highReVector[0], highReVector[1]), highReVector[2]));
    scalar lowRe(min(min(lowReVector[0], lowReVector[1]), lowReVector[2]));

    scalar m(1/(max(highRe-lowRe, SMALL)));
    scalar q(-m*lowRe);

    forAll(Re, i)
    {
        scalar unboundFRe(m*Re[i]+q);
        fRe_[i] = max(scalar(0), min(scalar(1), unboundFRe));
    }
}

Foam::tmp< Foam::volScalarField >
Foam::correlationPorousMedium::heatTransferCoefficient() const
{
    tmp<volVectorField> tHeatTransferCoefficient
    (
        new volVectorField
        (
            IOobject
            (
                "correlationPorousMedium::heatTransferCoefficient",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("", dimensionSet(1,0,-3,-1,0,0,0), vector::zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    tmp<volScalarField> tPr(new volScalarField(PrandtlNumber()));
    tmp<volScalarField> tRe(new volScalarField(ReynoldsNumber()));
    tmp<volVectorField> tUR(new volVectorField(rotate_ & U_));

    volVectorField& heatTransferCoefficient = tHeatTransferCoefficient.ref();
    volScalarField& Pr = tPr.ref();
    volScalarField& Re = tRe.ref();
    volVectorField& UR = tUR.ref();

    //here below I multiply and divide by U and magU only to make it so that the inner product between heatTransferCoefficient and UR gives a weighted average
    tmp<volVectorField> tReynoldsTurbVect(new volVectorField(reynoldsTurb_));   
    volVectorField& reynoldsTurbVect = tReynoldsTurbVect.ref();
    reynoldsTurbVect.replace(0,reynoldsTurbVect.component(0)*UR.component(0)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(0))));
    reynoldsTurbVect.replace(1,reynoldsTurbVect.component(1)*UR.component(1)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(1))));
    reynoldsTurbVect.replace(2,reynoldsTurbVect.component(2)*UR.component(2)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(2))));
    tmp<volScalarField> tReynoldsTurb(new volScalarField(max(mag((reynoldsTurbVect & UR) / max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR))),SMALL)));
    volScalarField& reynoldsTurb = tReynoldsTurb.ref();

    tmp<volVectorField> tReynoldsLamVect(new volVectorField(reynoldsLam_));
    volVectorField& reynoldsLamVect = tReynoldsLamVect.ref();
    reynoldsLamVect.replace(0,reynoldsLamVect.component(0)*UR.component(0)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(0))));
    reynoldsLamVect.replace(1,reynoldsLamVect.component(1)*UR.component(1)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(1))));
    reynoldsLamVect.replace(2,reynoldsLamVect.component(2)*UR.component(2)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(2))));
    tmp<volScalarField> tReynoldsLam(new volScalarField(max(mag((reynoldsLamVect & UR) / max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR))),SMALL)));
    volScalarField& reynoldsLam = tReynoldsLam.ref();

    tmp<volScalarField> nusseltTurbI(new volScalarField(nusseltConstTurb1_.component(0)*pow(max(reynoldsTurb,Re),nusseltReExpTurb_.component(0))*pow(Pr,nusseltPrExpTurb_.component(0))+nusseltConstTurb2_.component(0)));
    tmp<volScalarField> nusseltTurbJ(new volScalarField(nusseltConstTurb1_.component(1)*pow(max(reynoldsTurb,Re),nusseltReExpTurb_.component(1))*pow(Pr,nusseltPrExpTurb_.component(1))+nusseltConstTurb2_.component(1)));
    tmp<volScalarField> nusseltTurbK(new volScalarField(nusseltConstTurb1_.component(2)*pow(max(reynoldsTurb,Re),nusseltReExpTurb_.component(2))*pow(Pr,nusseltPrExpTurb_.component(2))+nusseltConstTurb2_.component(2)));

    tmp<volScalarField> tNusseltLamI(new volScalarField(nusseltConstLam1_.component(0)*pow(min(reynoldsLam,Re),nusseltReExpLam_.component(0))*pow(Pr,nusseltPrExpLam_.component(0))+nusseltConstLam2_.component(0)));
    tmp<volScalarField> tNusseltLamJ(new volScalarField(nusseltConstLam1_.component(1)*pow(min(reynoldsLam,Re),nusseltReExpLam_.component(1))*pow(Pr,nusseltPrExpLam_.component(1))+nusseltConstLam2_.component(1)));
    tmp<volScalarField> tNusseltLamK(new volScalarField(nusseltConstLam1_.component(2)*pow(min(reynoldsLam,Re),nusseltReExpLam_.component(2))*pow(Pr,nusseltPrExpLam_.component(2))+nusseltConstLam2_.component(2)));
    volScalarField& nusseltLamI = tNusseltLamI.ref();
    volScalarField& nusseltLamJ = tNusseltLamJ.ref();
    volScalarField& nusseltLamK = tNusseltLamK.ref();

    tmp<volScalarField> tReCoeff(new volScalarField((Re - reynoldsLam) / max(SMALL,(reynoldsTurb - reynoldsLam))));
    volScalarField& ReCoeff = tReCoeff.ref();

    tmp<volScalarField> heatTransferCoefficientI(new volScalarField((nusseltLamI +  max(0.0,min(1.0,ReCoeff)) * (nusseltTurbI - nusseltLamI))*thermo_.kappa()/hydraulicDiameter_));
    tmp<volScalarField> heatTransferCoefficientJ(new volScalarField((nusseltLamJ +  max(0.0,min(1.0,ReCoeff)) * (nusseltTurbJ - nusseltLamJ))*thermo_.kappa()/hydraulicDiameter_));
    tmp<volScalarField> heatTransferCoefficientK(new volScalarField((nusseltLamK +  max(0.0,min(1.0,ReCoeff)) * (nusseltTurbK - nusseltLamK))*thermo_.kappa()/hydraulicDiameter_));
    //here below I multiply and divide by U and magU only to make it so that the inner product betwee heatTransferCoefficient and UR gives a weighted average
    heatTransferCoefficient.replace(0,heatTransferCoefficientI*UR.component(0)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(0))));
    heatTransferCoefficient.replace(1,heatTransferCoefficientJ*UR.component(1)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(1))));
    heatTransferCoefficient.replace(2,heatTransferCoefficientK*UR.component(2)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(2))));

    tmp<volScalarField> scalarHeatTransferCoefficient(new volScalarField(
        max(dimensionedScalar("", dimensionSet(1,1,-4,-1,0,0,0),SMALL),(heatTransferCoefficient & UR))
        / max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR))));

    return scalarHeatTransferCoefficient;
}


Foam::tmp< Foam::volScalarField >
Foam::correlationPorousMedium::kepsilonConvergenceRate() const
{
    tmp<volVectorField> tKepsilonConvergenceRate
    (
        new volVectorField
        (
            IOobject
            (
                "correlationPorousMedium::kepsilonConvergenceRate",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("", dimensionSet(0,-1,0,0,0,0,0), vector::zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    tmp<volVectorField> tUR = rotate_ & U_;

    volVectorField& kepsilonConvergenceRate = tKepsilonConvergenceRate.ref();
    volVectorField& UR = tUR.ref();

    kepsilonConvergenceRate.replace(0,kepsilonConvergenceRate_.component(0)*UR.component(0)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(0))));
    kepsilonConvergenceRate.replace(1,kepsilonConvergenceRate_.component(1)*UR.component(1)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(1))));
    kepsilonConvergenceRate.replace(2,kepsilonConvergenceRate_.component(2)*UR.component(2)/max(dimensionedScalar("", dimensionSet(0,1,-1,0,0,0,0),SMALL),mag(UR.component(2))));

    tmp<volScalarField> scalarKepsilonConvergenceRate(new volScalarField(max(dimensionedScalar("", dimensionSet(0,0,-1,0,0,0,0),SMALL),kepsilonConvergenceRate & UR)));

    return scalarKepsilonConvergenceRate;
}

Foam::tmp<Foam::volScalarField>
Foam::correlationPorousMedium::fRe() const
{
    tmp<volScalarField> tfRe
    (
        new volScalarField(fRe_)
    );
    return tfRe;
}

// ************************************************************************* //

