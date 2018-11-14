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

#include "byZoneSodiumCorrelationPorousMedium.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(byZoneSodiumCorrelationPorousMedium, 0);

    addToRunTimeSelectionTable
    (
        porousMedium,
        byZoneSodiumCorrelationPorousMedium,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::byZoneSodiumCorrelationPorousMedium::byZoneSodiumCorrelationPorousMedium
(
    const volScalarField& rho,
    const volVectorField& U,
    const rhoThermo& thermo
)
:
    correlationPorousMedium
    (
        rho,
        U,
        thermo,
        dimensionedScalar("", dimless, 1.0),
        dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 1.0),
        dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 1.0),
        dimensionedTensor("", dimless, tensor::I),
        dimensionedScalar("", dimensionSet(1,0,-3,-1,0,0,0), SMALL),
        dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0), 0.0),
        dimensionedScalar("", dimensionSet(1,-1,-3,0,0,0,0), 0.0),
        dimensionedScalar("", dimensionSet(1,-1,-2,-1,0,0,0), 1000000.0),
        dimensionedVector("", dimless, vector::zero),
        dimensionedVector("", dimless, vector::zero),
        dimensionedVector("", dimless, vector::zero),
        dimensionedVector("", dimless, vector::zero),
        dimensionedVector("", dimless, vector::zero),
        dimensionedVector("", dimless, vector::zero),
        dimensionedVector("", dimless, vector::zero),
        dimensionedVector("", dimless, vector::zero),
        dimensionedVector("", dimless, vector::zero),
        dimensionedVector("", dimless, vector::zero),
        dimensionedVector("", dimless, vector::zero),
        dimensionedVector("", dimless, vector::zero),
        dimensionedVector("", dimless, vector::zero),
        dimensionedVector("", dimless, vector::zero),
        dimensionedScalar("", dimLength, 1.0),
        dimensionedScalar("", dimLength, 10.0),
        dimensionedVector("", dimensionSet(1,-2,-2,0,0,0,0), vector::zero),
        dimensionedVector("", dimless, vector::zero),
        dimensionedVector("", dimless, vector::zero),
        dimensionedVector("", dimless, vector::zero),
        dimensionedVector("", dimensionSet(0,-1,0,0,0,0,0), vector::zero)
    ),
    relRough_
    (
        IOobject(
            "porousMedium::relativeRoughness",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(IOdictionary::lookupOrDefault("relativeRoughness", 0.0)),
        zeroGradientFvPatchScalarField::typeName
    ),
    Np_
    (
        IOobject(
            "porousMedium::numerOfPinsInAssembly",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(IOdictionary::lookupOrDefault("NPins", 1)),
        zeroGradientFvPatchScalarField::typeName
    ),
    Dp_
    (
        IOobject(
            "porousMedium::pinDiameter",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(IOdictionary::lookupOrDefault("pinDiameter", 1e-2)),
        zeroGradientFvPatchScalarField::typeName
    ),
    Dw_
    (
        IOobject(
            "porousMedium::wireDiameter",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(IOdictionary::lookupOrDefault("wireDiameter", 1e-3)),
        zeroGradientFvPatchScalarField::typeName
    ),
    PWetWall_
    (
        IOobject(
            "porousMedium::wettedWrapperPerimeter",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(IOdictionary::lookupOrDefault("wettedWrapperPerimeter", 1)),
        zeroGradientFvPatchScalarField::typeName
    ),
    Lw_
    (
        IOobject(
            "porousMedium::wireLeadLength",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(IOdictionary::lookupOrDefault("wireLeadLength", 1e-1)),
        zeroGradientFvPatchScalarField::typeName
    ),
    Pp_
    (
        IOobject(
            "porousMedium::pinPitch",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(IOdictionary::lookupOrDefault("pinPitch", 1e-1)),
        zeroGradientFvPatchScalarField::typeName
    )
{
    PtrList<entry> entries(IOdictionary::lookup("zones"));

    forAll(entries, i)
    {
        const word& name = entries[i].keyword();
        dictionary& dict = entries[i].dict();

        label zoneId = mesh_.cellZones().findZoneID(name);

        if (zoneId < 0)
        {
            FatalIOErrorIn("byZoneSodiumCorrelationPorousMedium::"
                "byZoneSodiumCorrelationPorousMedium(...)", *this)
                << "Zone " << name  << " not found on mesh " << mesh_.name()
                << abort(FatalIOError);
        }

        frictionModels_.insert
        (
            name,
            word(dict.lookup("frictionModel"))
        );

        const UList<label>& addr = mesh_.cellZones()[zoneId];

        IndirectList<scalar>(gamma_.internalField(), addr)
            = scalar(dict.lookupOrDefault("voidFraction",1.0));

        IndirectList<scalar>(volumetricAreaSS_.internalField(), addr)
            = scalar(dict.lookupOrDefault("volumetricAreaSS",1.0));

        IndirectList<scalar>(volumetricAreaFuel_.internalField(), addr)
            = scalar(dict.lookupOrDefault("volumetricAreaFuel",1.0));

        IndirectList<tensor>(rotate_.internalField(), addr)
            = tensor(coordinateSystem("",vector::zero,dict.lookupOrDefault("localZaxis",vector(0,0,1)),dict.lookupOrDefault("localXaxis",vector(1,0,0))).R().R().T());

        IndirectList<scalar>(externalHeatTransferCoefficient_.internalField(), addr)
            = scalar(dict.lookupOrDefault("externalHeatTransferCoefficient",0.0)+SMALL);

        IndirectList<scalar>(externalT_.internalField(), addr)
            = scalar(dict.lookupOrDefault("externalT",0.0));

        IndirectList<scalar>(externalVolHeatSource_.internalField(), addr)
            = scalar(dict.lookupOrDefault("externalVolHeatSource",0.0));

        IndirectList<scalar>(externalRhoCp_.internalField(), addr)
            = scalar(dict.lookupOrDefault("externalRhoCp",1000000.0));

        if (frictionModels_[name] == "Darcy")
        {
            IndirectList<vector>(reynoldsTurb_.internalField(), addr)
            = vector(dict.lookupOrDefault("reynoldsTurb",vector(2.3e3,2.3e3,2.3e3)));

            IndirectList<vector>(reynoldsLam_.internalField(), addr)
            = vector(dict.lookupOrDefault("reynoldsLam",vector(1e3,1e3,1e3)));

            IndirectList<vector>(darcyConstTurb_.internalField(), addr)
            = vector(dict.lookupOrDefault("darcyConstTurb",vector::zero));

            IndirectList<vector>(darcyConstLam_.internalField(), addr)
                = vector(dict.lookupOrDefault("darcyConstLam",vector::zero));

            IndirectList<vector>(darcyExpTurb_.internalField(), addr)
                = vector(dict.lookupOrDefault("darcyExpTurb",vector::zero));

            IndirectList<vector>(darcyExpLam_.internalField(), addr)
                = vector(dict.lookupOrDefault("darcyExpLam",vector::zero));

            reynoldsTurb_.correctBoundaryConditions();
            reynoldsLam_.correctBoundaryConditions();
            darcyConstTurb_.correctBoundaryConditions();
            darcyConstLam_.correctBoundaryConditions();
            darcyExpTurb_.correctBoundaryConditions();
            darcyExpLam_.correctBoundaryConditions();
        }

        else if (frictionModels_[name] == "Churchill")
        {
            IndirectList<scalar>(relRough_.internalField(), addr)
            = scalar(readScalar(dict.lookup("relativeRoughness")));

            relRough_.correctBoundaryConditions();
        }

        else if (frictionModels_[name] == "Rehme")
        {
            IndirectList<scalar>(Np_.internalField(), addr)
            = scalar(readScalar(dict.lookup("NPins")));

            IndirectList<scalar>(Dp_.internalField(), addr)
            = scalar(readScalar(dict.lookup("pinDiameter")));

            IndirectList<scalar>(Dw_.internalField(), addr)
            = scalar(readScalar(dict.lookup("wireDiameter")));

            IndirectList<scalar>(PWetWall_.internalField(), addr)
            = scalar(readScalar(dict.lookup("wettedWrapperPerimeter")));

            IndirectList<scalar>(Lw_.internalField(), addr)
            = scalar(readScalar(dict.lookup("wireLeadLength")));

            IndirectList<scalar>(Pp_.internalField(), addr)
            = scalar(readScalar(dict.lookup("pinPitch")));

            Np_.correctBoundaryConditions();
            Dp_.correctBoundaryConditions();
            Dw_.correctBoundaryConditions();
            PWetWall_.correctBoundaryConditions();
            Lw_.correctBoundaryConditions();
            Pp_.correctBoundaryConditions();
        }

        IndirectList<vector>(nusseltConstTurb1_.internalField(), addr)
            = vector(dict.lookupOrDefault("nusseltConstTurb1",vector::zero));

        IndirectList<vector>(nusseltConstLam1_.internalField(), addr)
            = vector(dict.lookupOrDefault("nusseltConstLam1",vector::zero));

        IndirectList<vector>(nusseltConstTurb2_.internalField(), addr)
            = vector(dict.lookupOrDefault("nusseltConstTurb2",vector::zero));

        IndirectList<vector>(nusseltConstLam2_.internalField(), addr)
            = vector(dict.lookupOrDefault("nusseltConstLam2",vector::zero));

        IndirectList<vector>(nusseltReExpTurb_.internalField(), addr)
            = vector(dict.lookupOrDefault("nusseltReExpTurb",vector::zero));

        IndirectList<vector>(nusseltReExpLam_.internalField(), addr)
            = vector(dict.lookupOrDefault("nusseltReExpLam",vector::zero));

        IndirectList<vector>(nusseltPrExpTurb_.internalField(), addr)
            = vector(dict.lookupOrDefault("nusseltPrExpTurb",vector::zero));

        IndirectList<vector>(nusseltPrExpLam_.internalField(), addr)
            = vector(dict.lookupOrDefault("nusseltPrExpLam",vector::zero));

        IndirectList<scalar>(hydraulicDiameter_.internalField(), addr)
            = scalar(dict.lookupOrDefault("hydraulicDiameter",1.0));

        IndirectList<scalar>(hydraulicDiameterStructure_.internalField(), addr)
            = scalar(dict.lookupOrDefault("hydraulicDiameterStructure",10.0));

        IndirectList<vector>(pumpMomentumSource_.internalField(), addr)
            = vector(dict.lookupOrDefault("pumpMomentumSource",vector::zero));

        IndirectList<vector>(turbulenceIntensityConst_.internalField(), addr)
            = vector(dict.lookupOrDefault("turbulenceIntensityConst",vector::zero));

        IndirectList<vector>(turbulenceIntensityExp_.internalField(), addr)
            = vector(dict.lookupOrDefault("turbulenceIntensityExp",vector::zero));

        IndirectList<vector>(turbulenceLengthScaleConst_.internalField(), addr)
            = vector(dict.lookupOrDefault("turbulenceLengthScaleConst",vector::zero));

        IndirectList<vector>(kepsilonConvergenceRate_.internalField(), addr)
            = vector(dict.lookupOrDefault("kepsilonConvergenceRate",vector::zero));
    }

    gamma_.correctBoundaryConditions();
    volumetricAreaSS_.correctBoundaryConditions();
    volumetricAreaFuel_.correctBoundaryConditions();
    rotate_.correctBoundaryConditions();
    externalHeatTransferCoefficient_.correctBoundaryConditions();
    externalT_.correctBoundaryConditions();
    externalVolHeatSource_.correctBoundaryConditions();
    externalRhoCp_.correctBoundaryConditions();
    nusseltConstTurb1_.correctBoundaryConditions();
    nusseltConstLam1_.correctBoundaryConditions();
    nusseltConstTurb2_.correctBoundaryConditions();
    nusseltConstLam2_.correctBoundaryConditions();
    nusseltReExpTurb_.correctBoundaryConditions();
    nusseltReExpLam_.correctBoundaryConditions();
    nusseltPrExpTurb_.correctBoundaryConditions();
    nusseltPrExpLam_.correctBoundaryConditions();
    hydraulicDiameter_.correctBoundaryConditions();
    hydraulicDiameterStructure_.correctBoundaryConditions();
    pumpMomentumSource_.correctBoundaryConditions();
    turbulenceIntensityConst_.correctBoundaryConditions();
    turbulenceIntensityExp_.correctBoundaryConditions();
    turbulenceLengthScaleConst_.correctBoundaryConditions();
    kepsilonConvergenceRate_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::byZoneSodiumCorrelationPorousMedium::~byZoneSodiumCorrelationPorousMedium()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField >
Foam::byZoneSodiumCorrelationPorousMedium::rehmeFrictionFactor() const
{
    tmp<volVectorField> trehmeFrictionFactor
    (
        new volVectorField
        (
            IOobject
            (
                "correlationPorousMedium::rehmeFrictionFactor",
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
    volVectorField& rehmeFrictionFactor = trehmeFrictionFactor.ref();

    volScalarField Re(max(ReynoldsNumber(), SMALL));

    volScalarField F
    (
        pow(Pp_/Dp_, 0.5) +
        pow
        (
            (
                7.6*(Dp_+Dw_)*pow(Pp_/Dp_, 2)
                /
                Lw_
            ),
            2.16
        )
    );

    scalar pi(3.14159265359);
    volScalarField PWetBundle(Np_*pi*(Dp_+Dw_));

    volScalarField f
    (
        (PWetBundle/(PWetBundle+PWetWall_))
        *
        (
            pow(F, 0.5)*64/Re
            +
            pow(F, 0.9335)*0.0816/pow(Re, 0.133)
        )
    );

    rehmeFrictionFactor.replace(0, f);
    rehmeFrictionFactor.replace(1, f);
    rehmeFrictionFactor.replace(2, f);

    return trehmeFrictionFactor;
}

Foam::tmp<Foam::volVectorField >
Foam::byZoneSodiumCorrelationPorousMedium::churchillFrictionFactor() const
{
    tmp<volVectorField> tChurchillFrictionFactor
    (
        new volVectorField
        (
            IOobject
            (
                "correlationPorousMedium::churchillFrictionFactor",
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
    volVectorField& churchillFrictionFactor = tChurchillFrictionFactor.ref();

    volScalarField Re(max(ReynoldsNumber(), SMALL));

    volScalarField A
    (
        pow
        (
            -2.457*
            Foam::log
            (
                pow(7/Re, 0.9)+
                0.27*relRough_
            ), 
            16
        )
    );

    volScalarField B
    (
        pow
        (
            37530/Re,
            16
        )
    );

    volScalarField f
    (
        8*
        pow
        (
            pow(8/Re, 12) +
            scalar(1) / 
            (
                pow(A+B, 1.5)
            ),
            scalar(1)/12
        )
    );

    churchillFrictionFactor.replace(0, f);
    churchillFrictionFactor.replace(1, f);
    churchillFrictionFactor.replace(2, f);

    return tChurchillFrictionFactor;
}

Foam::tmp<Foam::volTensorField >
Foam::byZoneSodiumCorrelationPorousMedium::dragCoeffFromFF(tmp<volVectorField> tFF) const
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
    volTensorField& dragCoeff = tDragCoeff.ref();
    volVectorField& FF = tFF.ref();

    volVectorField UR(rotate_ & U_);
    dragCoeff.replace(0,0.5* rho_ * mag(UR)*FF.component(0)/hydraulicDiameter_);
    dragCoeff.replace(4,0.5* rho_ * mag(UR)*FF.component(1)/hydraulicDiameter_);
    dragCoeff.replace(8,0.5* rho_ * mag(UR)*FF.component(2)/hydraulicDiameter_);

    tmp<volTensorField> rotateBack(new volTensorField(rotate_.T()));

    dragCoeff = rotateBack & dragCoeff & rotate_;

    return tDragCoeff;
}

Foam::tmp<Foam::volTensorField >
Foam::byZoneSodiumCorrelationPorousMedium::dragCoeff() const
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
    volTensorField& dragCoeff = tDragCoeff.ref();
    
    const cellZoneMesh& cellZones(mesh_.cellZones());
    wordList cellZoneNames(cellZones.names());

    //Info << "Re min / avg / max : " << gMin(ReynoldsNumber()()) << " / " << gAverage(ReynoldsNumber()()) << " / "<< gMax(ReynoldsNumber()()) << endl;

    forAll(cellZones, i)
    {
        word cellZoneName = cellZones.names()[i];
        
        
        if      (frictionModels_[cellZoneName] == "Darcy")
        {
            volTensorField partialDragCoeff
            (
                dragCoeffFromFF
                (
                    darcyFrictionFactor()
                )
            );
            forAll(cellZones[i], j)
            {
                dragCoeff[cellZones[i][j]] = partialDragCoeff[cellZones[i][j]];
            }
        }
        else if (frictionModels_[cellZoneName] == "Churchill")
        {
            volTensorField partialDragCoeff
            (
                dragCoeffFromFF
                (
                    churchillFrictionFactor()
                )
            );
            forAll(cellZones[i], j)
            {
                dragCoeff[cellZones[i][j]] = partialDragCoeff[cellZones[i][j]];
            }
        }
        else if (frictionModels_[cellZoneName] == "Rehme")
        {
            volTensorField partialDragCoeff
            (
                dragCoeffFromFF
                (
                    rehmeFrictionFactor()
                )
            );
            forAll(cellZones[i], j)
            {
                dragCoeff[cellZones[i][j]] = partialDragCoeff[cellZones[i][j]];
            }
        }
    }

    dragCoeff.correctBoundaryConditions();

    return tDragCoeff;
}


// ************************************************************************* //
