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

#include "byZoneCorrelationPorousMedium.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(byZoneCorrelationPorousMedium, 0);

    addToRunTimeSelectionTable
    (
        porousMedium,
        byZoneCorrelationPorousMedium,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::byZoneCorrelationPorousMedium::byZoneCorrelationPorousMedium
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
            FatalIOErrorIn("byZoneCorrelationPorousMedium::"
                "byZoneCorrelationPorousMedium(...)", *this)
                << "Zone " << name  << " not found on mesh " << mesh_.name()
                << abort(FatalIOError);
        }

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
    reynoldsTurb_.correctBoundaryConditions();
    reynoldsLam_.correctBoundaryConditions();
    darcyConstTurb_.correctBoundaryConditions();
    darcyConstLam_.correctBoundaryConditions();
    darcyExpTurb_.correctBoundaryConditions();
    darcyExpLam_.correctBoundaryConditions();
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

Foam::byZoneCorrelationPorousMedium::~byZoneCorrelationPorousMedium()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
