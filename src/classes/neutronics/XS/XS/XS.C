/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2412                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2024 OpenCFD Ltd.          |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "XS.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(XS, 0);
    defineRunTimeSelectionTable(XS, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XS::XS
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    nuclearData_
    (
        IOobject
        (
            "nuclearData",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    isLowMemory_(nuclearData_.lookupOrDefault("isLowMemory", false)),
    energyGroups_(nuclearData_.lookupOrDefault("energyGroups", 1)),
    precGroups_(nuclearData_.lookupOrDefault("precGroups", 1)),
    legendreMoments_(1+nuclearData_.lookupOrDefault("legendreMoments", 0)),
    axialOrientation_(nuclearData_.lookupOrDefault("axialOrientation", vector(0.0, 0.0, 1.0))),
    ScNo_(nuclearData_.lookupOrDefault("ScNo", 1.0)),
    polyharmonicSplineMode_(nuclearData_.lookupOrDefault("polyharmonicSplineMode", 1)),
    xsVariablesDict_(nuclearData_.subDict("xsVariables")),
    nxsVariables_(xsVariablesDict_.toc().size()),
    xsVariableNames_(nxsVariables_),
    xsVariableTypes_(nxsVariables_),
    xsVariableFields_(nxsVariables_),
    xsVariableTransformedFields_(nxsVariables_),
    states_(nuclearData_.lookup("states")),
    referenceState_(states_.first().dict()),
    referenceZones_(referenceState_.lookup("zones")),
    zoneNumber_(referenceZones_.size()),
    IV_(energyGroups_),
    D_(energyGroups_),
    nuSigmaEff_(energyGroups_),
    sigmaPow_(energyGroups_),
    sigmaRemoval_(energyGroups_),
    sigmaFromTo_(legendreMoments_),
    chiPrompt_(energyGroups_),
    chiDelayed_(energyGroups_),
    Beta_(precGroups_),
    BetaTot_
    (
        IOobject
        (
            "BetaTot",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    lambda_(precGroups_),
    fuelFraction_
    (
        IOobject
        (
            "fuelFraction",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    secondaryPowerVolumeFraction_
    (
        IOobject
        (
            "secondaryPowerVolumeFraction",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    fractionToSecondaryPower_
    (
        IOobject
        (
            "fractionToSecondaryPower",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    discFactor_(energyGroups_),
    sigmaFromToYesNo_(energyGroups_),
    fuelFractionList_(zoneNumber_),
    secondaryPowerVolumeFractionList_(zoneNumber_),
    fractionToSecondaryPowerList_(zoneNumber_),
    dfAdjustList_(zoneNumber_),
    discFactorList_(zoneNumber_),
    integralFluxList_(zoneNumber_),
    adjustDiscFactors_(nuclearData_.lookupOrDefault("adjustDiscFactors", false)),
    useGivenDiscFactors_(nuclearData_.lookupOrDefault("useGivenDiscFactors", false)),
    groupsWoDF_(nuclearData_.lookupOrDefault<List<int>>("groupsWoDF", List<int>())),
    fastNeutrons_(nuclearData_.lookupOrDefault("fastNeutrons", true)),
    doNotParametrize_(nuclearData_.lookupOrDefault<List<int>>("doNotParametrize", List<int>())),
    disp_
    (
        IOobject
        (
            "disp",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("d_zero", dimLength, vector::zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    // radExp_
    // (
    //     IOobject
    //     (
    //         "radExp",
    //         mesh.time().timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
    //     zeroGradientFvPatchScalarField::typeName
    // ),
    // axExp_
    // (
    //     IOobject
    //     (
    //         "axExp",
    //         mesh.time().timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
    //     zeroGradientFvPatchScalarField::typeName
    // ),
    // logT_
    // (
    //     IOobject
    //     (
    //         "logT",
    //         mesh.time().timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
    //     zeroGradientFvPatchScalarField::typeName
    // ),
    // sqrtT_
    // (
    //     IOobject
    //     (
    //         "diffT",
    //         mesh.time().timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar("", dimensionSet(0,0,0,0.5,0,0,0), 0.0),
    //     zeroGradientFvPatchScalarField::typeName
    // ),
    IVList_(zoneNumber_),
    chiPromptList_(zoneNumber_),
    chiDelayedList_(zoneNumber_),
    BetaList_(zoneNumber_),
    BetaTotList_(zoneNumber_),
    lambdaList_(zoneNumber_),
    DList_(zoneNumber_),
    nuSigmaEffList_(zoneNumber_),
    sigmaPowList_(zoneNumber_),
    sigmaRemovalList_(zoneNumber_),
    sigmaFromToList_(zoneNumber_),
    CRmove_
    (
        IOobject
        (
            "CRmove",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    CRentries_(CRmove_.lookupOrDefault("zones", PtrList<entry>(0))),
    CRNumber_(CRentries_.size()),
    CRstart_(CRNumber_),
    CRfinish_(CRNumber_),
    CRspeed_(CRNumber_),
    CRFollowerName_(CRNumber_),
    CRinitialPosition_(CRNumber_),
    CRposition_(CRNumber_),
    initialDistanceFromMeshCR_(CRNumber_)
{
    #include "readNuclearData.H"
    #include "createXSfields.H"
    init();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XS::~XS()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::XS::correct
(
    // const volScalarField& Tfuel,
    // const volScalarField& Tclad,
    // const volScalarField& rhoCool,
    // const volScalarField& Tcool,
    // const volVectorField& Disp,
    // const volScalarField& TStructMech
)
{
    #include "setNeutronicsVariables.H"
}


void Foam::XS::init()
{
    // Build interpolation scheme
    forAll(referenceZones_, zoneI)
    {
        forAll(DList_[zoneI], energyI)
        {
            DList_[zoneI][energyI].build();
            nuSigmaEffList_[zoneI][energyI].build();
            sigmaPowList_[zoneI][energyI].build();
            sigmaRemovalList_[zoneI][energyI].build();
        }
        forAll(sigmaFromToList_[zoneI], momentI)
        {
            forAll(sigmaFromToList_[zoneI][momentI], energyJ)
            {
                forAll(sigmaFromToList_[zoneI][momentI][energyJ], energyI)
                {
                    sigmaFromToList_[zoneI][momentI][energyJ][energyI].build();
                }
            }
        }
    }

    #include "setNeutronicsConstants.H"
    #include "setPrecConst.H"
}


void Foam::XS::adjustDiscFactors(const PtrList<volScalarField>& fluxStar)
{
    if (adjustDiscFactors_)
    {
        #include "adjustDiscFactors.H"
    }
}


Foam::tmp<Foam::volScalarField> Foam::XS::sigmaFromTo
(
    label momentI,
    label energyJ,
    label energyI//,
    // const volScalarField& Tfuel,
    // const volScalarField& Tclad,
    // const volScalarField& rhoCool,
    // const volScalarField& Tcool,
    // const volVectorField& Disp,
    // const volScalarField& TStructMech
)
{
    if (!isLowMemory_)
    {
        return(sigmaFromTo_[momentI][energyJ][energyI]);
    }

    tmp<volScalarField> tsigmaFromTo
    (
        new volScalarField
        (
            IOobject
            (
                "XSLowMem::sigmaFromTo",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& sigmaFromTo(tsigmaFromTo.ref());

    // Check if XS is parametrized for the current energy group
    bool isParametrize = true;
    forAll(doNotParametrize_, groupI)
    {
        if (energyI == doNotParametrize_[groupI])
        {
            isParametrize = false;
            break;
        }
    }

    // Update cells for all zones
    scalarList variableValues(nxsVariables_);
    forAll(referenceZones_, zoneI)
    {
        label zone = zoneI;

        const word& name = referenceZones_[zoneI].keyword();

        label zoneId = mesh_.cellZones().findZoneID(name);

        forAll(mesh_.cellZones()[zoneId], cellIlocal)
        {
            label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];
            zone = zoneI;

            forAll(variableValues, varI)
            {
                if (xsVariableTypes_[varI] == "lin")
                {
                    variableValues[varI] = xsVariableFields_[varI][cellIglobal];
                }
                else
                {
                    variableValues[varI] = xsVariableTransformedFields_[varI][cellIglobal];
                }
            }

            sigmaFromTo[cellIglobal] =
                sigmaFromToList_[zone][momentI][energyJ][energyI].get
                (
                    variableValues,
                    isParametrize
                );
        }
    }

    return tsigmaFromTo;
}

// ************************************************************************* //
