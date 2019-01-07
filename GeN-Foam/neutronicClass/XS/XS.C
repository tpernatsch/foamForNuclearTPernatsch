/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    nuclearData
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
    nuclearDataRadialExp
    (
        IOobject
        (
            "nuclearDataRadialExp",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nuclearDataAxialExp
    (
        IOobject
        (
            "nuclearDataAxialExp",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nuclearDataFuelTemp
    (
        IOobject
        (
            "nuclearDataFuelTemp",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nuclearDataRhoCool
    (
        IOobject
        (
            "nuclearDataRhoCool",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nuclearDataTCool
    (
        IOobject
        (
            "nuclearDataTCool",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nuclearDataCladExp
    (
        IOobject
        (
            "nuclearDataCladExp",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    energyGroups_(nuclearData.lookupOrDefault("energyGroups",1)),
    precGroups_(nuclearData.lookupOrDefault("precGroups",1)),
    IV_(energyGroups_),
    D_(energyGroups_),
    nuSigmaEff_(energyGroups_),
    sigmaPow_(energyGroups_),
    sigmaDisapp_(energyGroups_),
    sigmaFromTo_(energyGroups_),
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
    discFactor_(energyGroups_),
    ScNo_(nuclearData.lookupOrDefault("ScNo",1.0)),
    zoneNumber(nuclearData.lookup("zones").size()),    
    fuelFractionList(zoneNumber),
    dfAdjustList(zoneNumber),
    discFactorList(zoneNumber),
    integralFluxList(zoneNumber),
    fastNeutrons(nuclearData.lookupOrDefault("fastNeutrons",true)),
    adjustDiscFactors_(nuclearData.lookupOrDefault("adjustDiscFactors", false)),
    useGivenDiscFactors(nuclearData.lookupOrDefault("useGivenDiscFactors", false)),
    groupsWoDF(nuclearData.lookupOrDefault<List<int> >("groupsWoDF", List<int>())),
    doNotParametrize(nuclearData.lookupOrDefault<List<int> >("doNotParametrize", List<int>())),
    entries(nuclearData.lookup("zones")),
    IVList(zoneNumber),
    chiPromptList(zoneNumber),
    chiDelayedList(zoneNumber),
    BetaList(zoneNumber),
    BetaTotList(zoneNumber),
    lambdaList(zoneNumber),
    DList(zoneNumber),
    nuSigmaEffList(zoneNumber),
    sigmaPowList(zoneNumber),
    sigmaDisappList(zoneNumber),
    sigmaFromToList(zoneNumber), 
    TfuelRef(nuclearDataFuelTemp.lookupOrDefault("TfuelRef",900.0)),
    TfuelPerturbed(nuclearDataFuelTemp.lookupOrDefault("TfuelPerturbed",1200.0)),
    fuelTempDList(zoneNumber),
    fuelTempNuSigmaEffList(zoneNumber),
    fuelTempSigmaPowList(zoneNumber),
    fuelTempSigmaDisappList(zoneNumber),
    fuelTempSigmaFromToList(zoneNumber),
    AxExp(nuclearDataAxialExp.lookupOrDefault("expansionFromNominal",1.0)),
    axialExpDList(zoneNumber),
    axialExpNuSigmaEffList(zoneNumber),
    axialExpSigmaPowList(zoneNumber),
    axialExpSigmaDisappList(zoneNumber),
    axialExpSigmaFromToList(zoneNumber),
    RadExp(nuclearDataRadialExp.lookupOrDefault("expansionFromNominal",1.0)),
    axialOrientation(nuclearDataRadialExp.lookupOrDefault("axialOrientation",vector(0.0, 0.0, 1.0))),
    radialExpDList(zoneNumber),
    radialExpNuSigmaEffList(zoneNumber),
    radialExpSigmaPowList(zoneNumber),
    radialExpSigmaDisappList(zoneNumber),
    radialExpSigmaFromToList(zoneNumber),
    rhoCoolRef(nuclearDataRhoCool.lookupOrDefault("rhoCoolRef",860.0)),
    rhoCoolPerturbed(nuclearDataRhoCool.lookupOrDefault("rhoCoolPerturbed",1.0)),
    rhoCoolDList(zoneNumber),
    rhoCoolNuSigmaEffList(zoneNumber),
    rhoCoolSigmaPowList(zoneNumber),
    rhoCoolSigmaDisappList(zoneNumber),
    rhoCoolSigmaFromToList(zoneNumber),
    TCoolRef(nuclearDataTCool.lookupOrDefault("TCoolRef",900.0)),
    TCoolPerturbed(nuclearDataTCool.lookupOrDefault("TCoolPerturbed",1200.0)),
    TCoolDList(zoneNumber),
    TCoolNuSigmaEffList(zoneNumber),
    TCoolSigmaPowList(zoneNumber),
    TCoolSigmaDisappList(zoneNumber),
    TCoolSigmaFromToList(zoneNumber),
    TcladRef(nuclearDataCladExp.lookupOrDefault("TcladRef",900.0)),
    TcladPerturbed(nuclearDataCladExp.lookupOrDefault("TcladPerturbed",1200.0)),
    cladExpDList(zoneNumber),
    cladExpNuSigmaEffList(zoneNumber),
    cladExpSigmaPowList(zoneNumber),
    cladExpSigmaDisappList(zoneNumber),
    cladExpSigmaFromToList(zoneNumber),
    CRmove
    (
        IOobject
        (
            "CRmove",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    CRentries(CRmove.lookup("zones")),
    CRNumber(CRmove.lookup("zones").size()),
    CRstart(CRNumber),
    CRfinish(CRNumber),
    CRspeed(CRNumber),
    CRFollowerName(CRNumber),
    CRinitialPosition(CRNumber),
    CRposition(CRNumber),
    initialDistanceFromMeshCR(CRNumber)        
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
    const volScalarField& Tfuel, 
    const volScalarField& Tclad, 
    const volScalarField& rhoCool, 
    const volScalarField& TCool,
    const volVectorField& Disp
)
{
    #include "setNeutronicsVariables.H"
}

void Foam::XS::init()
{
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




// ************************************************************************* //
