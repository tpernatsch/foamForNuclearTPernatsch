/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2406                                                  |
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

#if defined __has_include
#  if __has_include(<commDataLayer.H>)
#    include <commDataLayer.H>
#    define isCommDataLayerIncluded
#  endif
#endif

#ifdef isCommDataLayerIncluded

#include "nuclearFuelFMU.H"
#include "structure.H"
#include "addToRunTimeSelectionTable.H"
#include "SquareMatrix.H"
#include "LUscalarMatrix.H"
#include "commDataLayer.H"
#include "latticeMap.H"
#include "listConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerModels
{
    defineTypeNameAndDebug(nuclearFuelFMU, 0);
    addToRunTimeSelectionTable
    (
        powerModel,
        nuclearFuelFMU,
        powerModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModels::nuclearFuelFMU::nuclearFuelFMU
(
    structure& structureRef,
    const dictionary& dicts
)
:
    powerModel
    (
        structureRef,
        dicts
    ),
    fmiState_
    (
        IOobject
        (
            "nuclearFuelFMU",
            mesh_.time().timeName(),
            "uniform/fmiState",
            mesh_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
    Tfav_
    (
        IOobject
        (
            "Tfav."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tcav_
    (
        IOobject
        (
            "Tcav."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tsurface_
    (
        IOobject
        (
            "Tsurface."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    interpolationCorrectFactor_
    (
        IOobject
        (
            "interpolationCorrectFactor."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    isHeatFluxInput_(0),
    isRhoCpdTdtInput_(0),
    fractionOfPowerFromNeutronics_(0),
    fuelFraction_(0),
    isFuelTemperatureFieldFromFMU_(0),
    isCladTemperatureFieldFromFMU_(0),
    relaxationFactorTfluid_(0),
    relaxationFactorTstruct_(0),
    relaxationFactorProfilePowerDensity_(0),
    relaxationFactorHtc_(0),
    oldTfluid_(0),
    oldTstruct_(0),
    oldProfilePowerDensity_(0),
    oldHtc_(0),
    xPos_(0),
    yPos_(0),
    axialLoc_(0),
    fuelLength_(0),
    zMethod_(0),
    radialBasisFunctionMethod_(0),
    epsilon_(0),
    krigingType_(0),
    krigingA_(0),
    krigingC_(0),
    krigingC0_(0),
    krigingWeight_(0),
    krigingTypeFuel_(0),
    krigingAfuel_(0),
    krigingCfuel_(0),
    krigingC0fuel_(0),
    krigingWeightFuel_(0),
    krigingTypeClad_(0),
    krigingAclad_(0),
    krigingCclad_(0),
    krigingC0clad_(0),
    krigingWeightClad_(0),
    regionIndexToRegionName_(0),
    cellToRegion_(mesh_.cells().size(), 0)
{
    structure_.setRegionField(*this, structureRef.powerDensityNeutronics(), "initPowerDensity");

    this->setInterfacialArea();

    // Communicating with the FMU
    const Time& runTime = this->db().time();
    commDataLayer& data = commDataLayer::New(runTime);

    // Resize usefull objects for faster field reconstruction
    krigingWeight_.resize(this->toc().size());
    krigingWeightFuel_.resize(this->toc().size());
    krigingWeightClad_.resize(this->toc().size());
    invRBFmatrix_.resize(this->toc().size());
    // invRBFintegralMatrix_.resize(this->toc().size());
    invRBFmatrixFuel_.resize(this->toc().size());
    invRBFmatrixClad_.resize(this->toc().size());

    forAll(this->toc(), regioni)
    {
        const word region(this->toc()[regioni]);
        const dictionary& dict(this->subDict(region));

        // Setup cellToRegion_ mapping
        const labelList& regionCells(structure_.cellLists()[region]);
        forAll(regionCells, i)
        {
            const label celli(regionCells[i]);
            cellToRegion_[celli] = regioni;
        }

        // Read region dict entries
        scalar fractionOfPowerFromNeutronics
        (
            dict.lookupOrDefault<scalar>("fractionOfPowerFromNeutronics", 1.0)
        );
        scalar fuelFraction
        (
            dict.lookupOrDefault<scalar>("fuelFraction", 1.0)
        );
        scalarList xPos(0);
        scalarList yPos(0);
        const scalarField axialLoc(dict.get<scalarField>("axialLocations"));
        const scalar fuelLength(dict.lookupOrDefault<scalar>
        (
            "fuelLength", axialLoc[axialLoc.size()-1] - axialLoc[0]
        ));
        const word radialBasisFunctionMethod(dict.get<word>("radialBasisFunctionMethod"));

        // Get XY sampling positions
        {
            if (dict.found("xPos") && dict.found("yPos"))
            {
                Info<< "nuclearFuelFMU uses XY posisions from xPos and yPos in region "
                    << region << " (list mode)"
                    << endl;

                xPos = dict.get<scalarList>("xPos");
                yPos = dict.get<scalarList>("yPos");
            }
            else if (dict.found("xyPosLattice"))
            {
                Info<< "nuclearFuelFMU uses XY posisions from xyPosLattice in region "
                    << region << " (lattice mode)"
                    << endl;

                const latticeMap& latticeMap(dict.subDict("xyPosLattice"));
                xPos = latticeMap.getXpositions();
                yPos = latticeMap.getYpositions();
            }
            else
            {
                FatalErrorInFunction
                    << "None of (xPos, yPos) or xyPosLattice have been provided in region " << region
                    << ". It is mandatory to add XY-locations samples."
                    << exit(FatalError);
            }
        }

        // Fill in lists for this region
        fractionOfPowerFromNeutronics_.append(fractionOfPowerFromNeutronics);
        fuelFraction_.append(fuelFraction);
        const bool isFuelTemperatureFieldFromFMU(dict.found("TFuelNameFromFMU"));
        const bool isCladTemperatureFieldFromFMU(dict.found("TCladNameFromFMU"));
        const bool isHeatFluxInput(dict.found("heatFluxNameFromFMU"));
        const bool isRhoCpdTdtInput(dict.found("rhoCpdTdtNameFromFMU"));
        isFuelTemperatureFieldFromFMU_.append(isFuelTemperatureFieldFromFMU);
        isCladTemperatureFieldFromFMU_.append(isCladTemperatureFieldFromFMU);
        isHeatFluxInput_.append(isHeatFluxInput);
        isRhoCpdTdtInput_.append(isRhoCpdTdtInput);
        xPos_.append(xPos);
        yPos_.append(yPos);
        axialLoc_.append(axialLoc);
        fuelLength_.append(fuelLength);
        zMethod_.append((InterpolateTableBaseGF::interpolationMethodNames_[
            dict.lookupOrDefault<word>("axialPowerInterpolationMethod", "linear")
        ]));
        radialBasisFunctionMethod_.append(radialBasisFunctionMethod);

        // Compute recommanded RBF epsilon factors for Gaussian method
        if (radialBasisFunctionMethod == "gaussian")
        {
            epsilon_.append(dict.get<scalarList>("epsilon"));

            const scalar epsilonX(computeEspilonForRBF(xPos));
            const scalar epsilonY(computeEspilonForRBF(yPos));
            const scalar epsilonZ(computeEspilonForRBF(axialLoc));
            Info<< "Recommended RBF epsilon parameters in cellZone " << region
                << nl
                << "   Epsilon X = " << epsilonX << nl
                << "   Epsilon Y = " << epsilonY << nl
                << "   Epsilon Z = " << epsilonZ
                << endl;
        }
        else if (radialBasisFunctionMethod == "kriging")
        {
            const dictionary& krigingOptions(dict.subDict("krigingOptions"));

            const word krigingType(krigingOptions.get<word>("type"));
            if (krigingType == "spherical") {
                krigingType_.append(Foam::radialBasisFunctionInterpolation::spherical);
            } else if (krigingType == "gaussian") {
                krigingType_.append(Foam::radialBasisFunctionInterpolation::gaussian);
            } else if (krigingType == "exponential") {
                krigingType_.append(Foam::radialBasisFunctionInterpolation::exponential);
            }
            krigingA_.append(krigingOptions.get<scalar>("a"));
            krigingC_.append(krigingOptions.get<scalar>("c"));
            krigingC0_.append(krigingOptions.get<scalar>("c0"));

            if (isFuelTemperatureFieldFromFMU)
            {
                const dictionary& krigingOptionsFuel(dict.subDict("krigingOptionsFuel"));

                const word krigingTypeFuel(krigingOptionsFuel.get<word>("type"));
                if (krigingTypeFuel == "spherical") {
                    krigingTypeFuel_.append(Foam::radialBasisFunctionInterpolation::spherical);
                } else if (krigingTypeFuel == "gaussian") {
                    krigingTypeFuel_.append(Foam::radialBasisFunctionInterpolation::gaussian);
                } else if (krigingTypeFuel == "exponential") {
                    krigingTypeFuel_.append(Foam::radialBasisFunctionInterpolation::exponential);
                }
                krigingAfuel_.append(krigingOptionsFuel.get<scalar>("a"));
                krigingCfuel_.append(krigingOptionsFuel.get<scalar>("c"));
                krigingC0fuel_.append(krigingOptionsFuel.get<scalar>("c0"));
            }

            if (isCladTemperatureFieldFromFMU)
            {
                const dictionary& krigingOptionsClad(dict.subDict("krigingOptionsClad"));

                const word krigingTypeClad(krigingOptionsClad.get<word>("type"));
                if (krigingTypeClad == "spherical") {
                    krigingTypeClad_.append(Foam::radialBasisFunctionInterpolation::spherical);
                } else if (krigingTypeClad == "gaussian") {
                    krigingTypeClad_.append(Foam::radialBasisFunctionInterpolation::gaussian);
                } else if (krigingTypeClad == "exponential") {
                    krigingTypeClad_.append(Foam::radialBasisFunctionInterpolation::exponential);
                }
                krigingAclad_.append(krigingOptionsClad.get<scalar>("a"));
                krigingCclad_.append(krigingOptionsClad.get<scalar>("c"));
                krigingC0clad_.append(krigingOptionsClad.get<scalar>("c0"));
            }
        }

        //- Add to regionIndexToRegionName_ mapping
        regionIndexToRegionName_.append(region);

        // Set default values from phaseProperties
        const word axialLocNameToFMU(dict.get<word>("axialLocationsNameToFMU"));
        const wordList avgPowerDensityNameToFMU(dict.get<wordList>("avgPowerDensityNameToFMU"));
        const wordList axialProfilePowerDensityNameToFMU(dict.get<wordList>("axialProfilePowerDensityNameToFMU"));
        const wordList TstructNameToFMU(dict.lookupOrDefault<wordList>("TstructNameToFMU", {}));
        const wordList TfluidNameToFMU(dict.lookupOrDefault<wordList>("TfluidNameToFMU", {}));
        const wordList htcNameToFMU(dict.lookupOrDefault<wordList>("htcNameToFMU", {}));
        const wordList TFuelNameFromFMU
        (
            isFuelTemperatureFieldFromFMU
                ? dict.get<wordList>("TFuelNameFromFMU")
                : wordList(0)
        );
        const wordList TCladNameFromFMU
        (
            isCladTemperatureFieldFromFMU
                ? dict.get<wordList>("TCladNameFromFMU")
                : wordList(0)
        );
        wordList heatFluxNameFromFMU(0);
        wordList rhoCpdTdtNameFromFMU(0);

        const scalar initPowerDensityFromDict(dict.get<scalar>("initPowerDensity"));
        const scalar initTemperatureFromDict(dict.get<scalar>("initTstruct"));
        relaxationFactorTfluid_.append(dict.lookupOrDefault<scalar>("relaxationFactorTfluid", 1.0));
        relaxationFactorTstruct_.append(dict.lookupOrDefault<scalar>("relaxationFactorTstruct", 1.0));
        relaxationFactorProfilePowerDensity_.append(dict.lookupOrDefault<scalar>("relaxationFactorPowerDensity", 1.0));
        relaxationFactorHtc_.append(dict.lookupOrDefault<scalar>("relaxationFactorHtc", 1.0));

        // Init old value lists for relaxation
        List<List<scalar>> oldTfluidTemp_(0);
        List<List<scalar>> oldTstructTemp_(0);
        List<List<scalar>> oldProfilePowerDensityTemp_(0);
        List<List<scalar>> oldHtcTemp_(0);
        forAll(avgPowerDensityNameToFMU, nameI)
        {
            oldTfluidTemp_.append(scalarList(axialLoc.size(), initTemperatureFromDict));
            oldTstructTemp_.append(scalarList(axialLoc.size(), initTemperatureFromDict));
            oldProfilePowerDensityTemp_.append(scalarList(axialLoc.size(), 1.0));
            oldHtcTemp_.append(scalarList(axialLoc.size(), 0.0));
        }
        oldTfluid_.append(oldTfluidTemp_);
        oldTstruct_.append(oldTstructTemp_);
        oldProfilePowerDensity_.append(oldProfilePowerDensityTemp_);
        oldHtc_.append(oldHtcTemp_);

        scalar initHeatFluxFromDict(0);
        scalar initRhoCpdTdtFromDict(0);

        if (isHeatFluxInput)
        {
            heatFluxNameFromFMU = dict.get<wordList>("heatFluxNameFromFMU");
            initHeatFluxFromDict = dict.get<scalar>("initHeatFlux");
        }
        if (isRhoCpdTdtInput)
        {
            rhoCpdTdtNameFromFMU = dict.get<wordList>("rhoCpdTdtNameFromFMU");
            initRhoCpdTdtFromDict = dict.get<scalar>("initRhoCpdTdt");
        }

        // Initialize structure temperature
        forAll(regionCells, i)
        {
            const label celli(regionCells[i]);
            Tfav_[celli] = initTemperatureFromDict;
            Tcav_[celli] = initTemperatureFromDict;
            Tsurface_[celli] = initTemperatureFromDict;
        }

        // Check lengths of FMI port lists
        {
            const label nx(xPos.size());

            if (nx != yPos.size())
            {
                FatalErrorInFunction
                    << "xPos and yPos lengths are different in region " << region
                    << exit(FatalError);
            }
            if (nx != avgPowerDensityNameToFMU.size())
            {
                FatalErrorInFunction
                    << "xPos and avgPowerDensityNameToFMU list lengths are "
                    << "different in region " << region
                    << exit(FatalError);
            }
            if (nx != axialProfilePowerDensityNameToFMU.size())
            {
                FatalErrorInFunction
                    << "xPos and axialProfilePowerDensityNameToFMU list lengths are "
                    << "different in region " << region
                    << exit(FatalError);
            }
            if (!TstructNameToFMU.empty() && nx != TstructNameToFMU.size())
            {
                FatalErrorInFunction
                    << "xPos and TstructNameToFMU list lengths are different in "
                    << "region " << region
                    << exit(FatalError);
            }
            if (!TfluidNameToFMU.empty() && nx != TfluidNameToFMU.size())
            {
                FatalErrorInFunction
                    << "xPos and TfluidNameToFMU list lengths are different in "
                    << "region " << region
                    << exit(FatalError);
            }
            if (!htcNameToFMU.empty() && nx != htcNameToFMU.size())
            {
                FatalErrorInFunction
                    << "xPos and htcNameToFMU list lengths are different in "
                    << "region " << region
                    << exit(FatalError);
            }
            if (isHeatFluxInput && nx != heatFluxNameFromFMU.size())
            {
                FatalErrorInFunction
                    << "xPos and heatFluxNameFromFMU list lengths are different "
                    << "in region " << region
                    << exit(FatalError);
            }
            if (isRhoCpdTdtInput && nx != rhoCpdTdtNameFromFMU.size())
            {
                FatalErrorInFunction
                    << "xPos and rhoCpdTdtNameFromFMU list lengths are different "
                    << "in region " << region
                    << exit(FatalError);
            }
            if (isFuelTemperatureFieldFromFMU && nx != TFuelNameFromFMU.size())
            {
                FatalErrorInFunction
                    << "xPos and TFuelNameFromFMU list lengths are different "
                    << "in region " << region
                    << exit(FatalError);
            }
            if (isCladTemperatureFieldFromFMU && nx != TCladNameFromFMU.size())
            {
                FatalErrorInFunction
                    << "xPos and TCladNameFromFMU list lengths are different "
                    << "in region " << region
                    << exit(FatalError);
            }
        }

        // Axial location
        word axialLocStringified(Foam::listConversion::stringify<scalar>(axialLoc));
        data.storeObj
        (
            axialLocStringified,
            axialLocNameToFMU,
            commDataLayer::causality::out
        );

        // The initialization needs to be improved, not clean enough, need the
        // full distribution for each pin
        Info<< "nuclearFuelFMU region " << region << endl;
        forAll(avgPowerDensityNameToFMU, nameI)
        {
            Info<< "    " << xPos[nameI] << " " << yPos[nameI] << " : "
                << avgPowerDensityNameToFMU[nameI] << " "
                << axialProfilePowerDensityNameToFMU[nameI] << " ";
                // << TstructNameToFMU[nameI] << " "
                if (!heatFluxNameFromFMU.empty())
                {
                    Info<< heatFluxNameFromFMU[nameI] << " ";
                }
                if (!rhoCpdTdtNameFromFMU.empty())
                {
                    Info<< rhoCpdTdtNameFromFMU[nameI] << " ";
                }
                if (!TstructNameToFMU.empty())
                {
                    Info<< TstructNameToFMU[nameI] << " ";
                }
                if (!TfluidNameToFMU.empty())
                {
                    Info<< TfluidNameToFMU[nameI] << " ";
                }
                if (!htcNameToFMU.empty())
                {
                    Info<< htcNameToFMU[nameI] << " ";
                }
                if (isFuelTemperatureFieldFromFMU)
                {
                    Info<< TFuelNameFromFMU[nameI] << " ";
                }
                if (isCladTemperatureFieldFromFMU)
                {
                    Info<< TCladNameFromFMU[nameI];
                }
                Info<< endl;

            // Average power density
            data.storeObj
            (
                fmiState_.lookupOrDefault<scalar>
                (
                    avgPowerDensityNameToFMU[nameI], initPowerDensityFromDict
                ),
                avgPowerDensityNameToFMU[nameI],
                commDataLayer::causality::out
            );

            // Axial distribution of power density
            scalarList initPowerDensity
            (
                fmiState_.lookupOrDefault<scalarList>
                (
                    axialProfilePowerDensityNameToFMU[nameI],
                    oldProfilePowerDensityTemp_[nameI]
                    // scalarList(axialLoc.size(), 1.0)
                )
            );
            data.storeObj
            (
                Foam::listConversion::stringify(initPowerDensity),
                axialProfilePowerDensityNameToFMU[nameI],
                commDataLayer::causality::out
            );

            // Temperature
            if (!TstructNameToFMU.empty())
            {
                scalarList initTstruct
                (
                    fmiState_.lookupOrDefault<scalarList>
                    (
                        TstructNameToFMU[nameI],
                        oldTstructTemp_[nameI]
                        // scalarList(axialLoc.size(), initTemperatureFromDict)
                    )
                );
                data.storeObj
                (
                    // initTstruct,
                    Foam::listConversion::stringify(initTstruct),
                    TstructNameToFMU[nameI],
                    commDataLayer::causality::out
                );
            }
            if (!TfluidNameToFMU.empty())
            {
                scalarList initTfluid
                (
                    fmiState_.lookupOrDefault<scalarList>
                    (
                        TfluidNameToFMU[nameI],
                        oldTfluidTemp_[nameI]
                        // scalarList(axialLoc.size(), initTemperatureFromDict)
                    )
                );
                data.storeObj
                (
                    Foam::listConversion::stringify(initTfluid),
                    TfluidNameToFMU[nameI],
                    commDataLayer::causality::out
                );
            }
            if (!htcNameToFMU.empty())
            {
                scalarList initHtc
                (
                    fmiState_.lookupOrDefault<scalarList>
                    (
                        htcNameToFMU[nameI],
                        oldHtcTemp_[nameI]
                        // scalarList(axialLoc.size(), 0.0)
                    )
                );
                data.storeObj
                (
                    Foam::listConversion::stringify(initHtc),
                    htcNameToFMU[nameI],
                    commDataLayer::causality::out
                );
            }

            if (isFuelTemperatureFieldFromFMU)
            {
                // Fuel temperature
                scalarList initTFuel
                (
                    fmiState_.lookupOrDefault<scalarList>
                    (
                        TFuelNameFromFMU[nameI],
                        scalarList(axialLoc.size(), initTemperatureFromDict)
                    )
                );
                data.storeObj
                (
                    Foam::listConversion::stringify(initTFuel),
                    TFuelNameFromFMU[nameI],
                    commDataLayer::causality::in
                );
            }
            if (isCladTemperatureFieldFromFMU)
            {
                // Cladding temperature
                scalarList initTClad
                (
                    fmiState_.lookupOrDefault<scalarList>
                    (
                        TCladNameFromFMU[nameI],
                        scalarList(axialLoc.size(), initTemperatureFromDict)
                    )
                );
                data.storeObj
                (
                    Foam::listConversion::stringify(initTClad),
                    TCladNameFromFMU[nameI],
                    commDataLayer::causality::in
                );
            }

            // Heat flux
            if (isHeatFluxInput)
            {
                scalarList initHeatFlux
                (
                    fmiState_.lookupOrDefault<scalarList>
                    (
                        heatFluxNameFromFMU[nameI],
                        scalarList(axialLoc.size(), initHeatFluxFromDict)
                    )
                );
                data.storeObj
                (
                    Foam::listConversion::stringify(initHeatFlux),
                    heatFluxNameFromFMU[nameI],
                    commDataLayer::causality::in
                );
            }
            if (isRhoCpdTdtInput)
            {
                scalarList initRhoCpdTdt
                (
                    fmiState_.lookupOrDefault<scalarList>
                    (
                        rhoCpdTdtNameFromFMU[nameI],
                        scalarList(axialLoc.size(), initRhoCpdTdtFromDict)
                    )
                );
                data.storeObj
                (
                    Foam::listConversion::stringify(initRhoCpdTdt),
                    rhoCpdTdtNameFromFMU[nameI],
                    commDataLayer::causality::in
                );
            }
        }
        Info<< endl;
    }

    Tfav_.correctBoundaryConditions();
    Tcav_.correctBoundaryConditions();
    Tsurface_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerModels::nuclearFuelFMU::~nuclearFuelFMU()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::powerModels::nuclearFuelFMU::correct
(
    const volScalarField& HTSum,  // == SUM_j [htc_j*T_j*frac_j]
    const volScalarField& HSum    // == SUM_j [htc_j*frac_j]
)
{
    forAll(this->toc(), regioni)
    {
        // Update surface temperature from FMU using heat flux from FMUs
        correctHeatFluxInputsFromFMUs(HTSum, HSum, regioni);

        // Update surface temperature and power density sample lines for FMUs
        correctInputsForFMUs(HTSum, HSum, regioni);

        if (mesh_.time().writeTime())
        {
            fmiState_.set(this->toc()[regioni], axialLoc_[regioni]);
        }
    }
}

void Foam::powerModels::nuclearFuelFMU::correctHeatFluxInputsFromFMUs
(
    const volScalarField& HTSum,  // == SUM_j [htc_j*T_j*frac_j]
    const volScalarField& HSum,   // == SUM_j [htc_j*frac_j]
    const label regioni
)
{
    // Communicating with the FMU, extract default values stored in
    // commDataLayer
    const Time& runTime = this->db().time();
    commDataLayer& data = commDataLayer::New(runTime);

    // Volume of the cells
    const scalarList& V(mesh_.V());

    // Get dict in the cellZone power model
    word region(this->toc()[regioni]);
    const labelList& regionCells(structure_.cellLists()[region]);
    const dictionary& dict(this->subDict(region));

    const bool isHeatFluxInput(isHeatFluxInput_[regioni]);
    const bool isRhoCpdTdtInput(isRhoCpdTdtInput_[regioni]);
    const bool isFuelTemperatureFieldFromFMU(isFuelTemperatureFieldFromFMU_[regioni]);
    const bool isCladTemperatureFieldFromFMU(isCladTemperatureFieldFromFMU_[regioni]);
    const word radialBasisFunctionMethod(radialBasisFunctionMethod_[regioni]);

    const scalarField axialLoc(axialLoc_[regioni]);

    // Get field name for heat flux reconstruction
    wordList heatFluxFieldNameFromFMU(0);
    wordList rhoCpdTdtFieldNameFromFMU(0);
    if (isHeatFluxInput)
    {
        heatFluxFieldNameFromFMU = dict.get<wordList>("heatFluxNameFromFMU");
    }
    if (isRhoCpdTdtInput)
    {
        rhoCpdTdtFieldNameFromFMU = dict.get<wordList>("rhoCpdTdtNameFromFMU");
    }

    // Get field name for the fuel pin temperature (how to generalize it for
    // pebble bed?)
    wordList TFuelFieldNameFromFMU(0);
    wordList TCladFieldNameFromFMU(0);
    if (isFuelTemperatureFieldFromFMU)
    {
        TFuelFieldNameFromFMU = dict.get<wordList>("TFuelNameFromFMU");
    }
    if (isCladTemperatureFieldFromFMU)
    {
        TCladFieldNameFromFMU = dict.get<wordList>("TCladNameFromFMU");
    }

    // List of positions and values for the interpolation algorithm
    const label nFMIports(max(heatFluxFieldNameFromFMU.size(), rhoCpdTdtFieldNameFromFMU.size()));
    const label nAxialLoc(axialLoc.size());
    scalarList xPosList(nFMIports * nAxialLoc);
    scalarList yPosList(nFMIports * nAxialLoc);
    scalarList zPosList(nFMIports * nAxialLoc);
    scalarList heatFluxFieldFromFMUList(nFMIports * nAxialLoc);
    scalarList rhoCpdTdtFieldFromFMUList(nFMIports * nAxialLoc);
    // scalarList totalPowerNormalizedFieldFromFMUList
    // (
    //     isHeatFluxInput && isRhoCpdTdtInput
    //         ? nFMIports * nAxialLoc
    //         : 0,
    //     0.0
    // );
    scalarList TFuelFieldFromFMUList(nFMIports * nAxialLoc);
    scalarList TCladFieldFromFMUList(nFMIports * nAxialLoc);

    // Extract x, y, z, heat flux (or enthalpy) into list
    // forAll(heatFluxFieldNameFromFMU, nameI)
    for (label nameI=0; nameI < nFMIports; nameI++)
    {
        const scalar& xPos(xPos_[regioni][nameI]);
        const scalar& yPos(yPos_[regioni][nameI]);

        // x, y, z positions
        forAll(axialLoc, locI)
        {
            xPosList[nAxialLoc * nameI + locI] = xPos;
            yPosList[nAxialLoc * nameI + locI] = yPos;
            zPosList[nAxialLoc * nameI + locI] = axialLoc[locI];
        }

        // Extract heat flux or enthalpy
        word& heatFluxFieldFromFMU = data.getObj<word>
        (
            !heatFluxFieldNameFromFMU.empty()
                ? heatFluxFieldNameFromFMU[nameI]
                : rhoCpdTdtFieldNameFromFMU[nameI],
            commDataLayer::causality::in
        );
        word& rhoCpdTdtFieldFromFMU = data.getObj<word>
        (
            !rhoCpdTdtFieldNameFromFMU.empty()
                ? rhoCpdTdtFieldNameFromFMU[nameI]
                : heatFluxFieldNameFromFMU[nameI],
            commDataLayer::causality::in
        );

        // Dirty way but needs to be instantiated
        // To avoid errors, give the FMI from heat flux, but is is 'fine'
        // because the value is never used
        word& TFuelFieldFromFMU = data.getObj<word>
        (
            isFuelTemperatureFieldFromFMU
                ? TFuelFieldNameFromFMU[nameI]
                : (
                    !heatFluxFieldNameFromFMU.empty()
                        ? heatFluxFieldNameFromFMU[nameI]
                        : rhoCpdTdtFieldNameFromFMU[nameI]
                ),
            commDataLayer::causality::in
        );
        word& TCladFieldFromFMU = data.getObj<word>
        (
            isCladTemperatureFieldFromFMU
                ? TCladFieldNameFromFMU[nameI]
                : (
                    !heatFluxFieldNameFromFMU.empty()
                        ? heatFluxFieldNameFromFMU[nameI]
                        : rhoCpdTdtFieldNameFromFMU[nameI]
                ),
            commDataLayer::causality::in
        );

        // Split string to list of scalar
        word value;
        scalarList heatFluxFieldFromFMUListTemp(0);
        scalarList rhoCpdTdtFieldFromFMUListTemp(0);
        scalarList TFuelFieldFromFMUListTemp(0);
        scalarList TCladFieldFromFMUListTemp(0);

        if (isHeatFluxInput)
        {
            std::stringstream ssHeatFlux(heatFluxFieldFromFMU);
            while (getline(ssHeatFlux, value, ' '))
            {
                heatFluxFieldFromFMUListTemp.append(std::stod(value));
            }
        }
        if (isRhoCpdTdtInput)
        {
            std::stringstream ssRhoCpdTdt(rhoCpdTdtFieldFromFMU);
            while (getline(ssRhoCpdTdt, value, ' '))
            {
                rhoCpdTdtFieldFromFMUListTemp.append(std::stod(value));
            }
        }
        if (isFuelTemperatureFieldFromFMU)
        {
            std::stringstream ssTFuel(TFuelFieldFromFMU);
            while (getline(ssTFuel, value, ' '))
            {
                TFuelFieldFromFMUListTemp.append(std::stod(value));
            }
        }
        if (isCladTemperatureFieldFromFMU)
        {
            std::stringstream ssTClad(TCladFieldFromFMU);
            while (getline(ssTClad, value, ' '))
            {
                TCladFieldFromFMUListTemp.append(std::stod(value));
            }
        }

        // Fill heat flux (or enthalpy) list
        if (heatFluxFieldFromFMUListTemp.size() == 1)
        {
            if (isHeatFluxInput)
            {
                forAll(axialLoc, locI)
                {
                    heatFluxFieldFromFMUList[nAxialLoc * nameI + locI] = heatFluxFieldFromFMUListTemp[0];
                }
            }
            if (isRhoCpdTdtInput)
            {
                forAll(axialLoc, locI)
                {
                    rhoCpdTdtFieldFromFMUList[nAxialLoc * nameI + locI] = rhoCpdTdtFieldFromFMUListTemp[0];
                }
            }
            if (isFuelTemperatureFieldFromFMU)
            {
                forAll(axialLoc, locI)
                {
                    TFuelFieldFromFMUList[nAxialLoc * nameI + locI] = TFuelFieldFromFMUListTemp[0];
                }
            }
            if (isCladTemperatureFieldFromFMU)
            {
                forAll(axialLoc, locI)
                {
                    TCladFieldFromFMUList[nAxialLoc * nameI + locI] = TCladFieldFromFMUListTemp[0];
                }
            }
        }
        else if (heatFluxFieldFromFMUListTemp.size() == axialLoc.size())
        {
            if (isHeatFluxInput)
            {
                forAll(axialLoc, locI)
                {
                    heatFluxFieldFromFMUList[nAxialLoc * nameI + locI] = heatFluxFieldFromFMUListTemp[locI];
                }
            }
            if (isRhoCpdTdtInput)
            {
                forAll(axialLoc, locI)
                {
                    rhoCpdTdtFieldFromFMUList[nAxialLoc * nameI + locI] = rhoCpdTdtFieldFromFMUListTemp[locI];
                }
            }
            if (isFuelTemperatureFieldFromFMU)
            {
                forAll(axialLoc, locI)
                {
                    TFuelFieldFromFMUList[nAxialLoc * nameI + locI] = TFuelFieldFromFMUListTemp[locI];
                }
            }
            if (isCladTemperatureFieldFromFMU)
            {
                forAll(axialLoc, locI)
                {
                    TCladFieldFromFMUList[nAxialLoc * nameI + locI] = TCladFieldFromFMUListTemp[locI];
                }
            }
        }
        else
        {
            FatalErrorInFunction
                << "No heat flux list provided in region " << region
                << " for FMI port " << heatFluxFieldNameFromFMU[nameI]
                << exit(FatalError);
        }

        if (mesh_.time().writeTime()) // This .write() is slowing down the code (x20)!
        {
            if (isHeatFluxInput)
            {
                fmiState_.set
                (
                    heatFluxFieldNameFromFMU[nameI], heatFluxFieldFromFMUListTemp
                );
            }
            if (isRhoCpdTdtInput)
            {
                fmiState_.set
                (
                    rhoCpdTdtFieldNameFromFMU[nameI], rhoCpdTdtFieldFromFMUListTemp
                );
            }
            if (isFuelTemperatureFieldFromFMU)
            {
                fmiState_.set
                (
                    TFuelFieldNameFromFMU[nameI], TFuelFieldFromFMUListTemp
                );
            }
            if (isCladTemperatureFieldFromFMU)
            {
                fmiState_.set
                (
                    TCladFieldNameFromFMU[nameI], TCladFieldFromFMUListTemp
                );
            }
        }
    }


    // Compute total neutronics power in cellZone
    scalar totalPowerNeutronics(0);
    scalar Vtot(0);
    Info<< "fractionOfPowerFromNeutronics = " << fractionOfPowerFromNeutronics_[regioni]
        << "; "
        << "fuelFraction = " << fuelFraction_[regioni]
        << endl;
    forAll(regionCells, i)
    {
        const label celli(regionCells[i]);

        const scalar powerDensity
        (
            fractionOfPowerFromNeutronics_[regioni]
                * fuelFraction_[regioni]
                * structure_.powerDensityNeutronics()[celli]
        );

        totalPowerNeutronics += powerDensity * V[celli];

        Vtot += V[celli];
    }

    Info<< "Computed totalPowerNeutronics " << totalPowerNeutronics << " W" << endl;


    // Generate weight for RBF (only polyharmonicSpline and gaussian)
    scalarList interpolationWeightsHeatFlux(0);
    scalarList interpolationWeightsRhoCpdTdt(0);
    // scalarList interpolationWeightsCorrectionFactor(0);
    scalarList interpolationWeightsTFuel(0);
    scalarList interpolationWeightsTClad(0);
    if (radialBasisFunctionMethod == "polyharmonicSpline")
    {
        // if (isHeatFluxInput && isRhoCpdTdtInput)
        // {
        //     if (totalPowerNeutronics > 0)
        //     {
        //         const scalar iA(dict.get<scalar>("volumetricArea"));
        //         // const scalar alpha(dict.get<scalar>("volumeFraction"));

        //         forAll(totalPowerNormalizedFieldFromFMUList, i)
        //         {
        //             if (!isFuelTemperatureFieldFromFMU || TFuelFieldFromFMUList[i] > 0)
        //             {
        //                 totalPowerNormalizedFieldFromFMUList[i] =
        //                     mag(iA * heatFluxFieldFromFMUList[i] + /*alpha **/ rhoCpdTdtFieldFromFMUList[i]);
        //                     // /
        //                     // totalPowerNeutronics;
        //             }

        //             // Clipping to avoid huge disturbances in the correction of
        //             // the heat flux and rhocpdTdt RBFs
        //             // if (totalPowerNormalizedFieldFromFMUList[i] < 0.9 || 1.1 < totalPowerNormalizedFieldFromFMUList[i])
        //             // {
        //             //     totalPowerNormalizedFieldFromFMUList[i] = 1.0;
        //             // }
        //         }

        //         interpolationWeightsCorrectionFactor = Foam::radialBasisFunctionInterpolation::solvePolyharmonicSplineIntegral
        //         (
        //             xPosList, yPosList, zPosList, totalPowerNormalizedFieldFromFMUList,
        //             // 1.0, // Total integral has been normalized so 1
        //             totalPowerNeutronics, // Total integral = total power
        //             regionCells, mesh_, invRBFintegralMatrix_[regioni]
        //         );
        //     }
        // }
        if (isHeatFluxInput)
        {
            interpolationWeightsHeatFlux = Foam::radialBasisFunctionInterpolation::solvePolyharmonicSpline
            (
                xPosList, yPosList, zPosList, heatFluxFieldFromFMUList, invRBFmatrix_[regioni]
            );
        }
        if (isRhoCpdTdtInput)
        {
            interpolationWeightsRhoCpdTdt = Foam::radialBasisFunctionInterpolation::solvePolyharmonicSpline
            (
                xPosList, yPosList, zPosList, rhoCpdTdtFieldFromFMUList, invRBFmatrix_[regioni]
            );
        }
        if (isFuelTemperatureFieldFromFMU)
        {
            interpolationWeightsTFuel = Foam::radialBasisFunctionInterpolation::solvePolyharmonicSpline
            (
                xPosList, yPosList, zPosList, TFuelFieldFromFMUList, invRBFmatrix_[regioni]
            );
        }
        if (isCladTemperatureFieldFromFMU)
        {
            interpolationWeightsTClad = Foam::radialBasisFunctionInterpolation::solvePolyharmonicSpline
            (
                xPosList, yPosList, zPosList, TCladFieldFromFMUList, invRBFmatrix_[regioni]
            );
        }
    }
    else if (radialBasisFunctionMethod == "gaussian")
    {
        if (isHeatFluxInput)
        {
            interpolationWeightsHeatFlux = Foam::radialBasisFunctionInterpolation::solveGaussianRadialBasisFunction
            (
                xPosList, yPosList, zPosList, heatFluxFieldFromFMUList, epsilon_[regioni], invRBFmatrix_[regioni]
            );
        }
        if (isRhoCpdTdtInput)
        {
            interpolationWeightsRhoCpdTdt = Foam::radialBasisFunctionInterpolation::solveGaussianRadialBasisFunction
            (
                xPosList, yPosList, zPosList, rhoCpdTdtFieldFromFMUList, epsilon_[regioni], invRBFmatrix_[regioni]
            );
        }
        if (isFuelTemperatureFieldFromFMU)
        {
            interpolationWeightsTFuel = Foam::radialBasisFunctionInterpolation::solveGaussianRadialBasisFunction
            (
                xPosList, yPosList, zPosList, TFuelFieldFromFMUList, epsilon_[regioni], invRBFmatrix_[regioni]
            );
        }
        if (isCladTemperatureFieldFromFMU)
        {
            interpolationWeightsTClad = Foam::radialBasisFunctionInterpolation::solveGaussianRadialBasisFunction
            (
                xPosList, yPosList, zPosList, TCladFieldFromFMUList, epsilon_[regioni], invRBFmatrix_[regioni]
            );
        }
    }
    else if (radialBasisFunctionMethod == "kriging")
    {
        // Compute Kriging weights once for each cell of interest
        if (krigingWeight_[regioni].size() == 0)
        {
            krigingWeight_[regioni].resize(regionCells.size());

            if (isFuelTemperatureFieldFromFMU)
            {
                krigingWeightFuel_[regioni].resize(regionCells.size());
            }
            if (isCladTemperatureFieldFromFMU)
            {
                krigingWeightClad_[regioni].resize(regionCells.size());
            }

            forAll(regionCells, i)
            {
                const label celli(regionCells[i]);

                // Get cell position
                const scalar xCell(mesh_.C().internalField()[celli].x());
                const scalar yCell(mesh_.C().internalField()[celli].y());
                const scalar zCell(mesh_.C().internalField()[celli].z());

                krigingWeight_[regioni][i] = Foam::radialBasisFunctionInterpolation::solveKriging
                (
                    xPosList, yPosList, zPosList,
                    krigingType_[regioni],
                    krigingA_[regioni], krigingC_[regioni], krigingC0_[regioni],
                    xCell, yCell, zCell,
                    invRBFmatrix_[regioni]
                );

                if (isFuelTemperatureFieldFromFMU)
                {
                    krigingWeightFuel_[regioni][i] = Foam::radialBasisFunctionInterpolation::solveKriging
                    (
                        xPosList, yPosList, zPosList,
                        krigingTypeFuel_[regioni],
                        krigingAfuel_[regioni], krigingCfuel_[regioni], krigingC0fuel_[regioni],
                        xCell, yCell, zCell,
                        invRBFmatrixFuel_[regioni]
                    );
                }
                if (isCladTemperatureFieldFromFMU)
                {
                    krigingWeightClad_[regioni][i] = Foam::radialBasisFunctionInterpolation::solveKriging
                    (
                        xPosList, yPosList, zPosList,
                        krigingTypeClad_[regioni],
                        krigingAclad_[regioni], krigingCclad_[regioni], krigingC0clad_[regioni],
                        xCell, yCell, zCell,
                        invRBFmatrixClad_[regioni]
                    );
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << radialBasisFunctionMethod << " is an incorrect "
            << "radial basis function method provided in region " << region
            << ". Available methods: polyharmonicSpline, gaussian, kriging"
            << exit(FatalError);
    }

    scalar totalPowerFromHeatFlux(0), totalPowerEnthalpy(0);

    // Update all Tsurface cells in the region using heat flux interpolated
    forAll(regionCells, i)
    {
        const label celli(regionCells[i]);

        // Get coolant temperature
        const scalar HSumi(HSum[celli]);
        if (HSumi == 0)
        {
            FatalErrorInFunction
                << "No heat transfer model provided in region " << region
                << " (see phaseProperties.physicsModels.heatTransferModels)"
                << exit(FatalError);
        }
        const scalar Tcool(HTSum[celli] / max(HSumi, SMALL));

        // Get cell position
        const scalar xCell(mesh_.C().internalField()[celli].x());
        const scalar yCell(mesh_.C().internalField()[celli].y());
        const scalar zCell(mesh_.C().internalField()[celli].z());

        // Interpolate heat flux at x, y, z position
        scalar interpolatedValueHeatFlux(0);
        scalar interpolatedValueRhoCpdTdt(0);
        // scalar interpolatedValueTotalPower(1);
        scalar interpolatedValueTFuel(0);
        scalar interpolatedValueTClad(0);
        if (radialBasisFunctionMethod == "polyharmonicSpline")
        {
            if (isHeatFluxInput)
            {
                interpolatedValueHeatFlux = Foam::radialBasisFunctionInterpolation::polyharmonicSpline
                (
                    interpolationWeightsHeatFlux, xPosList, yPosList, zPosList,
                    xCell, yCell, zCell
                );
            }
            if (isRhoCpdTdtInput)
            {
                interpolatedValueRhoCpdTdt = Foam::radialBasisFunctionInterpolation::polyharmonicSpline
                (
                    interpolationWeightsRhoCpdTdt, xPosList, yPosList, zPosList,
                    xCell, yCell, zCell
                );
            }
            // if (isHeatFluxInput && isRhoCpdTdtInput && totalPowerNeutronics > 0)
            // {
            //     interpolatedValueTotalPower = Foam::radialBasisFunctionInterpolation::polyharmonicSplineIntegral
            //     (
            //         interpolationWeightsCorrectionFactor, xPosList, yPosList, zPosList,
            //         xCell, yCell, zCell
            //     );

            //     // if (interpolatedValueTotalPower > 0)
            //     // {
            //         const scalar k
            //         (
            //             interpolatedValueTotalPower
            //             /
            //             (iA_[celli] * interpolatedValueHeatFlux + interpolatedValueRhoCpdTdt)
            //         );


            //         if (0.5 < k && k < 1.5)
            //         {
            //             Info<< k << endl;
            //             interpolatedValueHeatFlux *= k;
            //             interpolatedValueRhoCpdTdt *= k;
            //         }
            //         else
            //         {
            //             interpolatedValueHeatFlux = 0.0;
            //             interpolatedValueRhoCpdTdt = 0.0;
            //         }
            //     // }
            // }

            if (isFuelTemperatureFieldFromFMU)
            {
                interpolatedValueTFuel = Foam::radialBasisFunctionInterpolation::polyharmonicSpline
                (
                    interpolationWeightsTFuel, xPosList, yPosList, zPosList,
                    xCell, yCell, zCell
                );
            }
            if (isCladTemperatureFieldFromFMU)
            {
                interpolatedValueTClad = Foam::radialBasisFunctionInterpolation::polyharmonicSpline
                (
                    interpolationWeightsTClad, xPosList, yPosList, zPosList,
                    xCell, yCell, zCell
                );
            }
        }
        else if (radialBasisFunctionMethod == "gaussian")
        {
            if (isHeatFluxInput)
            {
                interpolatedValueHeatFlux = Foam::radialBasisFunctionInterpolation::gaussianRadialBasisFunction
                (
                    interpolationWeightsHeatFlux, xPosList, yPosList, zPosList,
                    xCell, yCell, zCell,
                    epsilon_[regioni]
                );
            }
            else if (isRhoCpdTdtInput)
            {
                interpolatedValueRhoCpdTdt = Foam::radialBasisFunctionInterpolation::gaussianRadialBasisFunction
                (
                    interpolationWeightsRhoCpdTdt, xPosList, yPosList, zPosList,
                    xCell, yCell, zCell,
                    epsilon_[regioni]
                );
            }

            if (isFuelTemperatureFieldFromFMU)
            {
                interpolatedValueTFuel = Foam::radialBasisFunctionInterpolation::gaussianRadialBasisFunction
                (
                    interpolationWeightsTFuel, xPosList, yPosList, zPosList,
                    xCell, yCell, zCell,
                    epsilon_[regioni]
                );
            }
            if (isCladTemperatureFieldFromFMU)
            {
                interpolatedValueTClad = Foam::radialBasisFunctionInterpolation::gaussianRadialBasisFunction
                (
                    interpolationWeightsTClad, xPosList, yPosList, zPosList,
                    xCell, yCell, zCell,
                    epsilon_[regioni]
                );
            }
        }
        else if (radialBasisFunctionMethod == "kriging")
        {
            if (isHeatFluxInput)
            {
                interpolatedValueHeatFlux = Foam::radialBasisFunctionInterpolation::kriging
                (
                    krigingWeight_[regioni][i], heatFluxFieldFromFMUList
                );
            }
            else if (isRhoCpdTdtInput)
            {
                interpolatedValueRhoCpdTdt = Foam::radialBasisFunctionInterpolation::kriging
                (
                    krigingWeight_[regioni][i], rhoCpdTdtFieldFromFMUList
                );
            }

            if (isFuelTemperatureFieldFromFMU)
            {
                interpolatedValueTFuel = Foam::radialBasisFunctionInterpolation::kriging
                (
                    krigingWeightFuel_[regioni][i], TFuelFieldFromFMUList
                );
            }
            if (isCladTemperatureFieldFromFMU)
            {
                interpolatedValueTClad = Foam::radialBasisFunctionInterpolation::kriging
                (
                    krigingWeightClad_[regioni][i], TCladFieldFromFMUList
                );
            }
        }
        // else
        // {
        //     interpolatedValue = fractionOfPowerFromNeutronics_[regioni]
        //         * structure_.powerDensityNeutronics()[celli]
        //         * alpha_[celli] / iA_[celli];
        // }
        if (mag(interpolatedValueHeatFlux) > 1e20)
        {
            interpolatedValueHeatFlux = 0.0;
        }
        if (mag(interpolatedValueRhoCpdTdt) > 1e20)
        {
            interpolatedValueRhoCpdTdt = 0.0;
        }

        const scalar powerDensity
        (
            fractionOfPowerFromNeutronics_[regioni]
                * fuelFraction_[regioni]
                * structure_.powerDensityNeutronics()[celli]
        );

        if (isHeatFluxInput && !isRhoCpdTdtInput)
        {
            interpolatedValueRhoCpdTdt = (powerDensity - interpolatedValueHeatFlux * iA_[celli]); // / alpha_[celli];
        }
        if (isRhoCpdTdtInput && !isHeatFluxInput) // rho Cp dT/dt = interpolatedValue
        {
            interpolatedValueHeatFlux = (powerDensity - interpolatedValueRhoCpdTdt /** alpha_[celli]*/) / iA_[celli];
        }

        totalPowerFromHeatFlux += interpolatedValueHeatFlux * iA_[celli] * V[celli];
        totalPowerEnthalpy += interpolatedValueRhoCpdTdt /** alpha_[celli]*/ * V[celli];

        // Correct after the summation
        if (isHeatFluxInput && isRhoCpdTdtInput)
        {
            const scalar correctionFactor
            (
                interpolationCorrectFactor_.oldTime().internalField()[celli]
            );
            interpolatedValueHeatFlux *= correctionFactor;
            interpolatedValueRhoCpdTdt *= correctionFactor;
        }

        // Update Tsurface
        scalar TsurfaceTemp(Tcool + interpolatedValueHeatFlux / HSumi);
        if (0.0 < TsurfaceTemp || TsurfaceTemp < 10000.0)
        {
            Tsurface_[celli] = TsurfaceTemp;
        }
        else if (TsurfaceTemp <= 0.0)
        {
            Tsurface_[celli] = 10.0;

            WarningIn("Foam::nuclearFuelFMU::correctHeatFluxInputsFromFMUs()")
                << nl
                << "    Set Tsurface[" << celli << "] to 10.0 K. "
                << "Tcool = " << Tcool << " K, heatFlux = " << interpolatedValueHeatFlux
                << " W/m2, HSumi = " << HSumi
                << nl
                << endl;
        }
        else
        {
            Tsurface_[celli] = 10000.0;

            WarningIn("Foam::nuclearFuelFMU::correctHeatFluxInputsFromFMUs()")
                << nl
                << "    Set Tsurface[" << celli << "] to 10000.0 K. "
                << "Tcool = " << Tcool << " K, heatFlux = " << interpolatedValueHeatFlux
                << " W/m2, HSumi = " << HSumi
                << nl
                << endl;
        }

        // Update TClad and TFuel
        if (isFuelTemperatureFieldFromFMU)
        {
            Tfav_[celli] = interpolatedValueTFuel;

            // Update average fuel temperature used for coupling
            this->structureRef().TFuelAv()[celli] = Tfav_[celli];
        }
        if (isCladTemperatureFieldFromFMU)
        {
            Tcav_[celli] = interpolatedValueTClad;

            // Update average cladding temperature used for coupling
            this->structureRef().TCladAv()[celli] = Tcav_[celli];
        }
    }

    if (Vtot > 0)
    {
        Info<< "Average temperature in cellZone " << region << ":" << nl
            << "    Tfuel avg = " << fvc::domainIntegrate(Tfav_).value() / Vtot << " K" << nl
            << "    Tclad avg = " << fvc::domainIntegrate(Tcav_).value() / Vtot << " K" << nl
            << endl;
    }

    const scalar oldCorrectionFactor
    (
        interpolationCorrectFactor_.oldTime().internalField()[regionCells[0]]
    );
    Info<< "Applied correction factor in cellZone " << region << ":"
        << oldCorrectionFactor << endl;
    Info<< "Integrated power in cellZone " << region << ":" << nl
        << "    " << region << ":heat flux  = " << totalPowerFromHeatFlux * oldCorrectionFactor << " W" << nl
        << "    " << region << ":rhoCpdTdt  = " << totalPowerEnthalpy * oldCorrectionFactor << " W" << nl
        << "    " << region << ":neutronics = " << totalPowerNeutronics << " W"
        << endl;

    scalar newCorrectionFactor(1.0);

    if (totalPowerFromHeatFlux + totalPowerEnthalpy != 0 && totalPowerNeutronics != 0)
    {
        newCorrectionFactor = totalPowerNeutronics / (totalPowerFromHeatFlux + totalPowerEnthalpy);
    }
    // Cancel the correction if outside range, newCorrectionFactor should always
    // be close to 1
    if (1.0/1.5 > newCorrectionFactor || newCorrectionFactor > 1.5)
    {
        newCorrectionFactor = 1.0;
        Info<< "Correction factor out of range in cellZone " << region
            << ", put to 1.0." << endl;
    }
    forAll(regionCells, i)
    {
        const label celli(regionCells[i]);

        interpolationCorrectFactor_[celli] = newCorrectionFactor;
    }

    Info<< "New correction factor in cellZone " << region << ": "
        << newCorrectionFactor
        << endl;
}


void Foam::powerModels::nuclearFuelFMU::correctInputsForFMUs
(
    const volScalarField& HTSum,  // == SUM_j [htc_j*T_j*frac_j]
    const volScalarField& HSum,   // == SUM_j [htc_j*frac_j]
    const label regioni
)
{
    // Communicating with the FMU, extract default values stored in
    // commDataLayer
    const Time& runTime = this->db().time();
    commDataLayer& data = commDataLayer::New(runTime);

    // Get dictionary of the cellZone power model
    word region(this->toc()[regioni]);
    const dictionary& dict(this->subDict(region));

    const scalarField axialLoc(axialLoc_[regioni]);
    const scalar alphaTstruct(relaxationFactorTstruct_[regioni]);
    const scalar alphaTfluid(relaxationFactorTfluid_[regioni]);
    const scalar alphaHtc(relaxationFactorHtc_[regioni]);
    const scalar alphaProfilePow(relaxationFactorProfilePowerDensity_[regioni]);

    // Update axial locations if necessary
    const word axialLocNameToFMU(dict.get<word>("axialLocationsNameToFMU"));
    word axialLocToFMUtemp(Foam::listConversion::stringify<scalar>(axialLoc));
    word& axialLocToFMU = data.getObj<word>
    (
        axialLocNameToFMU,
        commDataLayer::causality::out
    );
    axialLocToFMU = axialLocToFMUtemp;

    // Get name of the FMI ports for output to the FMUs
    const wordList avgPowerDensityNameToFMU(dict.get<wordList>("avgPowerDensityNameToFMU"));
    const wordList axialProfilePowerDensityNameToFMU(dict.get<wordList>("axialProfilePowerDensityNameToFMU"));
    const wordList TstructNameToFMU(dict.lookupOrDefault<wordList>("TstructNameToFMU", {}));
    const wordList TfluidNameToFMU(dict.lookupOrDefault<wordList>("TfluidNameToFMU", {}));
    const wordList htcNameToFMU(dict.lookupOrDefault<wordList>("htcNameToFMU", {}));

    // Read region values
    // const scalar& fractionOfPowerFromNeutronics(fractionOfPowerFromNeutronics_[regioni]);

    // Loop over the fuel models
    forAll(avgPowerDensityNameToFMU, nameI)
    {
        const scalar& xPos(xPos_[regioni][nameI]);
        const scalar& yPos(yPos_[regioni][nameI]);

        // Create strings for list stringification
        scalarList TstructToFMUtemp(oldTstruct_[regioni][nameI]);
        scalarList TfluidToFMUtemp(oldTfluid_[regioni][nameI]);
        scalarList htcToFMUtemp(oldHtc_[regioni][nameI]);
        scalarList axialProfilePowerDensityToFMUtemp(oldProfilePowerDensity_[regioni][nameI]);

        // Loop over the list of axial locations to extract data
        // Is it necessary to findCell all the time ? Mesh deformation ?
        // Maybe create a labelListListList to store cellNumber for each locations
        scalarField profileData(0);
        forAll(axialLoc, locI)
        {
            const label cellNumber
            (
                mesh_.findCell(point(xPos, yPos, axialLoc[locI]))
            );

            // Accumulate surface temperature values
            TstructToFMUtemp[locI] =
                alphaTstruct * Tsurface_[cellNumber]
                + (1.0 - alphaTstruct) * TstructToFMUtemp[locI];
            TfluidToFMUtemp[locI] =
                alphaTfluid * HTSum[cellNumber] / max(HSum[cellNumber], SMALL)
                + (1.0 - alphaTfluid) * TfluidToFMUtemp[locI];
            htcToFMUtemp[locI] =
                alphaHtc * HSum[cellNumber]
                + (1.0 - alphaHtc) * htcToFMUtemp[locI];

            // Accumulate power values
            const scalar& qRef(structure_.powerDensityNeutronics()[cellNumber]);
            profileData.append(qRef);
        }

        // Rescale data
        scalarInterpolateTableGF zTable(axialLoc, profileData, zMethod_[regioni]);
        const scalar linPowerIntegral(zTable.integral(1));
        const scalar zLength(fuelLength_[regioni]);

        if (linPowerIntegral > 0.0)
        {
            // Normalization
            forAll(profileData, sampleI)
            {
                axialProfilePowerDensityToFMUtemp[sampleI] =
                    alphaProfilePow * profileData[sampleI] / linPowerIntegral
                    + (1.0 - alphaProfilePow) * axialProfilePowerDensityToFMUtemp[sampleI];
            }
        }
        else
        {
            forAll(profileData, sampleI)
            {
                axialProfilePowerDensityToFMUtemp[sampleI] = 1.0;
            }
            WarningIn("Foam::nuclearFuelFMU::correct()") << nl
                << "    Set normalized power density axial profile to 1.0 in "
                << avgPowerDensityNameToFMU[nameI] << "." << nl
                << "    The integral of the power density is 0.0." << nl
                << endl;
        }


        // Extract from commDataLayer to send to FMI
        if (!TstructNameToFMU.empty())
        {
            word& TstructToFMU = data.getObj<word>
            (
                TstructNameToFMU[nameI],
                commDataLayer::causality::out
            );

            // TstructToFMU = TstructToFMUtemp;
            TstructToFMU = Foam::listConversion::stringify<scalar>
            (
                TstructToFMUtemp
            );
        }
        if (!TfluidNameToFMU.empty())
        {
            word& TfluidToFMU = data.getObj<word>
            (
                TfluidNameToFMU[nameI],
                commDataLayer::causality::out
            );

            TfluidToFMU = Foam::listConversion::stringify<scalar>
            (
                TfluidToFMUtemp
            );
        }
        if (!htcNameToFMU.empty())
        {
            word& htcToFMU = data.getObj<word>
            (
                htcNameToFMU[nameI],
                commDataLayer::causality::out
            );

            htcToFMU = Foam::listConversion::stringify<scalar>
            (
                htcToFMUtemp
            );
        }
        scalar& avgPowerDensityToFMU = data.getObj<scalar>
        (
            avgPowerDensityNameToFMU[nameI],
            commDataLayer::causality::out
        );
        word& axialProfilePowerDensityToFMU = data.getObj<word>
        (
            axialProfilePowerDensityNameToFMU[nameI],
            commDataLayer::causality::out
        );

        // Update the values in FMI to FMUs
        avgPowerDensityToFMU = linPowerIntegral / zLength; // correct
        // axialProfilePowerDensityToFMU = axialProfilePowerDensityToFMUtemp;
        axialProfilePowerDensityToFMU = Foam::listConversion::stringify<scalar>
        (
            axialProfilePowerDensityToFMUtemp
        );

        // Update old values to the new one
        oldTstruct_[regioni][nameI] = TstructToFMUtemp;
        oldTfluid_[regioni][nameI] = TfluidToFMUtemp;
        oldHtc_[regioni][nameI] = htcToFMUtemp;
        oldProfilePowerDensity_[regioni][nameI] = axialProfilePowerDensityToFMUtemp;

        //- Write FMI state
        if (mesh_.time().writeTime())
        {
            if (!TstructNameToFMU.empty())
            {
                fmiState_.set
                (
                    TstructNameToFMU[nameI],
                    TstructToFMUtemp
                );
            }
            if (!TfluidNameToFMU.empty())
            {
                fmiState_.set
                (
                    TfluidNameToFMU[nameI],
                    TfluidToFMUtemp
                );
            }
            if (!htcNameToFMU.empty())
            {
                fmiState_.set
                (
                    htcNameToFMU[nameI],
                    htcToFMUtemp
                );
            }
            fmiState_.set
            (
                avgPowerDensityNameToFMU[nameI],
                avgPowerDensityToFMU
            );
            fmiState_.set
            (
                axialProfilePowerDensityNameToFMU[nameI],
                axialProfilePowerDensityToFMUtemp
            );
        }
    }
}


void Foam::powerModels::nuclearFuelFMU::correctT(volScalarField& T) const
{
    Info<< "Correct T nuclearFuelFMU" << endl;
    //- Set T to surface temperature
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        T[celli] = Tsurface_[celli];
    }
}


// * * * * * * * * * * * * * * Private Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::powerModels::nuclearFuelFMU::computeEspilonForRBF
(
    const scalarList list
) const
{
    scalar average(0);
    forAll(list, listI)
    {
        scalar dist(1e13);
        bool isUpdated(false);
        forAll(list, listJ)
        {
            const scalar value(max(list[listI] - list[listJ], list[listJ] - list[listI]));
            if (value > 0)
            {
                dist = min(value, dist);
                isUpdated = true;
            }
        }
        if (isUpdated)
        {
            average += dist;
        }
    }
    return(4.0 * log(1/0.97) / sqr(average / (list.size()-1)));
}


#endif

// ************************************************************************* //
