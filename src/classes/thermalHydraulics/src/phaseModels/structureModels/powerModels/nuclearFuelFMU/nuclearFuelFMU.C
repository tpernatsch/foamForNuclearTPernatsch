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
#include "radialBasisFunctionInterpolation.H"

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
    isHeatFluxInput_(0),
    fractionOfPowerFromNeutronics_(0),
    isFuelTemperatureFieldFromFMU_(0),
    xPos_(0),
    yPos_(0),
    axialLoc_(0),
    fuelLength_(0),
    zMethod_(0),
    radialBasisFunctionMethod_(0),
    epsilon_(0),
    regionIndexToRegionName_(0),
    cellToRegion_(mesh_.cells().size(), 0)
{
    structure_.setRegionField(*this, structureRef.powerDensityNeutronics(), "initPowerDensity");

    this->setInterfacialArea();

    // Communicating with the FMU
    const Time& runTime = this->db().time();
    commDataLayer& data = commDataLayer::New(runTime);

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
        const bool isFuelTemperatureFieldFromFMU
        (
            dict.found("TFuelNameFromFMU") && dict.found("TCladNameFromFMU")
            ? true
            : false
        );
        isFuelTemperatureFieldFromFMU_.append(isFuelTemperatureFieldFromFMU);
        const bool isHeatFluxInput(dict.found("heatFluxNameFromFMU"));
        isHeatFluxInput_.append(isHeatFluxInput);
        xPos_.append(xPos);
        yPos_.append(yPos);
        axialLoc_.append(axialLoc);
        fuelLength_.append(fuelLength);
        zMethod_.append((interpolateTableBase::interpolationMethodNames_[
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
        
        //- Add to regionIndexToRegionName_ mapping
        regionIndexToRegionName_.append(region);

        // Set default values from phaseProperties
        const word axialLocNameToFMU(dict.get<word>("axialLocationsNameToFMU"));
        const wordList avgPowerDensityNameToFMU(dict.get<wordList>("avgPowerDensityNameToFMU"));
        const wordList axialProfilePowerDensityNameToFMU(dict.get<wordList>("axialProfilePowerDensityNameToFMU"));
        const wordList TstructNameToFMU(dict.get<wordList>("TstructNameToFMU"));
        const wordList TFuelNameFromFMU
        (
            isFuelTemperatureFieldFromFMU
            ? dict.get<wordList>("TFuelNameFromFMU")
            : wordList(0)
        );
        const wordList TCladNameFromFMU
        (
            isFuelTemperatureFieldFromFMU
            ? dict.get<wordList>("TCladNameFromFMU")
            : wordList(0)
        );
        wordList heatFluxNameFromFMU(0);
        wordList rhoCpdTdtNameFromFMU(0);

        const scalar initPowerDensityFromDict(dict.get<scalar>("initPowerDensity"));
        const scalar initTemperatureFromDict(dict.get<scalar>("initTstruct"));
        scalar initHeatFluxFromDict(0);
        scalar initRhoCpdTdtFromDict(0);

        if (isHeatFluxInput)
        {
            heatFluxNameFromFMU = dict.get<wordList>("heatFluxNameFromFMU");
            initHeatFluxFromDict = dict.get<scalar>("initHeatFlux");
        }
        else
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
            if (nx != TstructNameToFMU.size())
            {
                FatalErrorInFunction
                    << "xPos and TstructNameToFMU list lengths are different in "
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
            if (!isHeatFluxInput && nx != rhoCpdTdtNameFromFMU.size())
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
            if (isFuelTemperatureFieldFromFMU && nx != TCladNameFromFMU.size())
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
        Info<< endl;
        Info<<"nuclearFuelFMU region " << region << endl;
        forAll(avgPowerDensityNameToFMU, nameI)
        {
            Info<< "    " << xPos[nameI] << " " << yPos[nameI] << " : "
                << avgPowerDensityNameToFMU[nameI] << " "
                << axialProfilePowerDensityNameToFMU[nameI] << " "
                << TstructNameToFMU[nameI] << " "
                << (
                    isHeatFluxInput 
                    ? heatFluxNameFromFMU[nameI]
                    : rhoCpdTdtNameFromFMU[nameI]
                );
                if (isFuelTemperatureFieldFromFMU)
                {
                    Info<< " "
                        << TFuelNameFromFMU[nameI]
                        << " "
                        << TCladNameFromFMU[nameI]
                        << endl;
                }
                else
                {
                    Info<< endl;
                }
    
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
                    scalarList(axialLoc.size(), 1.0)
                )
            );
            data.storeObj
            (
                Foam::listConversion::stringify(initPowerDensity),
                axialProfilePowerDensityNameToFMU[nameI],
                commDataLayer::causality::out
            );

            // Temperature
            scalarList initTstruct
            (
                fmiState_.lookupOrDefault<scalarList>
                (
                    TstructNameToFMU[nameI],
                    scalarList(axialLoc.size(), initTemperatureFromDict)
                )
            );
            data.storeObj
            (
                // initTstruct,
                Foam::listConversion::stringify(initTstruct),
                TstructNameToFMU[nameI],
                commDataLayer::causality::out
            );

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
            else
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
        correctInputsForFMUs(regioni);

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
    label regioni
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
    const bool isFuelTemperatureFieldFromFMU(isFuelTemperatureFieldFromFMU_[regioni]);
    const word radialBasisFunctionMethod(radialBasisFunctionMethod_[regioni]);

    const scalarField axialLoc(axialLoc_[regioni]);

    // Get field name for heat flux reconstruction
    const wordList fieldNameFromFMU(dict.get<wordList>
    (
        isHeatFluxInput ? "heatFluxNameFromFMU" : "rhoCpdTdtNameFromFMU"
    ));

    // Get field name for the fuel pin temperature (how to generalize it for
    // pebble bed?)
    wordList TFuelFieldNameFromFMU(0);
    wordList TCladFieldNameFromFMU(0);
    if (isFuelTemperatureFieldFromFMU)
    {
        TFuelFieldNameFromFMU = dict.get<wordList>("TFuelNameFromFMU");
        TCladFieldNameFromFMU = dict.get<wordList>("TCladNameFromFMU");
    }

    // List of positions and values for the interpolation algorithm
    scalarList xPosList(0);
    scalarList yPosList(0);
    scalarList zPosList(0);
    scalarList fieldFromFMUList(0);
    scalarList TFuelFieldFromFMUList(0);
    scalarList TCladFieldFromFMUList(0);
    
    // Extract x, y, z, heat flux (or enthalpy) into list
    forAll(fieldNameFromFMU, nameI)
    {
        const scalar& xPos(xPos_[regioni][nameI]);
        const scalar& yPos(yPos_[regioni][nameI]);

        // x, y, z positions
        forAll(axialLoc, locI)
        {
            xPosList.append(xPos);
            yPosList.append(yPos);
            zPosList.append(axialLoc[locI]);
        }

        // Extract heat flux or enthalpy
        word& fieldFromFMU = data.getObj<word>
        (
            fieldNameFromFMU[nameI],
            commDataLayer::causality::in
        );
        word& TFuelFieldFromFMU = data.getObj<word>
        (
            isFuelTemperatureFieldFromFMU
            ? TFuelFieldNameFromFMU[nameI]
            : fieldNameFromFMU[nameI],
            commDataLayer::causality::in
        );
        word& TCladFieldFromFMU = data.getObj<word>
        (
            isFuelTemperatureFieldFromFMU
            ? TCladFieldNameFromFMU[nameI]
            : fieldNameFromFMU[nameI],
            commDataLayer::causality::in
        );

        // Split string to list of scalar
        word value;
        scalarList fieldFromFMUListTemp(0);
        scalarList TFuelFieldFromFMUListTemp(0);
        scalarList TCladFieldFromFMUListTemp(0);

        std::stringstream ss(fieldFromFMU);
        while (getline(ss, value, ' '))
        {
            fieldFromFMUListTemp.append(std::stod(value));
        }

        if (isFuelTemperatureFieldFromFMU)
        {
            std::stringstream ssTFuel(TFuelFieldFromFMU);
            while (getline(ssTFuel, value, ' '))
            {
                TFuelFieldFromFMUListTemp.append(std::stod(value));
            }

            std::stringstream ssTClad(TCladFieldFromFMU);
            while (getline(ssTClad, value, ' '))
            {
                TCladFieldFromFMUListTemp.append(std::stod(value));
            }
        }

        // Fill heat flux (or enthalpy) list
        if (fieldFromFMUListTemp.size() == 1)
        {
            forAll(axialLoc, locI)
            {
                fieldFromFMUList.append(fieldFromFMUListTemp[0]);
            }
            if (isFuelTemperatureFieldFromFMU)
            {
                forAll(axialLoc, locI)
                {
                    TFuelFieldFromFMUList.append(TFuelFieldFromFMUListTemp[0]);
                    TCladFieldFromFMUList.append(TCladFieldFromFMUListTemp[0]);
                }
            }
        }
        else if (fieldFromFMUListTemp.size() == axialLoc.size())
        {
            forAll(axialLoc, locI)
            {
                fieldFromFMUList.append(fieldFromFMUListTemp[locI]);
            }
            if (isFuelTemperatureFieldFromFMU)
            {
                forAll(axialLoc, locI)
                {
                    TFuelFieldFromFMUList.append(TFuelFieldFromFMUListTemp[locI]);
                    TCladFieldFromFMUList.append(TCladFieldFromFMUListTemp[locI]);
                }
            }
        }
        else
        {
            FatalErrorInFunction
                << "No heat flux list provided in region " << region 
                << " for FMI port " << fieldNameFromFMU[nameI]
                << exit(FatalError);
        }
        
        if (mesh_.time().writeTime()) // This .write() is slowing down the code (x20)!
        {
            fmiState_.set
            (
                fieldNameFromFMU[nameI], fieldFromFMUListTemp
            );
            if (isFuelTemperatureFieldFromFMU)
            {
                fmiState_.set
                (
                    TFuelFieldNameFromFMU[nameI], TFuelFieldFromFMUListTemp
                );
                fmiState_.set
                (
                    TCladFieldNameFromFMU[nameI], TCladFieldFromFMUListTemp
                );
            }
        }
    }

    // Generate weight for RBF
    scalarList interpolationWeights(0);
    scalarList interpolationWeightsTFuel(0);
    scalarList interpolationWeightsTClad(0);
    if (radialBasisFunctionMethod == "polyharmonicSpline")
    {
        // if (!isSteadyStateMode_)
        // {
            interpolationWeights = Foam::radialBasisFunctionInterpolation::solvePolyharmonicSpline
            (
                xPosList, yPosList, zPosList, fieldFromFMUList, invRBFmatrix_
            );

            if (isFuelTemperatureFieldFromFMU)
            {
                interpolationWeightsTFuel = Foam::radialBasisFunctionInterpolation::solvePolyharmonicSpline
                (
                    xPosList, yPosList, zPosList, TFuelFieldFromFMUList, invRBFmatrix_
                );
                interpolationWeightsTClad = Foam::radialBasisFunctionInterpolation::solvePolyharmonicSpline
                (
                    xPosList, yPosList, zPosList, TCladFieldFromFMUList, invRBFmatrix_
                );
            }
        // }
        // else
        // {
        //     scalar totalPowerNeutronics(0); // Divided by the volumetric surface
        //     forAll(regionCells, celli)
        //     {
        //         totalPowerNeutronics += fractionOfPowerFromNeutronics_[regioni]
        //             * structure_.powerDensityNeutronics()[celli]
        //             * V[celli] / iA_[celli] * alpha_[celli];
        //     }

        //     interpolationWeights = Foam::radialBasisFunctionInterpolation::solvePolyharmonicSplineIntegral
        //     (
        //         xPosList, yPosList, zPosList, fieldFromFMUList, totalPowerNeutronics, regionCells, mesh_, invRBFmatrix_
        //     );
        // }
    }
    else if (radialBasisFunctionMethod == "gaussian")
    {
        interpolationWeights = Foam::radialBasisFunctionInterpolation::solveGaussianRadialBasisFunction
        (
            xPosList, yPosList, zPosList, fieldFromFMUList, epsilon_[regioni], invRBFmatrix_
        );
        
        if (isFuelTemperatureFieldFromFMU)
        {
            interpolationWeightsTFuel = Foam::radialBasisFunctionInterpolation::solveGaussianRadialBasisFunction
            (
                xPosList, yPosList, zPosList, TFuelFieldFromFMUList, epsilon_[regioni], invRBFmatrix_
            );
            interpolationWeightsTClad = Foam::radialBasisFunctionInterpolation::solveGaussianRadialBasisFunction
            (
                xPosList, yPosList, zPosList, TCladFieldFromFMUList, epsilon_[regioni], invRBFmatrix_
            );
        }
    }
    else if (radialBasisFunctionMethod == "kriging")
    {
        Foam::radialBasisFunctionInterpolation::solveKriging
        (
            xPosList, yPosList, zPosList, invRBFmatrix_
        );
    }
    else
    {
        FatalErrorInFunction
            << radialBasisFunctionMethod << " is an incorrect "
            << "radial basis function method provided in region " << region
            << ". Available methods: polyharmonicSpline, gaussian, kriging"
            << exit(FatalError);
    }

    scalar totalPowerFromHeatFlux(0), totalPowerEnthalpy(0), totalPowerNeutronics(0);

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
        scalar interpolatedValue(0);
        scalar interpolatedValueTFuel(0);
        scalar interpolatedValueTClad(0);
        if (radialBasisFunctionMethod == "polyharmonicSpline")
        {
            // if (!isSteadyStateMode_)
            // {
                interpolatedValue = Foam::radialBasisFunctionInterpolation::polyharmonicSpline
                (
                    interpolationWeights, 
                    xPosList, yPosList, zPosList, 
                    xCell, yCell, zCell
                );
                
                if (isFuelTemperatureFieldFromFMU)
                {
                    interpolatedValueTFuel = Foam::radialBasisFunctionInterpolation::polyharmonicSpline
                    (
                        interpolationWeightsTFuel, 
                        xPosList, yPosList, zPosList, 
                        xCell, yCell, zCell
                    );
                    interpolatedValueTClad = Foam::radialBasisFunctionInterpolation::polyharmonicSpline
                    (
                        interpolationWeightsTClad, 
                        xPosList, yPosList, zPosList, 
                        xCell, yCell, zCell
                    );
                }
            // }
            // else
            // {
            //     interpolatedValue = Foam::radialBasisFunctionInterpolation::polyharmonicSplineIntegral
            //     (
            //         interpolationWeights, 
            //         xPosList, yPosList, zPosList, 
            //         xCell, yCell, zCell
            //     );
            // }
        }
        else if (radialBasisFunctionMethod == "gaussian")
        {
            interpolatedValue = Foam::radialBasisFunctionInterpolation::gaussianRadialBasisFunction
            (
                interpolationWeights, 
                xPosList, yPosList, zPosList, 
                xCell, yCell, zCell, 
                epsilon_[regioni]
            );
            
            if (isFuelTemperatureFieldFromFMU)
            {
                interpolatedValueTFuel = Foam::radialBasisFunctionInterpolation::gaussianRadialBasisFunction
                (
                    interpolationWeightsTFuel, 
                    xPosList, yPosList, zPosList, 
                    xCell, yCell, zCell, 
                    epsilon_[regioni]
                );
                interpolatedValueTClad = Foam::radialBasisFunctionInterpolation::gaussianRadialBasisFunction
                (
                    interpolationWeightsTClad, 
                    xPosList, yPosList, zPosList, 
                    xCell, yCell, zCell, 
                    epsilon_[regioni]
                );
            }
        }
        else if (radialBasisFunctionMethod == "kriging")
        {
            interpolatedValue = Foam::radialBasisFunctionInterpolation::kriging
            (
                xPosList, yPosList, zPosList, fieldFromFMUList,
                xCell, yCell, zCell, 
                invRBFmatrix_
            );
            
            if (isFuelTemperatureFieldFromFMU)
            {
                interpolatedValueTFuel = Foam::radialBasisFunctionInterpolation::kriging
                (
                    xPosList, yPosList, zPosList, TFuelFieldFromFMUList,
                    xCell, yCell, zCell, 
                    invRBFmatrix_
                );
                interpolatedValueTClad = Foam::radialBasisFunctionInterpolation::kriging
                (
                    xPosList, yPosList, zPosList, TCladFieldFromFMUList,
                    xCell, yCell, zCell, 
                    invRBFmatrix_
                );
            }
        }
        // else
        // {
        //     interpolatedValue = fractionOfPowerFromNeutronics_[regioni]
        //         * structure_.powerDensityNeutronics()[celli]
        //         * alpha_[celli] / iA_[celli];
        // }

        scalar heatFlux(0);
        if (isHeatFluxInput)
        {
            heatFlux = interpolatedValue;
        }
        else // rho Cp dT/dt
        {
            const scalar powerDensity
            (
                fractionOfPowerFromNeutronics_[regioni]
                    * structure_.powerDensityNeutronics()[celli]
            );
            const scalar rhoCpdTdt(interpolatedValue);

            heatFlux = alpha_[celli] * (powerDensity - rhoCpdTdt) / iA_[celli];
            
            totalPowerEnthalpy += rhoCpdTdt * V[celli];
            totalPowerNeutronics += powerDensity * V[celli];
        }

        totalPowerFromHeatFlux += heatFlux * iA_[celli] * V[celli] / alpha_[celli];


        // Update Tsurface
        Tsurface_[celli] = Tcool + heatFlux / HSumi;

        // Update TClad and TFuel
        if (isFuelTemperatureFieldFromFMU)
        {
            Tcav_[celli] = interpolatedValueTFuel;
            Tfav_[celli] = interpolatedValueTClad;

            // Update average fuel and clad temp used for coupling
            this->structureRef().TFuelAv()[celli] = Tfav_[celli];
            this->structureRef().TCladAv()[celli] = Tcav_[celli];
        }
    }

    Info<< "Integrated power in cellZone " << region << ":" << nl 
        << "    heat flux  = " << totalPowerFromHeatFlux << " W" << nl
        << "    rhoCpdTdt  = " << totalPowerEnthalpy << " W" << nl
        << "    neutronics = " << totalPowerNeutronics << " W"
        << endl;
}

void Foam::powerModels::nuclearFuelFMU::correctInputsForFMUs(label regioni)
{
    // Communicating with the FMU, extract default values stored in 
    // commDataLayer 
    const Time& runTime = this->db().time();
    commDataLayer& data = commDataLayer::New(runTime);

    // Get dictionary of the cellZone power model
    word region(this->toc()[regioni]);
    const dictionary& dict(this->subDict(region));

    const scalarField axialLoc(axialLoc_[regioni]);

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
    const wordList TstructNameToFMU(dict.get<wordList>("TstructNameToFMU"));

    // Read region values
    const scalar& fractionOfPowerFromNeutronics(fractionOfPowerFromNeutronics_[regioni]);

    // Loop over the fuel models
    forAll(avgPowerDensityNameToFMU, nameI)
    {
        const scalar& xPos(xPos_[regioni][nameI]);
        const scalar& yPos(yPos_[regioni][nameI]);

        // Create strings for list stringification
        // word TstructToFMUtemp("");
        // word axialProfilePowerDensityToFMUtemp("");
        scalarList TstructToFMUtemp(0);
        scalarList axialProfilePowerDensityToFMUtemp(0);

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
            // TstructToFMUtemp += std::to_string(Tsurface_[cellNumber])+" ";
            TstructToFMUtemp.append(Tsurface_[cellNumber]);

            // Accumulate power values
            const scalar& qRef(structure_.powerDensityNeutronics()[cellNumber]);
            profileData.append(qRef * fractionOfPowerFromNeutronics);
        }

        // Rescale data
        scalarInterpolateTable zTable(axialLoc, profileData, zMethod_[regioni]);
        const scalar linPowerIntegral(zTable.integral(1));
        const scalar zLength(fuelLength_[regioni]); // axialLoc[axialLoc.size()-1] - axialLoc[0]);

        if (linPowerIntegral > 0.0)
        {
            // Normalization
            forAll(profileData, sampleI)
            {
                // axialProfilePowerDensityToFMUtemp += std::to_string
                // (
                //     profileData[sampleI] * zLength / linPowerIntegral
                // )+" ";
                axialProfilePowerDensityToFMUtemp.append
                (
                    // profileData[sampleI] * zLength / linPowerIntegral
                    profileData[sampleI] / linPowerIntegral
                );
            }
        }
        else
        {
            forAll(profileData, sampleI)
            {
                // axialProfilePowerDensityToFMUtemp += "1.0 ";
                axialProfilePowerDensityToFMUtemp.append(1.0);
            }
            WarningIn("Foam::nuclearFuelFMU::correct()") << nl
                << "    Set normalized power density axial profile to 1.0 in " 
                << avgPowerDensityNameToFMU[nameI] << "." << nl
                << "    The integral of the power density is 0.0." << nl
                << endl;
        }


        // Extract from commDataLayer to send to FMI
        word& TstructToFMU = data.getObj<word>
        (
            TstructNameToFMU[nameI], 
            commDataLayer::causality::out
        );
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
        // TstructToFMU = TstructToFMUtemp;
        TstructToFMU = Foam::listConversion::stringify<scalar>
        (
            TstructToFMUtemp
        );
        avgPowerDensityToFMU = linPowerIntegral / zLength; // correct
        // axialProfilePowerDensityToFMU = axialProfilePowerDensityToFMUtemp;
        axialProfilePowerDensityToFMU = Foam::listConversion::stringify<scalar>
        (
            axialProfilePowerDensityToFMUtemp
        );

        //- Write FMI state
        if (mesh_.time().writeTime())
        {
            fmiState_.set
            (
                TstructNameToFMU[nameI],
                TstructToFMUtemp
            );
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
