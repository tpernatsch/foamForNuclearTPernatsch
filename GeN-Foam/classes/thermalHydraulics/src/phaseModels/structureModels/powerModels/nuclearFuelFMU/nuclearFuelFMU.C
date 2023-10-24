/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2306                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2023 OpenCFD Ltd.         |
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
    /*Tfav_
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
    ),*/
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
    fractionOfPowerFromNeutronics_(0),
    xPos_(0),
    yPos_(0),
    axialLoc_(0),
    zMethod_(0),
    radialBasisFunctionMethod_(0),
    epsilon_(0),
    regionIndexToRegionName_(0),
    cellToRegion_(mesh_.cells().size(), 0)
{
    structure_.setRegionField(*this, structureRef.powerDensityNeutronics(), "initPowerDensity");

    this->setInterfacialArea();

    forAll(this->toc(), regioni)
    {
        word region(this->toc()[regioni]);
        const dictionary& dict(this->subDict(region));

        // Setup cellToRegion_ mapping
        const labelList& regionCells(structure_.cellLists()[region]);
        forAll(regionCells, i)
        {
            label celli(regionCells[i]);
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
        const word radialBasisFunctionMethod(dict.get<word>("radialBasisFunctionMethod"));

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
                << "None of xPos/yPos or xyPosLattice have been provided in region " << region
                << ". It is mandatory to add XY-locations samples."
                << exit(FatalError);
        }

        // Check lengths
        const label nx(xPos.size());
        if (nx != yPos.size())
        {
            FatalErrorInFunction
                << "xPos and yPos lengths are different in region " << region
                << exit(FatalError);
        }

        // Fill in lists for this region
        fractionOfPowerFromNeutronics_.append(fractionOfPowerFromNeutronics);
        xPos_.append(xPos);
        yPos_.append(yPos);
        axialLoc_.append(axialLoc);
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

        // Communicating with the FMU, set default values from phaseProperties
        const Time& runTime = this->db().time();
        commDataLayer& data = commDataLayer::New(runTime);
        const word axialLocNameToFMU(dict.get<word>("axialLocationsNameToFMU"));
        const wordList avgPowerDensityNameToFMU(dict.get<wordList>("avgPowerDensityNameToFMU"));
        const wordList axialProfilePowerDensityNameToFMU(dict.get<wordList>("axialProfilePowerDensityNameToFMU"));
        const wordList TstructNameToFMU(dict.get<wordList>("TstructNameToFMU"));
        const wordList heatFluxNameFromFMU(dict.get<wordList>("heatFluxNameFromFMU"));

        const scalar initPowerDensityFromDict(dict.get<scalar>("initPowerDensity"));
        const scalar initTemperatureFromDict(dict.get<scalar>("initTstruct"));
        const scalar initHeatFluxFromDict(dict.get<scalar>("initHeatFlux"));

        // Check lengths of FMI port lists
        {
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
            if (nx != heatFluxNameFromFMU.size())
            {
                FatalErrorInFunction
                    << "xPos and heatFluxNameFromFMU list lengths are different "
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

        forAll(avgPowerDensityNameToFMU, nameI)
        {
            // Average power density
            data.storeObj
            (
                initPowerDensityFromDict,
                avgPowerDensityNameToFMU[nameI],
                commDataLayer::causality::out
            );

            // Axial distribution of power density
            word initPowerDensity("");
            forAll(axialLoc, locI)
            {
                initPowerDensity += "1.0 ";
            }
            data.storeObj
            (
                initPowerDensity,
                axialProfilePowerDensityNameToFMU[nameI],
                commDataLayer::causality::out
            );

            // Temperature
            word initTstruct("");
            forAll(axialLoc, locI)
            {
                initTstruct += std::to_string(initTemperatureFromDict) + " ";
            }
            data.storeObj
            (
                initTstruct,
                TstructNameToFMU[nameI],
                commDataLayer::causality::out
            );

            // Heat flux
            word initHeatFlux("");
            forAll(axialLoc, locI)
            {
                initHeatFlux += std::to_string(initHeatFluxFromDict) + " ";
            }
            data.storeObj
            (
                initHeatFlux,
                heatFluxNameFromFMU[nameI],
                commDataLayer::causality::in
            );
        }
    }
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
    // Communicating with the FMU, extract default values stored in 
    // commDataLayer 
    const Time& runTime = this->db().time();
    commDataLayer& data = commDataLayer::New(runTime);

    const scalarList& V(mesh_.V());

    // Update surface temperature from FMU using heat flux from FMUs
    forAll(this->toc(), regioni)
    {
        word region(this->toc()[regioni]);
        const labelList& regionCells(structure_.cellLists()[region]);
        const dictionary& dict(this->subDict(region));

        const scalarField axialLoc(axialLoc_[regioni]);

        const wordList heatFluxNameFromFMU(dict.get<wordList>("heatFluxNameFromFMU"));

        // const vectorField gradPowerDensity(fvc::grad(structure_.powerDensityNeutronics()));

        scalarList xPosList(0);
        scalarList yPosList(0);
        scalarList zPosList(0);
        scalarList heatFluxList(0);
        // scalarList gradxPowerDensity(0);
        // scalarList gradyPowerDensity(0);
        
        // Extract x, y, z, heat flux into list
        forAll(heatFluxNameFromFMU, nameI)
        {
            const scalar& xPos(xPos_[regioni][nameI]);
            const scalar& yPos(yPos_[regioni][nameI]);

            // x, y, z positions
            forAll(axialLoc, locI)
            {
                xPosList.append(xPos);
                yPosList.append(yPos);
                zPosList.append(axialLoc[locI]);

                // label cellNumber = mesh_.findCell(point(xPos, yPos, axialLoc[locI]));
                
                // gradxPowerDensity.append(0.217334 / 114.0869 * gradPowerDensity[cellNumber].x());
                // gradyPowerDensity.append(0.217334 / 114.0869 * gradPowerDensity[cellNumber].y());
            }

            // Extract heat flux
            word& heatFluxFromFMU = data.getObj<word>
            (
                heatFluxNameFromFMU[nameI],
                commDataLayer::causality::in
            );

            // Split string to list of scalar
            word value;
            std::stringstream ss(heatFluxFromFMU);
            scalarList heatFluxListTemp(0);
            while (getline(ss, value, ' '))
            {
                heatFluxListTemp.append(std::stod(value));
            }

            // Fill heat flux list
            if (heatFluxListTemp.size() > 0)
            {
                if (heatFluxListTemp.size() == axialLoc.size())
                {
                    forAll(axialLoc, locI)
                    {
                        heatFluxList.append(heatFluxListTemp[locI]);
                    }
                }
                else
                {
                    forAll(axialLoc, locI)
                    {
                        heatFluxList.append(heatFluxListTemp[0]);
                    }
                }
            }
            else
            {
                FatalErrorInFunction
                    << "No heat flux list provided in region " << region 
                    << " for FMI port " << heatFluxNameFromFMU[nameI]
                    << exit(FatalError);
            }
        }

        // Generate weight for RBF
        scalarList w(0); 
        if (radialBasisFunctionMethod_[regioni] == "polyharmonicSpline")
        {
            w = Foam::radialBasisFunctionInterpolation::solvePolyharmonicSpline
            (
                xPosList, yPosList, zPosList, heatFluxList, invRBFmatrix_
            );
        }
        else if (radialBasisFunctionMethod_[regioni] == "gaussian")
        {
            w = Foam::radialBasisFunctionInterpolation::solveGaussianRadialBasisFunction
            (
                xPosList, yPosList, zPosList, heatFluxList, epsilon_[regioni], invRBFmatrix_
            );
        }
        else
        {
            FatalErrorInFunction
                << radialBasisFunctionMethod_[regioni] << " is an incorrect "
                << "radial basis function method provided in region " << region
                << ". Available methods: polyharmonicSpline, gaussian"
                << exit(FatalError);
        }
        /*
        const scalarList w = solvePolyharmonicSplineDerivative
        (
            xPosList, yPosList, zPosList, heatFluxList, gradxPowerDensity, gradyPowerDensity
        );
        */
        /*
        const labelList& regionCells(structure_.cellLists()[region]);
        const scalarList w = solvePolyharmonicSplineIntegral
        (
            xPosList, yPosList, zPosList, heatFluxList, 300e6 / 114.08690675627066, regionCells, mesh_
        );
        */

        scalar totalPowerFromHeatFlux(0);

        // Update all Tsurface cells in the region using heat flux interpolated 
        forAll(regionCells, i)
        {
            label celli(regionCells[i]);
            
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
            /*
            polyharmonicSplineDerivative
            (
                w, 
                xPosList, yPosList, zPosList, 
                xCell, yCell, zCell
            )
            */
            /*
            polyharmonicSplineIntegral
            (
                w, 
                xPosList, yPosList, zPosList, 
                xCell, yCell, zCell
            )
            */

            scalar heatFlux(0); 
            if (radialBasisFunctionMethod_[regioni] == "polyharmonicSpline")
            {
                heatFlux = Foam::radialBasisFunctionInterpolation::polyharmonicSpline
                (
                    w, 
                    xPosList, yPosList, zPosList, 
                    xCell, yCell, zCell
                );
            }
            else if (radialBasisFunctionMethod_[regioni] == "gaussian")
            {
                heatFlux = Foam::radialBasisFunctionInterpolation::gaussianRadialBasisFunction
                (
                    w, 
                    xPosList, yPosList, zPosList, 
                    xCell, yCell, zCell, 
                    epsilon_[regioni]
                );
            }

            totalPowerFromHeatFlux += heatFlux * iA_[celli] * V[celli] / alpha_[celli];
                
            // Update Tsurface
            Tsurface_[celli] = Tcool + heatFlux / HSumi;
        }

        Info<< "Integrated heat flux in cellZone " << region << " = " 
            << totalPowerFromHeatFlux << " W"
            << endl;
    }


    // Update surface temperature and power density sample lines for FMU
    forAll(this->toc(), regioni)
    {
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
            word TstructToFMUtemp("");
            word axialProfilePowerDensityToFMUtemp("");

            // Loop over the list of FMI temperature ports to extract data
            scalarField profileData(0);
            forAll(axialLoc, locI)
            {
                const label cellNumber
                (
                    mesh_.findCell(point(xPos, yPos, axialLoc[locI]))
                );
                
                // Accumulate surface temperature values
                TstructToFMUtemp += std::to_string(Tsurface_[cellNumber])+" ";
                // TstructToFMUtemp += std::to_string(HTSum[cellNumber] / HSum[cellNumber])+" ";

                // Accumulate power values
                const scalar& qRef(structure_.powerDensityNeutronics()[cellNumber]);
                profileData.append(qRef * fractionOfPowerFromNeutronics);
            }

            // Rescale data
            scalarInterpolateTable zTable(axialLoc, profileData, zMethod_[regioni]);
            const scalar linPowerIntegral(zTable.integral(1));

            if (linPowerIntegral > 0.0)
            {
                forAll(profileData, sampleI)
                {
                    axialProfilePowerDensityToFMUtemp += std::to_string
                    (
                        profileData[sampleI] / linPowerIntegral
                    )+" ";
                }
            }
            else
            {
                forAll(profileData, sampleI)
                {
                    axialProfilePowerDensityToFMUtemp += "1.0 ";
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

            const scalar zLength(axialLoc[axialLoc.size()-1] - axialLoc[0]);

            // Update the values in FMI
            TstructToFMU = TstructToFMUtemp;
            avgPowerDensityToFMU = linPowerIntegral / zLength;
            axialProfilePowerDensityToFMU = axialProfilePowerDensityToFMUtemp;
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
