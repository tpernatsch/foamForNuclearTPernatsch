/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2212                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2022 OpenCFD Ltd.         |
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
#include "commDataLayer.H"

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

        // Fill in lists for this region
        fractionOfPowerFromNeutronics_.append(fractionOfPowerFromNeutronics);
        xPos_.append(dict.get<scalarList>("xPos"));
        yPos_.append(dict.get<scalarList>("yPos"));
        axialLoc_.append(dict.get<scalarField>("axialLocations"));
        epsilon_.append(dict.get<scalarList>("epsilon"));
        zMethod_.append((interpolateTableBase::interpolationMethodNames_[
            dict.lookupOrDefault<word>("axialInterpolationMethod", "linear")
        ]));

        if (xPos_[-1].size() != xPos_[-1].size())
        {
            FatalErrorInFunction
                << "xPos and yPos lengths are different in region " << region
                << exit(FatalError);
        }

        {
            const scalar epsilonX = computeEspilonForRBF(xPos_[xPos_.size()-1]);
            const scalar epsilonY = computeEspilonForRBF(yPos_[yPos_.size()-1]);
            const scalar epsilonZ = computeEspilonForRBF(axialLoc_[axialLoc_.size()-1]);
            Info<< "Recommended RBF epsilon parameters in cellZone " << region << nl
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
        const word axialLocNameToFMU = dict.get<word>("axialLocationsNameToFMU");
        const wordList avgPowerDensityNameToFMU = dict.get<wordList>("avgPowerDensityNameToFMU");
        const wordList axialProfilePowerDensityNameToFMU = dict.get<wordList>("axialProfilePowerDensityNameToFMU");
        const wordList TstructNameToFMU = dict.get<wordList>("TstructNameToFMU");
        const wordList surfacePowerNameFromFMU = dict.get<wordList>("surfacePowerNameFromFMU");

        const scalar initPowerFromDict(dict.get<scalar>("initPowerDensity"));
        const scalar initTemperatureFromDict(dict.get<scalar>("initTstruct"));
        const scalar initSurfacePowerFromDict(dict.get<scalar>("initSurfacePower"));
        
        const scalarField axialLoc(axialLoc_[regioni]);

        // Axial location
        word axialLocStringified(stringify(axialLoc));
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
                initPowerFromDict,
                avgPowerDensityNameToFMU[nameI],
                commDataLayer::causality::out
            );

            // Axial distribution of power density
            word initPowerDensity("");
            forAll(axialLoc, locI)
            {
                initPowerDensity += "1.0 "; // std::to_string(initPowerFromDict) + " ";
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

            // Surface power
            word initSurfacePower("");
            forAll(axialLoc, locI)
            {
                initSurfacePower += std::to_string(initSurfacePowerFromDict) + " ";
            }
            data.storeObj
            (
                initSurfacePower,
                surfacePowerNameFromFMU[nameI],
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

    // Update surface temperature from FMU using surface power
    forAll(this->toc(), regioni)
    {
        word region(this->toc()[regioni]);
        const dictionary& dict(this->subDict(region));

        const scalarField axialLoc(axialLoc_[regioni]);

        const wordList surfacePowerNameFromFMU = dict.get<wordList>("surfacePowerNameFromFMU");

        scalarList xPosList(0);
        scalarList yPosList(0);
        scalarList zPosList(0);
        scalarList surfacePowerList(0);
        
        // Extract x, y, z, surface power into list
        forAll(surfacePowerNameFromFMU, nameI)
        {
            // x, y, z positions
            forAll(axialLoc, locI)
            {
                xPosList.append(xPos_[regioni][nameI]);
                yPosList.append(yPos_[regioni][nameI]);
                zPosList.append(axialLoc[locI]);
            }

            // Extract surface power
            word& surfacePowerFromFMU = data.getObj<word>
            (
                surfacePowerNameFromFMU[nameI],
                commDataLayer::causality::in
            );

            // Split string to list of scalar
            word value;
            std::stringstream ss(surfacePowerFromFMU);
            scalarList surfacePowerListTemp(0);
            while (getline(ss, value, ' '))
            {
                surfacePowerListTemp.append(std::stod(value));
            }

            // Fill surface power list
            if (surfacePowerListTemp.size() > 0)
            {
                if (surfacePowerListTemp.size() == axialLoc.size())
                {
                    forAll(axialLoc, locI)
                    {
                        surfacePowerList.append(surfacePowerListTemp[locI]);
                    }
                }
                else
                {
                    forAll(axialLoc, locI)
                    {
                        surfacePowerList.append(surfacePowerListTemp[0]);
                    }
                }
            }
            else
            {
                FatalErrorInFunction
                    << "No surface power list provided in region " << region
                    << " by FMU"
                    << exit(FatalError);
            }
        }

        // Generate weight for RBF
        const scalarList w = solveRadialBasisFunction
        (
            xPosList, yPosList, zPosList, surfacePowerList, epsilon_[regioni]
        );

        const labelList& regionCells(structure_.cellLists()[region]);
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

            // Interpolate surface power at x, y, z position
            const scalar surfacePower
            (
                radialBasisFucntion
                (
                    w, 
                    xPosList, yPosList, zPosList, 
                    xCell, yCell, zCell, 
                    epsilon_[regioni]
                )
            );

            // Update Tstruct
            Tsurface_[celli] = Tcool + surfacePower / HSumi;
        }
    }

    // Update surface temperature and power density sample lines for FMU
    forAll(this->toc(), regioni)
    {
        word region(this->toc()[regioni]);
        const dictionary& dict(this->subDict(region));

        const scalarField axialLoc(axialLoc_[regioni]);

        // Update axial locations if necessary
        const word axialLocNameToFMU = dict.get<word>("axialLocationsNameToFMU");
        word axialLocToFMUtemp(stringify(axialLoc));
        word& axialLocToFMU = data.getObj<word>
        (
            axialLocNameToFMU, 
            commDataLayer::causality::out
        );
        axialLocToFMU = axialLocToFMUtemp;


        const wordList avgPowerDensityNameToFMU = dict.get<wordList>("avgPowerDensityNameToFMU");
        const wordList axialProfilePowerDensityNameToFMU = dict.get<wordList>("axialProfilePowerDensityNameToFMU");
        const wordList TstructNameToFMU = dict.get<wordList>("TstructNameToFMU");

        // Read region values
        const scalar& fractionOfPowerFromNeutronics(fractionOfPowerFromNeutronics_[regioni]);

        // Loop over the pins
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
                label cellNumber = mesh_.findCell(point(xPos, yPos, axialLoc[locI]));
                
                // Accumulate surface temperature values
                TstructToFMUtemp += std::to_string(Tsurface_[cellNumber])+" ";
                // TstructToFMUtemp += std::to_string(HTSum[cellNumber] / HSum[cellNumber])+" ";

                // Accumulate power values
                const scalar& qRef(structure_.powerDensityNeutronics()[cellNumber]);
                // axialProfilePowerDensityToFMUtemp += std::to_string(qRef * fractionOfPowerFromNeutronics)+" ";
                profileData.append(qRef * fractionOfPowerFromNeutronics);
            }

            // Rescale data
            scalarInterpolateTable zTable(axialLoc, profileData, zMethod_[regioni]);
            const scalar linPowerIntegral = zTable.integral(1);

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
                    << "    Set normalized power density axial profile to 1.0." << nl
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

Foam::scalarList Foam::powerModels::nuclearFuelFMU::solveRadialBasisFunction
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalarList vList,
    const scalarList eps
) const
{
    SquareMatrix<scalar> A(xList.size(), 0.0);
    scalarList w(xList.size(), 0.0);

    forAll(xList, i)
    {
        forAll(xList, j)
        {
            A[i][j] = exp(
                -eps[0] * sqr(xList[i]-xList[j]) 
                -eps[1] * sqr(yList[i]-yList[j])
                -eps[2] * sqr(zList[i]-zList[j])
            );
        }
    }

    solve(w, A, vList);

    return(w);
}


Foam::scalar Foam::powerModels::nuclearFuelFMU::radialBasisFucntion
(
    const scalarList w,
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalar x,
    const scalar y,
    const scalar z,
    const scalarList eps
) const
{
    scalar res(0);
    forAll(w, i)
    {
        res += w[i] * exp(
            -eps[0] * sqr(x-xList[i])
            -eps[1] * sqr(y-yList[i]) 
            -eps[2] * sqr(z-zList[i])
        );
    }
    return(res);
}


Foam::word Foam::powerModels::nuclearFuelFMU::stringify
(
    const scalarList list
) const
{
    word wordTemp;
    forAll(list, i)
    {
        wordTemp += std::to_string(list[i]) + " ";
    }
    return(wordTemp);
}


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
