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
        const scalarList xPos(dict.get<scalarList>("xPos"));
        const scalarList yPos(dict.get<scalarList>("yPos"));
        const scalarField axialLoc(dict.get<scalarField>("axialLocations"));
        const word radialBasisFunctionMethod(dict.get<word>("radialBasisFunctionMethod"));

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
                    << "xPos and avgPowerDensityNameToFMU lengths are "
                    << "different in region " << region
                    << exit(FatalError);
            }
            if (nx != axialProfilePowerDensityNameToFMU.size())
            {
                FatalErrorInFunction
                    << "xPos and axialProfilePowerDensityNameToFMU lengths are "
                    << "different in region " << region
                    << exit(FatalError);
            }
            if (nx != TstructNameToFMU.size())
            {
                FatalErrorInFunction
                    << "xPos and TstructNameToFMU lengths are different in "
                    << "region " << region
                    << exit(FatalError);
            }
            if (nx != heatFluxNameFromFMU.size())
            {
                FatalErrorInFunction
                    << "xPos and heatFluxNameFromFMU lengths are different "
                    << "in region " << region
                    << exit(FatalError);
            }
        }

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
                
                // gradxPowerDensity.append(gradPowerDensity[cellNumber].x());
                // gradyPowerDensity.append(gradPowerDensity[cellNumber].y());
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
                    << " by FMU"
                    << exit(FatalError);
            }
        }

        // Generate weight for RBF
        scalarList w(0); 
        if (radialBasisFunctionMethod_[regioni] == "polyharmonicSpline")
        {
            w = solvePolyharmonicSpline
            (
                xPosList, yPosList, zPosList, heatFluxList
            );
        }
        else if (radialBasisFunctionMethod_[regioni] == "gaussian")
        {
            w = solveGaussianRadialBasisFunction
            (
                xPosList, yPosList, zPosList, heatFluxList, epsilon_[regioni]
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
        const scalarList w = solvePolyharmonicSplineIntegral
        (
            xPosList, yPosList, zPosList, heatFluxList, region
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
                heatFlux = polyharmonicSpline
                (
                    w, 
                    xPosList, yPosList, zPosList, 
                    xCell, yCell, zCell
                );
            }
            else if (radialBasisFunctionMethod_[regioni] == "gaussian")
            {
                heatFlux = gaussianRadialBasisFunction
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
        word axialLocToFMUtemp(stringify(axialLoc));
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

Foam::scalarList Foam::powerModels::nuclearFuelFMU::solveGaussianRadialBasisFunction
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalarList vList,
    const scalarList eps
)
{
    Info << "Use Gaussian radial basis function" << endl;
    
    const label nx(xList.size());

    scalarList w(nx, 0.0);

    if (invRBFmatrix_.m() != nx)
    {
        invRBFmatrix_.resize(nx);

        SquareMatrix<scalar> A(nx, 0.0);

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

        LUscalarMatrix Atemp(A);

        Info<< "Inverse RBF matrix" << endl;
        Atemp.inv(invRBFmatrix_);

        // solve(w, A, vList);
    }

    w = invRBFmatrix_ * vList;

    return(w);
}


Foam::scalar Foam::powerModels::nuclearFuelFMU::gaussianRadialBasisFunction
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


Foam::scalarList Foam::powerModels::nuclearFuelFMU::solvePolyharmonicSpline
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalarList vList
)
{
    Info << "Use polyharmonic spline RBF" << endl;
    const label nx(xList.size());
    scalarList w(nx+4, 0.0);

    if (invRBFmatrix_.m() != nx+4)
    {
        invRBFmatrix_.resize(nx+4);

        SquareMatrix<scalar> A(nx+4, 0.0);

        forAll(xList, i)
        {
            forAll(xList, j)
            {
                if (i != j)
                {
                    const scalar r(sqrt(
                        sqr(xList[i]-xList[j]) 
                        + sqr(yList[i]-yList[j])
                        + sqr(zList[i]-zList[j])
                    ));
                    A[i][j] = sqr(r) * log(r);
                }
                else
                {
                    A[i][j] = 0;
                }
            }
            // Polynomial correction
            A[i][nx] = 1;
            A[i][nx+1] = xList[i];
            A[i][nx+2] = yList[i];
            A[i][nx+3] = zList[i];
            A[nx][i] = 1;
            A[nx+1][i] = xList[i];
            A[nx+2][i] = yList[i];
            A[nx+3][i] = zList[i];
        }

        // Inverse the matrix once and store it for later iterations
        LUscalarMatrix Atemp(A);

        Info<< "Inverse RBF PHS matrix" << endl;
        Atemp.inv(invRBFmatrix_);

        // solve(w, A, vListTemp);
    }

    scalarList vListTemp(vList);
    for (int i = 0; i < 4; i++)
    {
        vListTemp.append(0);
    }

    w = invRBFmatrix_ * vListTemp;

    return(w);
}


Foam::scalar Foam::powerModels::nuclearFuelFMU::polyharmonicSpline
(
    const scalarList w,
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    const label nx(xList.size());
    scalar res(0);
    forAll(xList, i)
    {
        const scalar r(sqrt(
            sqr(x-xList[i]) 
            + sqr(y-yList[i])
            + sqr(z-zList[i])
        ));
        res += w[i] * sqr(r) * log(r);
    }
    res += w[nx] + w[nx+1]*x + w[nx+2]*y + w[nx+3]*z;
    return(res);
}


Foam::scalarList Foam::powerModels::nuclearFuelFMU::solvePolyharmonicSplineDerivative
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalarList vList,
    const scalarList dvxList,
    const scalarList dvyList
) const
{
    const label nx(xList.size());
    SquareMatrix<scalar> A(3*nx+4, 0.0);
    scalarList w(3*nx+4, 0.0);

    /*
    Info<< "Size: " << nl
        << "xList : " << xList.size() << nl
        << "yList : " << yList.size() << nl
        << "zList : " << zList.size() << nl
        << "vList : " << vList.size() << nl
        << "dvxList : " << dvxList.size() << nl
        << "dvyList : " << dvyList.size() << nl
        << endl;
    */

    forAll(xList, i)
    {
        const scalar xi(xList[i]);
        const scalar yi(yList[i]);
        const scalar zi(zList[i]);
        forAll(xList, j)
        {
            const scalar xj(xList[j]);
            const scalar yj(yList[j]);
            const scalar dx(xi-xj);
            const scalar dy(yi-yj);
            const scalar r(sqrt(
                sqr(dx) + sqr(dy) + sqr(zi-zList[j])
            ));
            const scalar rd(sqrt(0*sqr(dx) + sqr(dy)));
            const scalar phir(sqr(r) * log(r));
            if (r > 0)
            {

                // wi with fi (A)
                A[i][j] = phir;
                // vix with fi (C)
                A[i][nx+j] = dx * phir;
                // viy with fi (H)
                A[i][2*nx+j] = dy * phir;
            }
            if (rd > 0)
            {
                const scalar dphix((2.0 * dx * log(rd) + dx)*0);
                const scalar dphiy(2.0 * dy * log(rd) + dy);
                // wi with df/dx (D)
                A[nx+i][j] = dphix;
                // vix with df/dx (E)
                A[nx+i][nx+j] = phir + dx * dphix;
                // viy with df/dx (F)
                A[nx+i][2*nx+j] = dy * dphix;
                
                // wi with df/dy (I)
                A[2*nx+i][j] = dphiy;
                // vix with df/dy (K)
                A[2*nx+i][nx+j] = dx * dphiy * 0;
                // viy with df/dy (J)
                A[2*nx+i][2*nx+j] = phir + dy * dphiy;
            }
        }
        // Polynomial correction (B)
        A[i][3*nx] = 1;
        A[i][3*nx+1] = xi;
        A[i][3*nx+2] = yi;
        A[i][3*nx+3] = zi;
        A[3*nx][i] = 1;
        A[3*nx+1][i] = xi;
        A[3*nx+2][i] = yi;
        A[3*nx+3][i] = zi;

        // Polynome for df/dx (G)
        A[nx+i][3*nx+1] = 1;
        // Polynome for df/dy (L)
        A[2*nx+i][3*nx+2] = 1;
    }

    scalarList vListTemp(vList);
    forAll(dvxList, i)
    {
        vListTemp.append(dvxList[i] / 114.0869 * 0.217334 * 0);
    }
    forAll(dvyList, i)
    {
        vListTemp.append(dvyList[i] / 114.0869 * 0.217334);
    }
    for (int i = 0; i < 4; i++)
    {
        vListTemp.append(0);
    }

    solve(w, A, vListTemp);

    return(w);
}


Foam::scalar Foam::powerModels::nuclearFuelFMU::polyharmonicSplineDerivative
(
    const scalarList w,
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    const label nx(xList.size());
    scalar res(0);
    forAll(xList, i)
    {
        const scalar r(sqrt(
            sqr(x-xList[i]) 
            + sqr(y-yList[i])
            + sqr(z-zList[i])
        ));
        if (r > 0)
        {
            res += (
                w[i] + w[i+nx] * (x-xList[i]) + w[i+2*nx] * (y-yList[i])
            ) * sqr(r) * log(r);
        }
    }
    res += w[3*nx] + w[3*nx+1]*x + w[3*nx+2]*y + w[3*nx+3]*z;
    return(res > 0 ? res : 0);
}


Foam::scalarList Foam::powerModels::nuclearFuelFMU::solvePolyharmonicSplineIntegral
(
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalarList vList,
    const word region
) const
{
    const label nx(xList.size());
    SquareMatrix<scalar> A(nx+1, 0.0);
    scalarList w(nx+1, 0.0);

    const labelList& regionCells(structure_.cellLists()[region]);
    const scalarList& V(mesh_.V());

    forAll(xList, i)
    {
        forAll(xList, j)
        {
            if (i != j)
            {
                const scalar r(sqrt(
                    sqr(xList[i]-xList[j]) 
                    + sqr(yList[i]-yList[j])
                    + sqr(zList[i]-zList[j])
                ));
                A[i][j] = sqr(r) * log(r);
            }
            else
            {
                A[i][j] = 0;
            }
        }
        // Polynomial correction
        A[i][nx] = 1;
        // A[i][nx+1] = xList[i];
        // A[i][nx+2] = yList[i];
        // A[i][nx+3] = zList[i];
        // A[i][nx+1] = zList[i];
        // A[nx][i] = 1;
        // A[nx+1][i] = xList[i];
        // A[nx+2][i] = yList[i];
        // A[nx+3][i] = zList[i];
        // A[nx+1][i] = zList[i];

        scalar totalVolume(0), totalPhi(0);
        forAll(regionCells, celli)
        {
            const scalar xCell(mesh_.C().internalField()[celli].x());
            const scalar yCell(mesh_.C().internalField()[celli].y());
            const scalar zCell(mesh_.C().internalField()[celli].z());
            const scalar xPos(xCell-xList[i]);
            const scalar yPos(yCell-yList[i]);
            const scalar zPos(zCell-zList[i]);
            label cellNumber = mesh_.findCell(point(xPos, yPos, zPos));
            // totalPower += //alpha_[celli]
            //     /***/ fractionOfPowerFromNeutronics_[regioni]
            //     * structure_.powerDensityNeutronics()[celli]
            //     * V[celli];
            totalVolume += V[celli];
            const scalar rsqr(sqr(xPos)+sqr(yPos)+sqr(zPos));
            totalPhi += (rsqr * log(sqrt(rsqr))) * V[cellNumber];
        }

        // Info<< totalPhi << " " << totalVolume << endl;

        A[nx][i] = totalPhi;
        A[nx][nx] = totalVolume;
    }

    scalarList vListTemp(vList);
    vListTemp.append(300e6 / 114.08690675627066);
    // vListTemp.append(0);

    solve(w, A, vListTemp);

    return(w);
}


Foam::scalar Foam::powerModels::nuclearFuelFMU::polyharmonicSplineIntegral
(
    const scalarList w,
    const scalarList xList,
    const scalarList yList,
    const scalarList zList,
    const scalar x,
    const scalar y,
    const scalar z
) const
{
    const label nx(xList.size());
    scalar res(0);
    forAll(xList, i)
    {
        const scalar r(sqrt(
            sqr(x-xList[i]) 
            + sqr(y-yList[i])
            + sqr(z-zList[i])
        ));
        res += w[i] * sqr(r) * log(r);
    }
    res += w[nx];
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
