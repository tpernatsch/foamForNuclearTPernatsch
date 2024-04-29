/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2312                                                  |
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::customPimpleControl::correct()
{
    return pimpleControl::correct();
}

inline bool Foam::customPimpleControl::correct(bool interrupt)
{
    setFirstIterFlag();

    nCorrPISOInPrevPIMPLE_ = corrPISO_;

    ++corrPISO_;

    if (debug)
    {
        Info<< algorithmName_ << " correct: corrPISO = " << corrPISO_ << endl;
    }

    interrupt = corrPISOUntilConvergence_ and interrupt;

    if (corrPISO_ <= nCorrPISO_ and !interrupt)
    {
        // Implicit FMI condition
        #ifdef isCommDataLayerIncluded
        if (/*corrPISO_ > 1 &&*/ isFMIcoupledCorrection_) // Skip the first correction step
        {
            commDataLayer& data = commDataLayer::New(runTime_);
            scalar& isPressureCorrectorConverged = data.getObj<scalar>
            (
                "isPressureCorrectorConverged", 
                commDataLayer::causality::out
            );
            isPressureCorrectorConverged = 0.0;
            Info<< "Pressure corrector converged: " << isPressureCorrectorConverged
                << endl;
            // Need to execute the functionObj FMUController
            label FMUControllerLabel(-1);
            forAll(runTime_.functionObjects(), labelI)
            {
                if (runTime_.functionObjects()[labelI].type() == "FMUController")
                {
                    FMUControllerLabel = labelI;
                    break; // FMUController is unique
                }
            }
            if (FMUControllerLabel != -1)
            {
                runTime_.functionObjects()[FMUControllerLabel].execute();
            }
        }
        #endif

        return true;
    }

    corrPISO_ = 0;

    setFirstIterFlag();

    // FMI loop correction
    #ifdef isCommDataLayerIncluded
    if (isFMIcoupledCorrection_)
    {
        commDataLayer& data = commDataLayer::New(runTime_);
        scalar& isPressureCorrectorConverged = data.getObj<scalar>
        (
            "isPressureCorrectorConverged", 
            commDataLayer::causality::out
        );
        isPressureCorrectorConverged = 1.0;
        Info<< "Pressure corrector converged: " << isPressureCorrectorConverged
            << endl;
        // Need to execute the functionObj FMUController to send the info back to 
        // the main script 
        label FMUControllerLabel(-1);
        forAll(runTime_.functionObjects(), labelI)
        {
            if (runTime_.functionObjects()[labelI].type() == "FMUController")
            {
                FMUControllerLabel = labelI;
                break; // FMUController is unique
            }
        }
        if (FMUControllerLabel != -1)
        {
            runTime_.functionObjects()[FMUControllerLabel].execute();
        }
    }
    #endif

    return false;
}

// ************************************************************************* //
