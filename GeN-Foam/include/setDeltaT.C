/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "setDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::setDeltaT(Time& runTime,  PtrList<solver>& solvers)
{
    runTime.setDeltaT
    (
        min
        (
            runTime.deltaTValue(), 
            runTime.controlDict().get<scalar>("deltaT")
        )
    );
}


void Foam::adjustDeltaT(Time& runTime,  PtrList<solver>& solvers)
{
    // Update the time-step limited by the solvers maxDeltaT
    if (runTime.controlDict().lookupOrDefault("adjustTimeStep", false))
    {
        bool transient = false;
        scalar deltaT = VGREAT;

        forAll(solvers, i)
        {
            if (solvers[i].transient())
            {
                transient = true;
                deltaT = min(deltaT, solvers[i].maxDeltaT());
                // Info <<"New delta T is "<< solvers[i].maxDeltaT()<< endl;
            }
        }

        if (transient && deltaT < ROOTVGREAT)
        {
            runTime.setDeltaT
            (
               min
                (
                   solver::deltaTFactor*runTime.deltaTValue(),
                    deltaT
                )
            );
            Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
        }
    }
}


// ************************************************************************* //
