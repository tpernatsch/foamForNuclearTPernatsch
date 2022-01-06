/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "FFPair.H"
#include "latentHeatModel.H"
#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::latentHeatModel> Foam::latentHeatModel::New
(
    const phaseChangeModel& pcm,
    const dictionary& dict, 
    const objectRegistry& objReg
)
{
    word type(dict.lookup("type"));

    Info<< "Selecting latentHeatModel: "
        << type << endl;

    auto* ctorPtr = latentHeatModelsConstructorTable(type);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown latentHeatModel type "
            << type << endl << endl
            << "Valid latentHeatModel types are : " << endl
            << latentHeatModelsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return
        autoPtr<latentHeatModel>
        (
            ctorPtr(pcm, dict, objReg)
        );
}


// ************************************************************************* //
