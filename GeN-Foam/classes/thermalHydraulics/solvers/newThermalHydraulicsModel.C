/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "thermalHydraulicsModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::thermalHydraulicsModel> Foam::thermalHydraulicsModel::New
(
    Time& time,
    fvMesh& mesh,
    customPimpleControl& pimple,
    fv::options& fvOpt
)
{
    const dictionary& dict
    (
        IOdictionary
        (
            IOobject
            (
                "phaseProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    word type(dict.lookup("thermalHydraulicsType"));

    Info<< "Selecting thermalHydraulicsModel of type: " << type << endl;

    auto* ctorPtr = thermalHydraulicsModelsConstructorTable(type);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown thermalHydraulicsModel of type: "
            << type << endl << endl
            << "Valid thermalHydraulicsModel types are: " << endl
            << thermalHydraulicsModelsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return
        autoPtr<thermalHydraulicsModel>
        (
            ctorPtr(time, mesh, pimple, fvOpt)
        );

}


// ************************************************************************* //
