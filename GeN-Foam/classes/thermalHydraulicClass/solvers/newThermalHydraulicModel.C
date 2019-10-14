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

#include "thermalHydraulicModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::thermalHydraulicModel> Foam::thermalHydraulicModel::New
(
    const fvMesh& mesh,
    pimpleControl& pimple,
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

    word type(dict.lookup("thermalHydraulicType"));

    Info<< "Selecting thermalHydraulicModel of type: " << type << endl;

    thermalHydraulicModelsConstructorTable::iterator cstrIter =
        thermalHydraulicModelsConstructorTablePtr_->find(type);

    if (cstrIter == thermalHydraulicModelsConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown thermalHydraulicModel of type: "
            << type << endl << endl
            << "Valid thermalHydraulicModel types are: " << endl
            << thermalHydraulicModelsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(mesh, pimple, fvOpt);
}


// ************************************************************************* //
