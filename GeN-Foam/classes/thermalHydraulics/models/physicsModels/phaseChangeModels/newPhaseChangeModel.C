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

#include "phaseChangeModel.H"
#include "fluid.H"
#include "saturationModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseChangeModel> Foam::phaseChangeModel::New
(
    const dictionary& dict,
    const pimpleControl& pimple, 
    const fluid& fluid1,
    const fluid& fluid2,
    const volScalarField& p,
    const volScalarFieldTable& htcs,
    volScalarField& dmdt,
    volScalarField& iT,
    volScalarField& iA
)
{
    word type(dict.lookup("type"));

    Info<< "Selecting phaseChangeModel: "
        << type << endl;

    phaseChangeModelsConstructorTable::iterator cstrIter =
        phaseChangeModelsConstructorTablePtr_->find(type);

    if (cstrIter == phaseChangeModelsConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown phaseChangeModel type "
            << type << endl << endl
            << "Valid phaseChangeModel types are : " << endl
            << phaseChangeModelsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(dict, pimple, fluid1, fluid2, p, htcs, dmdt, iT, iA);
}


// ************************************************************************* //
