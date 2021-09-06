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

#include "pump.H"

//- From forward declarations
#include "structure.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pump, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pump::pump
(
    const fvMesh& mesh,
    const dictionary& dict,
    const labelList& cellList
)
:
    IOdictionary
    (
        IOobject
        (
            typeName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(mesh),
    cellList_(cellList),
    pumpValue_(this->get<vector>("momentumSource")),
    t0_(0.0),
    timeDependent_(false)
{
    Info << "Creating pump in " << dict.dictName() << endl;

    word timeProfileDictName("momentumSourceTimeProfile");
    if (this->found(timeProfileDictName))
    {
        const dictionary& timeProfileDict(dict.subDict(timeProfileDictName));
        word type
        (
            timeProfileDict.get<word>("type")
        );
        timeProfilePtr_.reset        
        (
            Function1<scalar>::New
            (
                type,
                timeProfileDict,
                type
            )
        );
        timeDependent_ = true;
        t0_ = timeProfileDict.lookupOrDefault("startTime", 0.0);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pump::correct
( 
    volVectorField& momentumSource
)
{
    vector pumpValue(pumpValue_);
    if (timeDependent_)
    {
        scalar t(mesh_.time().timeOutputValue()-t0_);
        pumpValue *= timeProfilePtr_->value(t);
    }

    forAll(cellList_, i)
    {
        const label& celli(cellList_[i]);
        momentumSource[celli] = pumpValue;
    }

    momentumSource.correctBoundaryConditions();
}

// ************************************************************************* //
