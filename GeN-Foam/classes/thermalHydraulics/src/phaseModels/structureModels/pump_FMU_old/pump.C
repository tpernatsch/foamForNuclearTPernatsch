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
#include "commDataLayer.H"

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
    pumpDirection_(this->get<vector>("pumpDirection")),
    momentumSourceInit_(this->get<scalar>("momentumSourceInit"))
{

    Info << "Creating pump in " << dict.dictName() << endl;

    // Communicating with the FMU
    const Time& runTime = this->db().time();
    commDataLayer& data = commDataLayer::New(runTime);
    
    data.storeObj(momentumSourceInit_,"momentumSourceCoupled",commDataLayer::causality::in);
 
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


/*
void Foam::pump::correct
( 
    volVectorField& momentumSource
)
{

    const Time& runTime = this->db().time();
    commDataLayer& data = commDataLayer::New(runTime);

     momentumSource =
        data.getObj<scalar>("momentumSource",commDataLayer::causality::in);

    vector pumpValue(pumpDirection_ * momentumSource); // pumpDirection_ tbd in .H

    forAll(cellList_, i)
    {
        const label& celli(cellList_[i]);
        momentumSource[celli] = pumpValue;
    }

    momentumSource.correctBoundaryConditions();
}
*/

void Foam::pump::correct
( 
    volVectorField& momentumSource
)
{

    const Time& runTime = this->db().time();
    commDataLayer& data = commDataLayer::New(runTime);

    const scalar momentumSourceCoupled =
        data.getObj<scalar>("momentumSourceCoupled",commDataLayer::causality::in);


    //updates the vector field by adjusting the magnitude
    vector pumpValue(momentumSourceCoupled * pumpDirection_);
    //pumpValue = momentumSource * pumpDirection_; // pumpDirection_ tbd in .H


//    momentum source is not a vector field here so I create a vector field whose magnitude gets updated

     forAll(cellList_, i)
     {
        const label& celli(cellList_[i]);
        momentumSource[celli] = pumpValue;
     }
     momentumSource.correctBoundaryConditions();

}



// ************************************************************************* //
