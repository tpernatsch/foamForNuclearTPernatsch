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

Class
    Foam::regimeMapModel::noRegimeMap

Description
    No regime map. A default regime still needs to be created as the logic of
    the solvers that use this library (e.g. FFSEulerFoam) require it

Autor
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    noRegimeMap.C

\*---------------------------------------------------------------------------*/

#ifndef noRegimeMap_H
#define noRegimeMap_H

#include "regimeMapModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace regimeMapModels
{

/*---------------------------------------------------------------------------*\
                           Class noRegimeMap Declaration
\*---------------------------------------------------------------------------*/

class noRegimeMap
:
    public regimeMapModel
{
private:

public:

    //- Runtime type information
    TypeName("none");

    //- Constructors

        noRegimeMap
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& physicsModelsDict
        );

        //- Destructor
        virtual ~noRegimeMap(){}

    //- Member Functions

        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regimeMapModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //