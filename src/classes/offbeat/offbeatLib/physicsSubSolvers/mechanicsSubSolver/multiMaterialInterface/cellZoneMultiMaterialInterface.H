/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::cellZoneMultiMaterialInterface

Description
    Materials are defined by cellZones, and interfaces between materials lie
    at the interface between cellZones.

SourceFiles
    cellZoneMultiMaterialInterface.C

\*---------------------------------------------------------------------------*/

#ifndef cellZoneMultiMaterialInterface_H
#define cellZoneMultiMaterialInterface_H

#include "uniformMultiMaterialInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class cellZoneMultiMaterialInterface Declaration
\*---------------------------------------------------------------------------*/

class cellZoneMultiMaterialInterface
:
    public uniformMultiMaterialInterface
{
    // Private data
    
        //- Weights to apply at the interface between cellZones
        scalar interfaceWeight_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellZoneMultiMaterialInterface(const cellZoneMultiMaterialInterface&);

        //- Disallow default bitwise assignment
        void operator=(const cellZoneMultiMaterialInterface&);


public:

    //- Runtime type information
    TypeName("cellZone");


    // Constructors

        //- Construct from dictionary
        cellZoneMultiMaterialInterface
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    ~cellZoneMultiMaterialInterface();


    // Member Functions

        //- Correct the face weights
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
