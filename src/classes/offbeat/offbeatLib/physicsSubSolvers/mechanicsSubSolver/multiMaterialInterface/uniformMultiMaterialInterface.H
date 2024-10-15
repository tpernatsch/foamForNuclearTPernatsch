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
    Foam::uniformMultiMaterialInterface

Description
    Interfaces between materials are defined by a uniform.

SourceFiles
    uniformMultiMaterialInterface.C

\*---------------------------------------------------------------------------*/

#ifndef uniformMultiMaterialInterface_H
#define uniformMultiMaterialInterface_H

#include "multiMaterialInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class uniformMultiMaterialInterface Declaration
\*---------------------------------------------------------------------------*/

class uniformMultiMaterialInterface
:
    public multiMaterialInterface
{
    // Private data
    
        //- Default weight to apply uniformly
        scalar defaultWeight_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        uniformMultiMaterialInterface(const uniformMultiMaterialInterface&);

        //- Disallow default bitwise assignment
        void operator=(const uniformMultiMaterialInterface&);


public:

    //- Runtime type information
    TypeName("uniform");


    // Constructors

        //- Construct from Dictionary
        uniformMultiMaterialInterface
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    ~uniformMultiMaterialInterface();


    // Member Functions

        //- Return true if this interface applies a multi-material correction
        virtual bool active() const
        {
            return true;
        }
        
        //- Correct the face weights
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
