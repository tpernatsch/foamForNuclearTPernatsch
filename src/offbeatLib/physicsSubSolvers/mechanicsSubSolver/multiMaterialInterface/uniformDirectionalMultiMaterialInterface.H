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
    Foam::uniformDirectionalMultiMaterialInterface

Description
    Interfaces between materials are defined by a directionally-
    dependent weight which is optionally defined in rotated
    coordinates.

SourceFiles
    uniformDirectionalMultiMaterialInterface.C

\*---------------------------------------------------------------------------*/

#ifndef uniformDirectionalMultiMaterialInterface_H
#define uniformDirectionalMultiMaterialInterface_H

#include "uniformMultiMaterialInterface.H"
#include "coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class uniformDirectionalMultiMaterialInterface Declaration
\*---------------------------------------------------------------------------*/

class uniformDirectionalMultiMaterialInterface
:
    public uniformMultiMaterialInterface
{
    // Private data

        //- True if weights are rotated
        bool rotate_;
        
        //- Rotational tensor defining principle directions
        autoPtr<coordinateRotation> rotation_;
        
        //- Default weight to apply uniformly for the
        //  additional component of the strain
        vector defaultWeightEps_;
    
        //- Default weight to apply uniformly for the
        //  explicit components of the twoMuLambda*gradD term
        vector defaultWeightGrad_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        uniformDirectionalMultiMaterialInterface(const uniformDirectionalMultiMaterialInterface&);

        //- Disallow default bitwise assignment
        void operator=(const uniformDirectionalMultiMaterialInterface&);


public:

    //- Runtime type information
    TypeName("uniformDirectional");


    // Constructors

        //- Construct from Dictionary
        uniformDirectionalMultiMaterialInterface
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    ~uniformDirectionalMultiMaterialInterface();


    // Member Functions

        //- Correct the face weights
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
