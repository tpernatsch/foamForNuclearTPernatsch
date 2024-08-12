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
    Foam::multiMaterialInterface

Description
    Base class for defining the set of faces that constitute the interfaces
    between different materials.

SourceFiles
    multiMaterialInterface.C

\*---------------------------------------------------------------------------*/

#ifndef multiMaterialInterface_H
#define multiMaterialInterface_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "faceSet.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class multiMaterialInterface Declaration
\*---------------------------------------------------------------------------*/

class multiMaterialInterface
{
    // Private data
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiMaterialInterface(const multiMaterialInterface&);

        //- Disallow default bitwise assignment
        void operator=(const multiMaterialInterface&);


protected:

    // Protected data

        //- Reference to mesh
        const fvMesh& mesh_;
        
        //- Face correction weights
        surfaceScalarField weights_;
        
        
    // Protected Member Functions

public:

    // Static data members

        //- Runtime type information
        TypeName("none");

        // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            multiMaterialInterface,
            dictionary,
            (const fvMesh& mesh, const dictionary& dict),
            (mesh, dict)
        );


    // Constructors

        //- Construct from mesh and dictionary
        multiMaterialInterface
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Selectors

        //- Select from dictionary
        static autoPtr<multiMaterialInterface> New
        (
            const fvMesh& mesh,
            const dictionary& dict 
        );


    //- Destructor
    virtual ~multiMaterialInterface();


    // Member Functions
    
        //- Return a surfaceScalarField of interface weights
        inline const surfaceScalarField& weights() const
        {
            return weights_;
        }
        
        //- Return true if this interface applies a multi-material correction
        virtual bool active() const
        {
            return false;
        }

        //- Correct the face weights
        virtual void correct()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
