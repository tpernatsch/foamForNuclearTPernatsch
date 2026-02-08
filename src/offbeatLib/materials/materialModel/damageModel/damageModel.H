/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::damageModel

Description
    yield stress mother class
    To use with the vonMises rheology

SourceFiles
    yieldStress.C

\mainauthor
    Matthieu Reymond - CEA/EPFL
    Alessandro Scoalro - EPFL

\contribution

\date 
    Feb. 2023


\*---------------------------------------------------------------------------*/

#ifndef damageModel_H
#define damageModel_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class "damageModel" Declaration
\*---------------------------------------------------------------------------*/

class damageModel
{
    // Private data    


    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        damageModel(const damageModel&);

        //- Disallow default bitwise assignment
        void operator=(const damageModel&);
        
protected:
        
        //- Reference to mesh
        const fvMesh& mesh_;

        // Input/options dictionary
        const dictionary damageDict_;

public:
 
    //- Runtime type information
    TypeName("none");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        damageModel,
        dictionary,
        (
            const fvMesh& mesh, 
            const dictionary& damageDict,
            const dictionary& baseDict
        ),
        (mesh, damageDict, baseDict)
    );

    // Constructors

        //- Construct from dictionary
        damageModel
        (
            const fvMesh& mesh,
            const dictionary& damageDict,
            const dictionary& baseDict
        );

    // Selectors

        //- Select from dictionary
        static autoPtr<damageModel> New
        (   
            const fvMesh& mesh,
            const dictionary& matDict
        );

    //- Destructor
   virtual ~damageModel();


    // Member Functions    

        //- Update damage and apply to stress
        virtual void correct
        (
            volScalarField& sigmaHyd,
            volSymmTensorField& sigmaDev,
            volSymmTensorField& epsilonEl,
            const labelList& addr
        ){};

        //- Soften E and nu
        virtual void soften
        (
            scalarField& E,
            scalarField& nu,
            const labelList& addr
        ){};


    // Access Functions

        //- Return a constant reference to the damage field
        virtual const volScalarField& d() const
        {
            return volScalarField::null();
        };


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
