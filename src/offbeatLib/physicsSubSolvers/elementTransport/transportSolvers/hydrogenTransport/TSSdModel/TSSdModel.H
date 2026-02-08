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
    Foam::TSSdModel

Description
    yield stress mother class
    To use with the vonMises rheology
SourceFiles
    yieldStress.C

\mainauthor
    Matthieu Reymond - EPFL

\contribution

\date 
    Feb. 2023


\*---------------------------------------------------------------------------*/

#ifndef TSSdModel_H
#define TSSdModel_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "zeroGradientFvPatchField.H"
#include "materials.H"
#include "primitiveFields.H"
#include "fuelMaterial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class "TSSdModel" Declaration
\*---------------------------------------------------------------------------*/

class TSSdModel
{
    // Private data    


    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        TSSdModel(const TSSdModel&);

        //- Disallow default bitwise assignment
        void operator=(const TSSdModel&);
        
protected:

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Const reference to the materials class
        const materials& mat_;

        //- Input/options dictionary
        const dictionary lawDict_;

        //- terminal solubility dissolution field
        volScalarField& TSSd_;

        //scalar TSSd0;

public:
 
    //- Runtime type information
    TypeName("TSSdModel");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        TSSdModel,
        dictionary,
        (
            const fvMesh& mesh, 
            const materials& mat,
            const dictionary& lawDict
        ),
        (mesh, mat, lawDict)
    );
    // Constructors

        //- Construct from dictionary
        TSSdModel
        (
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& lawDict
        );

    // Selectors

        //- Select from dictionary
        static autoPtr<TSSdModel> New
        (   
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& lawDict
        );

    //- Destructor
   virtual ~TSSdModel();


    //- Member Functions    
    //- Update Yield Stress Value
   virtual void correctTSSd
   (
        const labelList& addr
   ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
