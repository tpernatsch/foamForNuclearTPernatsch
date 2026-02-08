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
    Foam::HDiffusionCoefficientModel

Description
    Base class for the diffusion coefficient of hydrogen ions
    
SourceFiles
    HDiffusionCoefficientModel.C

\mainauthor
    Matthieu Reymond - EPFL

\contribution

\date 
    Feb. 2023


\*---------------------------------------------------------------------------*/

#ifndef HDiffusionCoefficientModel_H
#define HDiffusionCoefficientModel_H

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
                        Class "HDiffusionCoefficientModel" Declaration
\*---------------------------------------------------------------------------*/

class HDiffusionCoefficientModel
{
    // Private data    


    // Private Member Functions
        //- Disallow default bitwise copy construct
        HDiffusionCoefficientModel(const HDiffusionCoefficientModel&);

        //- Disallow default bitwise assignment
        void operator=(const HDiffusionCoefficientModel&);
        
protected:

        //- Reference to mesh
        const fvMesh& mesh_;
    
        //- Const reference to the materials class
        const materials& mat_;

        //- Input/options dictionary
        const dictionary lawDict_;

        //- diffusion coefficient field
        volScalarField& diffH_;

        bool effectiveDiffCoeff_;

        virtual void correctEffectiveDiffusionCoefficient(const labelList& addr);

public:
 
    //- Runtime type information
    TypeName("HDiffusionCoefficientModel");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        HDiffusionCoefficientModel,
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
        HDiffusionCoefficientModel
        (
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& lawDict
        );

    // Selectors

        //- Select from dictionary
        static autoPtr<HDiffusionCoefficientModel> New
        (   
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& lawDict
        );

    //- Destructor
   virtual ~HDiffusionCoefficientModel();


    //- Member Functions    
    //- Update diffusion coefficient Value
   virtual void correctDH
   (
        const labelList& addr
   ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
