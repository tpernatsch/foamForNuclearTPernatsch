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
    Foam::yieldStressModel

Description
	The `yieldStressModel` class is the mother class for all the types of yield
	stresses implemented in OFFBEAT.
	It is intended to be used with the 'misesPlasticity' constitutive mechanical
	behaviour law.

SourceFiles
    yieldStress.C

\mainauthor
    Matthieu Reymond - EPFL

\contribution

\date 
    Feb. 2023


\*---------------------------------------------------------------------------*/

#ifndef yieldStressModel_H
#define yieldStressModel_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class "yieldStressModel" Declaration
\*---------------------------------------------------------------------------*/

class yieldStressModel
{
    // Private data    


    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        yieldStressModel(const yieldStressModel&);

        //- Disallow default bitwise assignment
        void operator=(const yieldStressModel&);
        
protected:

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Input/options dictionary
        const dictionary lawDict_;
        
        //- yield stress field
        volScalarField& sigmaY_;
        
        //- Hardening modulus
        volScalarField& Hp_;

public:
 
    //- Runtime type information
    TypeName("YieldStressModel");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        yieldStressModel,
        dictionary,
        (
            const fvMesh& mesh, 
            const dictionary& lawDict
        ),
        (mesh, lawDict)
    );
    // Constructors

        //- Construct from dictionary
        yieldStressModel
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );

    // Selectors

        //- Select from dictionary
        static autoPtr<yieldStressModel> New
        (   
            const fvMesh& mesh,
            const dictionary& lawDict
        );

    //- Destructor
   virtual ~yieldStressModel();


    //- Member Functions    
    //- Update Yield Stress Value
   virtual void correctYieldStress
   (
        const labelList& addr
   ) = 0;

   // Access functions

   //- Return yield stress
   const volScalarField& sigmaY() const
   {
        return sigmaY_;
   }
    //- Return yield stress, non-const
   volScalarField& sigmaY()
   {
        return sigmaY_;
   }

   //- Return Hp
   virtual const volScalarField& Hp() const
   {
        return Hp_;
   }   
   
   //- Return Hp
   virtual volScalarField& Hp()
   {
        return Hp_;
   }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
