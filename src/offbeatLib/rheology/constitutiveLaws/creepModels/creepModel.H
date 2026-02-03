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
    Foam::creepModel

Description
    Creep model mother class.

    \htmlinclude creepModel_creepModel.html

SourceFiles
    creepModel.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef creepModel_H
#define creepModel_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "zeroGradientFvPatchField.H"
#include "userParameters.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class "creepModel" Declaration
\*---------------------------------------------------------------------------*/

class creepModel
{
    // Private data    


    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        creepModel(const creepModel&);

        //- Disallow default bitwise assignment
        void operator=(const creepModel&);
        
protected:
        
        //- Reference to mesh
        const fvMesh& mesh_;

        // Input/options dictionary
        const dictionary lawDict_;

        //- Creep Strain 
        volSymmTensorField& epsilonCreep_;

        //- Increment of creep Strain 
        volSymmTensorField& DepsilonCreep_;
        
        //- Equivalent Creep Strain 
        volScalarField& epsilonCreepEq_;
        
        //- Increment of equivalent Creep Strain 
        volScalarField& DepsilonCreepEq_;

public:
 
    //- Runtime type information
    TypeName("fromLatestTime");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            creepModel,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& lawDict
            ),
            (mesh, lawDict)
        );


    // Constructors

        //- Construct from dictionary
        creepModel
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );


    // Selectors

        //- Select from dictionary
        static autoPtr<creepModel> New
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );


    //- Destructor
   virtual ~creepModel();

    // Access functions

        virtual const volScalarField& DepsilonCreepEq() const
        {
            return DepsilonCreepEq_;
        }

    // Member Functions    

        virtual void correctCreep
        (
            volSymmTensorField& epsilonEl,
            const labelList& addr
        );  

        virtual const volSymmTensorField& additionalStrain() 
        {
            return epsilonCreep_;
        }; 

        //- correctAxialStrain in case of plane stress
        virtual void correctAxialStrain
        (
            volSymmTensorField& epsilon,
            const labelList& addr
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
