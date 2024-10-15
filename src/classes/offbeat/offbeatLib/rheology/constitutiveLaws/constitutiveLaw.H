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
    Foam::constitutiveLaw

Description
    General class that handles the evolution of the constituveLaw (mechanical behaviour)

    It handles the specific constitutive laws, e.g. elasticity, plasticity etc.
    by modifying the stress tensor field.

    The parent class is not runTime selectable.

    \htmlinclude constituveLaw_constituveLaw.html

SourceFiles
    constitutiveLaw.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef constitutiveLaw_H
#define constitutiveLaw_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "zeroGradientFvPatchField.H"
#include "damageModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermomechanicsSubSolver Declaration
\*---------------------------------------------------------------------------*/

class constitutiveLaw
{
    // Private data    


    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        constitutiveLaw(const constitutiveLaw&);

        //- Disallow default bitwise assignment
        void operator=(const constitutiveLaw&);
        
protected:
        
        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to the law dictionary
        const dictionary lawDict_;

        //- autoPtr to damageModel (none by default)
        autoPtr<damageModel> damageModel_;

public:
 
    //- Runtime type information
    TypeName("constitutiveLaw");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            constitutiveLaw,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& lawDict
            ),
            (mesh, lawDict)
        );


    // Constructors

        //- Construct from mesh, materials , dict
        constitutiveLaw
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );


    // Selectors

        //- Selector
        static autoPtr<constitutiveLaw> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
   virtual ~constitutiveLaw();


    // Member Functions
    
        //- Update hydrostatic and deviatoric component of stress
        virtual void correct
        (
            volScalarField& sigmaHyd, 
            volSymmTensorField& sigmaDev, 
            volSymmTensorField& epsEl,
            const labelList& addr
        ) = 0;

        //- Update stress values by summing hydrostatic and deviatoric component   
        virtual void updateStress
        (
            volSymmTensorField& sigma, 
            volScalarField& sigmaHyd, 
            volSymmTensorField& sigmaDev,
            const labelList& addr
        );

        //- Correct the elastic strain tesnor removing oldTime copmonents of 
        //  plasticstrains  
        virtual void correctEpsilonEl
        (
            volSymmTensorField& epsilonEl,
            const labelList& addr
        ){};

        //- Correct the additionalStrain field given through argument
        virtual void correctAdditionalStrain
        (
            volSymmTensorField& additionalStrain,
            const labelList& addr
        ){};

        //- Update total fields (if the law is incremental e.g. plasticity)
        virtual void updateTotalFields
        (
            const labelList& addr
        ) {};

        //- correctAxialStrain in case of plane stress
        virtual void correctAxialStrain
        (
            volSymmTensorField& epsilon,
            const labelList& addr
        ) {};

        //- Return nextDeltaT according to rheology law
        virtual scalar nextDeltaT(const labelList& addr)
        {
            return GREAT;
        }; 

        //- Return nextMacroTime according to rheology law
        virtual scalar nextMacroTime()
        {
            return GREAT;
        }; 
          
        //- Return lastMacroTime according to rheology law
        virtual scalar lastMacroTime()
        {
            return GREAT;
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
