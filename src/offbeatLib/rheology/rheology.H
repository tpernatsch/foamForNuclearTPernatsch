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
    Foam::rheology

Description
    Parent class that handles the evolution of rheological models.

    It handles the specific constitutive laws, e.g. elasticity, plasticity etc.
    by modifying the stress tensor field.

    The parent class is not runTime selectable.

SourceFiles
    rheology.C

\todo
    Merge with materials?
    Make sigma and epsilon members of rheology?

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef rheology_H
#define rheology_H

#include "IOdictionary.H"
#include "globalOptions.H"
#include "materials.H"

#include "primitiveFields.H"
#include "volFields.H"
#include "fvMesh.H"

#include "constitutiveLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mechanicsSubSolver;
/*---------------------------------------------------------------------------*\
                         Class rheology Declaration
\*---------------------------------------------------------------------------*/

class rheology
{
    // Private data
    

    // Private Member Functions

        //- Disallow default bitwise copy construct
        rheology(const rheology&);

        //- Disallow default bitwise assignment
        void operator=(const rheology&);

protected:
        
        //- Reference to mesh
        const fvMesh& mesh_;

        // Non const reference to the materials class
        materials& mat_;

        // Rheology dictionary
        const dictionary rheologyDict_;
        
        //- Thermal strain
        volSymmTensorField epsTh_;
        
        //- Elastic strain
        volSymmTensorField epsEl_;

        //- Hydrostatic stress
        volScalarField sigmaHyd_;

        //- Gradient of hydrostatic stress
        volVectorField gradSigmaHyd_;

        //- Deviatoric stress
        volSymmTensorField sigmaDev_;

        //- Additional strain components from various models
        volSymmTensorField additionalStrain_;
        
        //- Allow thermal expansion or not
        Switch thermalExpansion_;
        
        //- Are we solving a 2D plane stress problem?
        Switch planeStress_;
        
        //- Are we solving a modified plane strain problem?
        Switch modifiedPlaneStrain_;

        //- Solve the pressure equation (decreases oscillation in certain cases)
        Switch solvePressureEqn_;

        //- Smoothing coefficient scale factor for the pressure equation
        scalar pressureSmoothingScaleFactor_;
    
        //- List of constitutive laws
        PtrList<constitutiveLaw> lawsList_;
        
        //- Initial residual of current inner iteration
        scalar residual_;
        
public:

    //- Runtime type information
    TypeName("none");
    
    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            rheology,
            dictionary,
            (
                const fvMesh& mesh,
                materials& mat,
                const dictionary& rheologyDict
            ),
            (mesh, mat, rheologyDict)
        );


    // Constructors

        //- Construct from mesh, materials and dictionary
        rheology
        (
            const fvMesh& mesh,
            materials& mat,
            const dictionary& rheologyDict
        ); 

        //- Return a pointer to a new rheology
        static autoPtr<rheology> New
        (
            const fvMesh& mesh,
            materials& mat,
            const dictionary& solverDict
        );


    //- Destructor
    virtual ~rheology();


    //- Access functions
        
        //- Return reference to constitutive laws list
        inline const PtrList<constitutiveLaw>& lawsList() const
        {
            return lawsList_;
        }


    // Member Functions

        //- Correct/update the properties
        virtual void correct
        (
            volSymmTensorField& sigma, 
            volSymmTensorField& epsilon,
            volVectorField& D
        ) = 0;

        //- Solve pressure equation
        virtual void solvePressureEqn(volScalarField& sigmaHydExplicit);

        //- Correct/update the properties
        virtual void getAdditionalStrain()
        {
            additionalStrain_ = mat_.additionalStrains();
        };

        //- Update total fields at the end of time step
        virtual void updateTotalFields()
        {};
        
        //- Return planeStress_ switch
        inline const Switch& planeStress() const
        {
            return planeStress_;
        }
        
        //- Return modifiedPlaneStrain switch
        inline const Switch& modifiedPlaneStrain() const
        {
            return modifiedPlaneStrain_;
        }

        //- Return nextDeltaT according to rheology
        virtual scalar nextDeltaT()
        {
            return GREAT;
        }; 

        //- Return nextMacroTime according to rheology
        virtual scalar nextMacroTime()
        {
            return GREAT;
        }; 
          
        //- Return maximum time according to rheology
        virtual scalar lastMacroTime()
        {
            return GREAT;
        };

        //- Return true if converged
        virtual bool converged();

        //- Return a const reference to residual_
        inline const scalar& residual() 
        {
                return residual_;
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
