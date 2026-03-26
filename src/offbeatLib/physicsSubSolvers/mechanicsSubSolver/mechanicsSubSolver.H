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
    Foam::mechanicsSubSolver

Description
    When using the `constant` mechanics solver class, the displacement field `D`
    remains unchanged from its initial configuration throughout the simulation.

    \note
    While the `constant` solver keeps the `D` field fixed by default, it is possible
    to alter the `D` field during the simulation by leveraging advanced OpenFOAM
    features like `#codeStream#`, `funkySetFields`, or specific function objects.

    \warning
    For backward-compatibility, `constant` can be still used as the typename
    for the `constant` mechanics solver. However, its use is deprecated.

________________________________________________________________________________

\par Options
	As a `mechanicsSubSolver` class, `constant` **allows** the user to specify
	a few additional parameters in the `mechanicsOptions` sub-dictionary located
	within the `solverDict`.

    Parameters in `mechanicsOptions`:

	- <b>`forceSummary`</b> - It prints a summary of boundary-forces at each iteration.
    **It is	set to false by default**.
	- <b>`cylindricalStress`</b> - Creates cylindrical stress, strain and displacement field.
    It is useful for visualization and data-analysis purposes. **It is set to false by
    default**.
	- <b>`sphericalStress`</b> - Creates spherical stress, strain and displacement field.
    It is useful for visualization and data-analysis purposes. **It is set to false by
    default**.

________________________________________________________________________________

Usage
	Here is a code snippet of the `solverDict` to be used for activating the
	`constant` mechanics solver class:

    \verbatim
        mechanicsSolver constant;

        mechanicsOptions
        {
            // Print summary of boundary-forces at each iteration
            forceSummary false;

            // Creates cylindrical stress and displacement field (for debugging)
            cylindricalStress false;
        }
    \endverbatim

SourceFiles
    mechanicsSubSolver.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef mechanicsSubSolver_H
#define mechanicsSubSolver_H

#include "fvMesh.H"
#include "IOdictionary.H"

#include "regionCoupledOFFBEATFvPatch.H"
#include "implicitGapContactFvPatchVectorField.H"
#ifdef OPENFOAMESI
#include "volPointInterpolation.H"
#include "processorPolyPatch.H"
#endif

#include "materials.H"
#include "rheology.H"
#include "pointIOField.H"
#include "multiMaterialInterface.H"
#include "fvMatrices.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class mechanicsSubSolver Declaration
\*---------------------------------------------------------------------------*/

class mechanicsSubSolver
{
    // Private data


    // Private Member Functions

        //- Disallow default bitwise copy construct
        mechanicsSubSolver(const mechanicsSubSolver&);

        //- Disallow default bitwise assignment
        void operator=(const mechanicsSubSolver&);

protected:

        //- Relax the displacement equation with the supplied relaxation
        //- factor, using the old- and old-old field values
        static void relax
        (
            fvVectorMatrix& fvm,
            const vectorField& psi0,
            const vectorField& psi00,
            const scalar alpha
        );

        //- Relax the displacement equation with the fvSolution relaxation
        //- factor, using the old- and old-old field values
        static void relax
        (
            fvVectorMatrix& fvm,
            const vectorField& psi0,
            const vectorField& psi00
        );

        //- Reference to mesh
        fvMesh& mesh_;

        // Non const reference to the materials class
        const materials& mat_;

        //- Material constitutive laws
        rheology& rheo_;

        //- Input/options dictionary
        const dictionary mechanicsDict_;

        //- Displacement field
        volVectorField D_;

        //- Gradient of displacement field
        volTensorField gradD_;

        //- Stress tensor
        //- TODO: it should be moved to rheology
        volSymmTensorField sigma_;

        //- Equivalent or Von Mises stress
        volScalarField sigmaEq_;

        //- Explicit component of the stress tensor
        //- TODO: it should be moved to rheology
        volTensorField sigmaExp_;

        //- Divergence of the explicit component of the stress tensor
        //- TODO: it should be moved to rheology
        volVectorField divSigmaExp_;

        //- Strain Tensor
        //- TODO: it should be moved to rheology
        volSymmTensorField epsilon_;

        //- Interface pressure (for contact cases)
        volScalarField interfaceP_;

        //- Gap width (for contact cases)
        volScalarField gapWidth_;

        //- Apply momentum stabilistation correction
        Switch RhieChowCorrection_;

        //- Final residual of latest inner iteration
        vector residual_;

        //- Final residual of firs inner loop iteration
        vector initialResidual_;

        //- Different residual based on the relative change of D
        scalar relResidual_;          
        
        //- Final relative residual of firs inner loop iteration
        scalar initialRelResidual_;  
        
        //- Absolute relative errod for D or DD fields
        scalar absErr_;    
        
        //- Final abs err of firs inner loop iteration
        scalar initialAbsErr_;  
        
        //- Reference cell to fix if model is not fixed
        label refCell_;

        //- Enable multi-material corrections
        Switch multiMaterialCorrection_;

        //- Selection of faces for explicit multi-material corrections
        autoPtr<multiMaterialInterface> materialInterface_;

        //- Flag to determine if force summary should be printed
        Switch forceSummary_;

        //- Tranform to and write stresses, strains and displacements in
        //- cylindrical coordinates
        Switch cylindricalStress_;

        //- Tranform to and write stresses, strains and displacements in
        //- spherical coordinates
        Switch sphericalStress_;

        //- Stress tensor in cylindrical coordinates
        //- TODO: it should be moved to rheology
        autoPtr<volSymmTensorField> sigmaCyl_;

        //- Displacement in cylindrical coordinates
        autoPtr<volSymmTensorField> epsilonCyl_;

        //- Displacement in cylindrical coordinates
        autoPtr<volVectorField> DCyl_;

        //- Stress tensor in spherical coordinates
        //- TODO: it should be moved to rheology
        autoPtr<volSymmTensorField> sigmaSph_;

        //- Displacement in spherical coordinates
        autoPtr<volSymmTensorField> epsilonSph_;

        //- Displacement in spherical coordinates
        autoPtr<volVectorField> DSph_;

     // Protected member functions

        //- Find nearest cell to volumetric centre of mesh
        void findRefCell();

        //- Read and return number of inner correctors
        int nCorrectors(const word);
        
        //- Return true if a reference cell is set
        inline bool hasReference() const
        {
            return refCell_ != -1;
        }

        //- Find volumetric average of a field
        template<class T>
        inline T domainAverage(const Field<T>& psi)
        {
            const scalarField& V = mesh_.V().field();
            return gSum(psi*V)/gSum(V);
        }

        //- Disable solution in the out-of-plane direction for wedge cases
        virtual void checkWedges() const;

        //- Calculate gapWidth
        virtual tmp<volScalarField> calculateInitialGapWidth() const;

        //- Print surface force summary
        void forceSummary();

        //- Calculate stress and displacement in cylindrical coordinates
        void correctCylindrical();

        //- Calculate stress and displacement in spherical coordinates
        void correctSpherical();

        //- Return $2*\mu + \lambda$ at the faces including corrections
        tmp<surfaceScalarField> correctedTwoMuLambdaF
        (
            const volScalarField& twoMuLambda
        ) const;

        //- Return the gradient of displacement grad(D) including corrections
        tmp<volTensorField> correctedGradD
        (
            const volScalarField& twoMuLambda,
            const volVectorField& D,
            const volTensorField& gradD,
            const volTensorField& sigmaExp
        ) const;

        //- Return the gradient of incremental displacement grad(DD) including corrections
        tmp<volTensorField> correctedGradDD
        (
            const volScalarField& twoMuLambda,
            const volVectorField& DD
        ) const;

        //- Return the fvm laplacian of the displacement (i.e. the divergence
        //  of the implicit component of the traction)
        tmp<fvMatrix<vector>> fvmLaplacianD
        (
            const surfaceScalarField& twoMuLambda_f,
            const volVectorField& D,
            const word& laplacianName
        ) const;  

        //- Return the fvm laplacian of the displacement (i.e. the divergence
        //  of the implicit component of the traction)
        tmp<fvMatrix<vector>> fvmLaplacianD
        (
            const surfaceTensorField& twoMuLambda_f,
            const volVectorField& D,
            const word& laplacianName
        ) const;    

        //- Return the corrected explicit component of the traction
        template<class Type>
        tmp<volTensorField> correctedSigmaExp
        (
            const GeometricField<Type, fvPatchField, volMesh>& sigma,
            const volScalarField& twoMuLambda,
            const volTensorField& gradD
        ) const;

        //- Return the divergence of the explicit component of the traction
        template<class Type>
        tmp<volVectorField> correctedDivSigmaExp
        (
            const volTensorField& sigmaExp,
            const GeometricField<Type, fvPatchField, volMesh>& sigma,
            const volScalarField& twoMuLambda,
            const surfaceScalarField& twoMuLambda_f,
            const volVectorField& D,
            const volTensorField& gradD
        ) const;

public:

    //- Runtime type information
    TypeName("constant");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            mechanicsSubSolver,
            dictionary,
            (
                fvMesh& mesh,
                const materials& materials,
                rheology& rheo,
                const dictionary& mechanicsDict
            ),
            (mesh, materials, rheo, mechanicsDict)
        );

    // Constructors

        //- Construct from mesh, materials, rheology and dictionary
        mechanicsSubSolver
        (
            fvMesh& mesh,
            const materials& materials,
            rheology& rheo,
            const dictionary& mechanicsDict
        );

        //- Return a pointer to a new thermoMechanicalProperties
        static autoPtr<mechanicsSubSolver> New
        (
            fvMesh& mesh,
            const materials& materials,
            rheology& rheo,
            const dictionary& solverDict
        );

    //- Destructor
        virtual ~mechanicsSubSolver();


    // Member Functions

        //- Correct/update the properties
        virtual void correct()
        {}

        //- Update total fields at the end of time step
        virtual void updateTotalFields();

        //- Update incremental fields at the end of time step
        virtual void updateIncrementalFields()
        {}

        //- Return true if converged
        virtual bool converged
        (
            volVectorField& D, 
            const vector& residual,
            const scalar& relResidual,
            const scalar& absErr
        );

        //- Return true if converged
        virtual bool converged(volVectorField& D);

        //- Return true if converged
        virtual bool converged()
        {
            return converged(D_);
        };

        //- Update mesh
        virtual void updateMesh(){};

    //- Access functions

        //- Return mesh
        inline const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return stress tensor field
        inline volSymmTensorField& sigma()
        {
            return sigma_;
        }

        //- Return stress tensor field
        inline const volSymmTensorField& sigma() const
        {
            return sigma_;
        }

        //- Return non-const displacement field
        inline volVectorField& D()
        {
            return D_;
        }

        //- Return displacement field
        inline const volVectorField& D() const
        {
            return D_;
        }

        //- Return total or incremental displacement depending on
        //- which one should be used to move mesh points (useful when
        //- not knowing if the instance is an updateMesh solver or not
        virtual inline volVectorField& DorDD()
        {
            return D_;
        }

        //- Return total or incremental displacement depending on
        //- which one should be used to move mesh points (useful when
        //- not knowing if the instance is an updateMesh solver or not
        virtual inline const volVectorField& DorDD() const
        {
            return D_;
        }

        //- Return non-const gradient field
        inline volTensorField& gradD()
        {
            return gradD_;
        }

        //- Return gradient field
        inline const volTensorField& gradD() const
        {
            return gradD_;
        }

        //- Return gap width field
        inline const volScalarField& gapWidth() const
        {
            return gapWidth_;
        }

        //- Return strain tensor field
        inline volSymmTensorField& epsilon()
        {
            return epsilon_;
        }

        //- Return strain tensor field
        inline const volSymmTensorField& epsilon() const
        {
            return epsilon_;
        }

        //- Return initialResidual_
        inline vector residual() const
        {
                return residual_;
        };

        //- Return firstResidual_
        inline vector initialResidual() const
        {
                return initialResidual_;
        };

        //- Return relativeRes_
        inline scalar relResidual() const
        {
                return relResidual_;
        };

        //- Return a reference to the powerModule
        inline rheology& rheo()
        {
            return rheo_;
        }

        //- Return a reference to the materials properties
        inline const materials& mat() const
        {
            return mat_;
        }

        //- Return a reference to the mechanicsSolver dictionary
        inline const dictionary& mechanicsDict() const
        {
            return mechanicsDict_;
        }

        //- Relax main field
        virtual void relax()
        {
            D_.relax();
        };

        virtual scalar nextDeltaT()
        {
            return rheo_.nextDeltaT();
        };

        virtual scalar nextMacroTime()
        {
            return GREAT;
        };

        //- Return maximum time
        virtual scalar lastMacroTime()
        {
            return GREAT;
        };

        //- Couple or uncouple implicit contact boundaries
        virtual void coupleBoundaries
        (
            volVectorField& D,
            bool coupleOrUnCouple
        );

     // Public static member functions

        //- Update mesh
        static void updateMesh
        (
            fvMesh& mesh,
            const volVectorField& movingD,
            const bool& keepFlatEmptyFaces = true
        );

        //- Disable solution in the out-of-plane direction for wedge cases
        static void checkWedges
        (
            fvMesh& mesh
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mechanicsSubSolverTemplates.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
