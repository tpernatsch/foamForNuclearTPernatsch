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
    Foam::largeStrainUpdLag

Description
    This class olves for incremental displacement `DD` in an updated Lagrangian configuration 
    (i.e. the mesh is updated at every time step), and employs a finite- (or large-) 
    strain mechanics framework.
    
    The geometric non-linearity (i.e. the inter-dependence between displacement
    and geometric configuration) is taken care of by multiplying the stress 
    tensor (Cauchy) by relJ*relF.inv(), where relF is the relative deformation 
    gradient and relJ is its determinant.

    By default, the strain tensor takes into account the second order deformations:

    $$
	\begin{aligned}
	E = \frac{1}{2} (\nabla D + \nabla^T D + \nabla D \cdot \nabla^T D)
	\end{aligned}
	$$
     
    However, with the "linearStrainTensor" keyword, the user can decide to force
    a linear strain tensor:

	$$
	\begin{aligned}
	\epsilon = \frac{1}{2} (\nabla D + \nabla^T D)
	\end{aligned}
	$$
	
________________________________________________________________________________

\par Options	
	As a `mechanicsSubSolver` class, `largeStrainUpdLag` **allows** the user to 
    specify a few additional parameters in the `mechanicsOptions` sub-dictionary located 
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
    - <b>`strainTensor`</b> - The strain tensor can be selected as `linear`, `EulerAlmansi`
    **(default choice)**, `Hencky`, `GreenLagrange`.
    - <b>`RhieChowCorrection`</b> - Apply a correction to the explicit component of the 
    divergence of the stress, similar  to a Rhie & Chow correction, to stabilize the momentum equation. 
    **It is set to true by default unless the `multiMaterialCorrection` is activated.** 
    It is strongly suggested to use a RhieChowCorrection unless the `multiMaterialCorrection` 
    itself is activated.  
    - <b>`RhieChowScaleFactor`</b> - Scaling factor for the RhieChowCorrection. **A value of 1 
    is used by default** and corresponds to applying the full correction. A value of 0 applies 
    no correction at all.  
    - <b>`multiMaterialCorrection`</b> - Additional subdictionary that applies a custom discretization 
    scheme for both the implicit and explicit components of the stress divergence. 
    This means that the gradient scheme for the stress and the laplacian scheme 
    for the displacement fields are neglected. **The correction is not activated by default.**
    - <b>`inertialRelaxation`</b> - Optionally apply inertial relaxation to the
    solution of the displacement equation. This may be activated in models containing
    regions that are not fixed, e.g. free pellets, and are therefore numerically
    ill-posed.
    - <b>`inertialRelaxationDeltaT`</b> - If inertial relaxation is activated, a
    virtual time step is defined to reduce the rate of solution change with each
    iteration. A time step of 1e-7 has been found to work well in most cases.

________________________________________________________________________________

Usage
	Here is a code snippet of the `solverDict` to be used for activating the
	`largeStrainUpdLag` mechanics solver class:
    
    \verbatim
        mechanicsSolver largeStrainUpdLag;

        mechanicsOptions
        {
            // Print summary of boundary-forces at each iteration
            forceSummary false;
            
            // Creates cylindrical stress and displacement field
            cylindricalStress false;
            
            // Creates spherical stress and displacement field
            sphericalStress false;

            // Type of strain tensor
            strainTensor     Hencky;
            
            // Apply RhieChowCorrection to stabilize momentum equation (it should
            // be activated in most situations)
            RhieChowCorrection true;
            
            // Optional subDict for multi-material cases (it has priority over
            // RhieChowCorrection when activated)
            multiMaterialCorrection
            {
                // other correction types are possible
                type none;
            }

            // Activate inertial relaxation (optional, default=false)
            inertialRelaxation false;

            // Inertial relaxation time step (optional, default=1e-7)
            inertialRelaxationDeltaT 1e-7;
        }
    \endverbatim
    
SourceFiles
    largeStrainUpdLag.C

\mainauthor
    E. Brunetto, A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    C. Fiorina - EPFL\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef largeStrainUpdLagUpdLag_H
#define largeStrainUpdLag_H

#include "mechanicsSubSolver.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                Class largeStrainUpdLag Declaration
\*---------------------------------------------------------------------------*/

class largeStrainUpdLag
:
        public mechanicsSubSolver
{
    // Private data

    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        largeStrainUpdLag(const largeStrainUpdLag&);

        //- Disallow default bitwise assignment
        void operator=(const largeStrainUpdLag&);

protected:

        //- Word to select the strain tensor definition
        word strainTensorName_;

        //- Increment of displacement field DD = D - D.oldTime()
        volVectorField DD_;

        //- Gradient of the displacement increment
        volTensorField gradDD_;

        //- Total deformation gradient
        volTensorField F_;

        //- Inverse total deformation gradient
        volTensorField Finv_;

        //- Jacobian of total deformation gradient
        volScalarField J_;

        //- Relative deformation gradient
        volTensorField relF_;

        //- Relative Jacobian: Jacobian of the relative deformation gradient
        volTensorField relFinv_;

        //- Relative Jacobian: Jacobian of the relative deformation gradient
        volScalarField relJ_;

        //- Reference mesh read from the initial points. The remaining elements 
        //- of the mesh (i.e. the connectivity) can be derived from the main 
        //- mesh.
        fvMesh referenceMesh_;
        pointIOField initialPoints_;

        //- Has the prevIter been stored?
        bool storedPrevIter_;

public:

    //- Runtime type information
    TypeName("largeStrainUpdLag");
        
    // Constructors

        //- Construct from mesh, materials, rheoloy and dictionary
        largeStrainUpdLag
        (
            fvMesh& mesh,
            const materials& materials,
            rheology& rheo,
            const dictionary& mechanicsDict
        );    

    //- Destructor
     ~largeStrainUpdLag();

     //- Access

        
    //- Return non-const displacement field
        inline volVectorField& DD() 
        {
            return DD_;
        }
        
        //- Return total or incremental displacement depending on
        //- which one should be used to move mesh points (useful when 
        //- not knowing if the instance is an updateMesh solver or not
        virtual inline volVectorField& DorDD() 
        {
            return DD_;
        }

        //- Return total or incremental displacement depending on
        //- which one should be used to move mesh points (useful when 
        //- not knowing if the instance is an updateMesh solver or not
        virtual inline const volVectorField& DorDD() const
        {
            return DD_;
        }

    // Member Functions
                
        //- Correct/update the properties
        virtual void correct();          

        //- Update mesh
        virtual void updateMesh()
        {
            //- Set to false so that mesh.oldPoints() gives the current points
            mesh_.moving(false);

            mechanicsSubSolver::updateMesh
            (
                mesh_, 
                DD_,
                !rheo_.modifiedPlaneStrain()
            );
        }; 

        //- Return true if converged
        virtual bool converged()
        {
            return mechanicsSubSolver::converged(DD_);
        };
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
