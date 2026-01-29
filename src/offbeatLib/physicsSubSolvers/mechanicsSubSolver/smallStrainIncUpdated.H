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
    Foam::smallStrainIncUpdated

Description
    This class olves for incremental displacement `DD` in an updated configuration 
    (i.e. the mesh is updated at the start of the time step), but employs the 
    small-strain approximation.

    This means that the strain tensor is still assumed to be linear linear, i.e.:
	    
	$$
	\begin{aligned}
	\epsilon = \frac{1}{2} (\nabla D + \nabla^T D)
	\end{aligned}
	$$

    and that while the configuration is updated at every time step, the geometric non-linearity 
    between start of time-step and end-of time step configuration is neglected, i.e. the 
    start of time-step or reference configuration and updated or end of time-step
    configuration are considered to be the same (within that time-step). 
	
________________________________________________________________________________

\par Options	
	As a `mechanicsSubSolver` class, `smallStrainIncrementalUpdated` **allows** the user to 
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
    - <b>`updateMesh`</b> - Decides wheter the configuration is updated at every time-step
    or not. **It is set to true by default**. If set to false, the solver becomes
    effectively a small strain incremental solver (with no update).
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

________________________________________________________________________________

Usage
	Here is a code snippet of the `solverDict` to be used for activating the
	`smallStrainIncrementalUpdated` mechanics solver class:
    
    \verbatim
        mechanicsSolver smallStrainIncrementalUpdated;

        mechanicsOptions
        {
            // Print summary of boundary-forces at each iteration
            forceSummary false;
            
            // Creates cylindrical stress and displacement field
            cylindricalStress false;
            
            // Creates spherical stress and displacement field
            sphericalStress false;

            // Update the mesh at every time-step
            updateMesh on;
            
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
        }
    \endverbatim

SourceFiles
    smallStrainIncUpdated.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef smallStrainIncUpdated_H
#define smallStrainIncUpdated_H

#include "mechanicsSubSolver.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                Class smallStrainIncUpdated Declaration
\*---------------------------------------------------------------------------*/

class smallStrainIncUpdated
:
        public mechanicsSubSolver
{
    // Private data

    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        smallStrainIncUpdated(const smallStrainIncUpdated&);

        //- Disallow default bitwise assignment
        void operator=(const smallStrainIncUpdated&);

protected:
        //- Increment of displacement field DD = D - D.oldTime()
        volVectorField DD_;

        //- Gradient of the displacement increment
        volTensorField gradDD_;

        //- Incremental stress tensor
        //- TODO: it should be moved to rheology
        volSymmTensorField DSigma_;

        //- Incremental explicit component of the stress tensor
        //- TODO: it should be moved to rheology
        volTensorField DSigmaExp_;

        //- Incremental divergence of the explicit component of the stress tensor
        //- TODO: it should be moved to rheology
        volVectorField divDSigmaExp_;

        //- Final residual of latest inner iteration
        vector residualDD_;

        //- Final residual of firs inner loop iteration
        vector initialResidualDD_;

        //- Update mesh or not
        const Switch updateMesh_;

        //- Reference mesh read from the initial points. The remaining elements 
        //- of the mesh (i.e. the connectivity) can be derived from the main 
        //- mesh.
        autoPtr<fvMesh> referenceMesh_;

        //- Has the prevIter been stored?
        bool storedPrevIter_;

public:

    //- Runtime type information
    TypeName("smallStrainIncrementalUpdated");
        
    // Constructors

        //- Construct from mesh, materials, rheoloy and dictionary
        smallStrainIncUpdated
        (
            fvMesh& mesh,
            const materials& materials,
            rheology& rheo,
            const dictionary& mechanicsDict
        );    

    //- Destructor
     ~smallStrainIncUpdated();

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
            if(referenceMesh_.valid())
            {
                return DD_;
            }
            else
            {
                return D_;
            }
        }
    
        //- Return total or incremental displacement depending on
        //- which one should be used to move mesh points (useful when 
        //- not knowing if the instance is an updateMesh solver or not)
        virtual inline const volVectorField& DorDD() const
        {
            if(referenceMesh_.valid())
            {
                return DD_;
            }
            else
            {
                return D_;
            }
        }


    // Member Functions
        
        //- Update total fields at the end of time step
        virtual void updateIncrementalFields();
                
        //- Correct/update the properties
        virtual void correct();     

        //- Update mesh
        virtual void updateMesh()
        {
            if(referenceMesh_.valid())
            {
                //- Set to false so that mesh.oldPoints() gives the current points
                mesh_.moving(false);

                mechanicsSubSolver::updateMesh
                (
                    mesh_, 
                    DD_,
                    !rheo_.modifiedPlaneStrain()
                );
            }
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
