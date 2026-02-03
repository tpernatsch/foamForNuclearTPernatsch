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
    Foam::smallStrain

Description
    This class solves for total displacement `D` in a total Lagrangian configuration 
    (i.e. the mesh is NOT updated), and employs the small-strain approximation.

    This means that the strain tensor is linear, i.e.:
	    
	$$
	\begin{aligned}
	\epsilon = \frac{1}{2} (\nabla D + \nabla^T D)
	\end{aligned}
	$$

    and that the geometric non-linearity is neglected, i.e. the reference and updated 
    configurations are considered to be the same. 
	
________________________________________________________________________________

\par Options	
	As a `mechanicsSubSolver` class, `smallStrain` **allows** the user to specify
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
        mechanicsSolver smallStrain;

        mechanicsOptions
        {
            // Print summary of boundary-forces at each iteration
            forceSummary false;
            
            // Creates cylindrical stress and displacement field
            cylindricalStress false;
            
            // Creates spherical stress and displacement field
            sphericalStress false;
            
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
    smallStrain.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef smallStrain_H
#define smallStrain_H

#include "mechanicsSubSolver.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class smallStrain Declaration
\*---------------------------------------------------------------------------*/

class smallStrain
:
        public mechanicsSubSolver
{
    // Private data

    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        smallStrain(const smallStrain&);

        //- Disallow default bitwise assignment
        void operator=(const smallStrain&);

protected:

        //- Has the prevIter been stored?
        bool storedPrevIter_;

public:

    //- Runtime type information
    TypeName("smallStrain");
        
    // Constructors

        //- Construct from mesh, materials, rheology and dictionary
        smallStrain
        (
            fvMesh& mesh,
            const materials& materials,
            rheology& rheo,
            const dictionary& mechanicsDict
        );    

    //- Destructor
     ~smallStrain();


    // Member Functions
                
        //- Correct/update the properties
        virtual void correct();      
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
