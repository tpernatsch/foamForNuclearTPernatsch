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
    Finite- (or large-) strain mechanics sub-solver.
    
    It solves for incremental displacement 'DD' in an updated Lagrangian 
    configuration (i.e. the mesh is updated at the start of the time step). 

    The geometric non-linearity (i.e. the inter-dependence between displacement
    and geometric configuration) is taken care of by multiplying the stress 
    tensor (Cauchy) by relJ*relF.inv(), where relF is the relative deformation 
    gradient and relJ is its determinant.

    By default, the strain tensor takes into account the second order deformations:
    \f[
                        epsilon = E = 0.5*(gradD + gradD.T() + gradD*gradD.T())
    \f]
     
    However, with the "linearStrainTensor" keyword, the user can decide to force
    a linear strain tensor:
    \f[
                        epsilon = 0.5*(gradD + gradD.T())
    \f]
     
Usage
    In solverDict file:
    \verbatim
    mechanicsSolver largeStrainUpdLag;

    mechanicsOptions
    {
        // Print summary of boundary-forces at each iteration
        forceSummary false;

        // Force linear strain tensor
        linearStrainTensor false;
        
        // Creates cylindrical stress and displacement field (for debugging)
        cylindricalStress false;
        
        // Apply RhieChowCorrection to stabilize momentum equation (it should
        // be activated in most situations)
        RhieChowCorrection true;
        
        // Additional subDict for multi-material cases
        multiMaterialCorrection
        {
            ...
        }
        
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
        vectorField DispPrev;
        vectorField DispPrevPrev;   

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
        
        //- Return total or incremental displacement (useful when not knowing 
        //  if the instance is a total or incremental solver)
        virtual inline volVectorField& DorDD() 
        {
            return DD_;
        }

        //- Return total or incremental displacement (useful when not knowing 
        //  if the instance is a total or incremental solver)
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
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
