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
    Foam::largeStrainTotLag

Description
    Finite- (or large-) strain mechanics sub-solver.
    
    It solves for total displacement 'D' in a total Lagrangian configuration 
    (i.e. the mesh is not updated). 

    The geometric non-linearity (i.e. the inter-dependence between displacement
    and geometric configuration) is taken care of by multiplying the stress 
    tensor (Cauchy) by J*F.inv(), where F is the deformation gradient and J is 
    the determinant of F.

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
    mechanicsSolver largeStrainTotLag;

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
    largeStrainTotLag.C

\mainauthor
    E. Brunetto, A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef largeStrainTotLag_H
#define largeStrainTotLag_H

#include "mechanicsSubSolver.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class largeStrainTotLag Declaration
\*---------------------------------------------------------------------------*/

class largeStrainTotLag
:
        public mechanicsSubSolver
{
    // Private data     

        vectorField DispPrev;

        vectorField DispPrevPrev;

        //- Word to select the strain tensor definition
        word strainTensorName_;


        //- Total deformation gradient
        volTensorField F_;

        //- Inverse of the total deformation gradient
        volTensorField Finv_;

        //- Relative deformation gradient
        //- Notice: added field here even if not used by the mechanical 
        //- subsolver because it is needed by some constitutive laws 
        //- (e.g. neoHookeanMisesPlasticity)
        volTensorField relF_;

        //- Jacobian of total deformation gradient
        volScalarField J_;

    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        largeStrainTotLag(const largeStrainTotLag&);

        //- Disallow default bitwise assignment
        void operator=(const largeStrainTotLag&);

public:

    //- Runtime type information
    TypeName("largeStrainTotLag");
        
    // Constructors

        //- Construct from mesh, materials, rheology and dictionary
        largeStrainTotLag
        (
            fvMesh& mesh,
            const materials& materials,
            rheology& rheo,
            const dictionary& mechanicsDict
        );    

    //- Destructor
     ~largeStrainTotLag();


    // Member Functions
                
        //- Correct/update the properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
