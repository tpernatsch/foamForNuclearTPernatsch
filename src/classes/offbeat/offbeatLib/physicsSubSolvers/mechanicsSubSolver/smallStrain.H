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
    Small-strain mechanics sub-solver.
    
    It solves for total displacement 'D' in a total Lagrangian 
    configuration (i.e. the mesh is NOT updated). 

    The geometric non-linearity (i.e. the inter-dependence between displacement
    and geometric configuration) is not taken care of following the small-strain
    approximation: reference and updated configuration are considered to be the
    same.

    The strain tensor is linear:
    \f[
                        epsilon = 0.5*(gradD + gradD.T())
    \f]
     
Usage
    In solverDict file:
    \verbatim
    mechanicsSolver smallStrain;

    mechanicsOptions
    {
        // Print summary of boundary-forces at each iteration
        forceSummary false;
        
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
        vectorField DispPrev;
        vectorField DispPrevPrev;        
        Switch inertialRelaxation_;
        scalar inertialRelaxationDeltaT_;

    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        smallStrain(const smallStrain&);

        //- Disallow default bitwise assignment
        void operator=(const smallStrain&);

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
