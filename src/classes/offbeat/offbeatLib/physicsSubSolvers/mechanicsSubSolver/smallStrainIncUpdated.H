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
    Mechanics sub-solver solving for incremental displacement field. Small 
    strain approximation is used for the strain tensor definition. The mesh is 
    updated once per time-step.

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

        //- Final residual of latest inner iteration
        vector residualDD_;

        //- Final residual of firs inner loop iteration
        vector initialResidualDD_;

        //- Reference mesh read from the initial points. The remaining elements 
        //- of the mesh (i.e. the connectivity) can be derived from the main 
        //- mesh.
        fvMesh referenceMesh_;

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
