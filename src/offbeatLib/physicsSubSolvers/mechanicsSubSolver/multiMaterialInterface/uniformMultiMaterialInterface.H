/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::uniformMultiMaterialInterface

Description
    This class applies a uniform multi-material correction to the solution of the 
    momentum balance equation in mechanicsSubSolver.

________________________________________________________________________________

\par Options
    The parameters for this class are to be set in the `multiMaterialCorrection` 
    sub-dictionary within the `mechanicsOptions` subdictionary, which in turn is
    located within the `solverDict`.

    Parameters in `multiMaterialCorrection`:

    - <b>`defaultWeights`</b> - Parameter of the `multiMaterialCorrection` subdictionary; 
    defines the fraction of the full correction applied, between 1 and 0. Applying full 
    correction in all directions for every cell may negatively impact convergence. **We 
    recommend setting `defaultWeights` between 0.9 and 1 for a good balance (1 if possible).**
    - <b>`defaultWeightsGrad`</b> - Defines the fraction of the full correction applied 
    to the gradient contribution in sigmaExp, ranging from 0 to 1. **By default it is set
    equal to defaultWeights**. Applying full correction in all directions for every cell may negatively impact convergence. **We 
    recommend setting `defaultWeights` between 0.9 and 1 for a good balance (1 if possible).**
    
________________________________________________________________________________

Usage
	Here is a code snippet of the `solverDict` to be used for activating the
	`uniform` multi-material correction:
    
    \verbatim
        // It is just an example, you can select any mechanicsSolver
        mechanicsSolver smallStrain;
        
        // Rest of the solverDict etc...

        mechanicsOptions
        {
            // Other options if needed

            // This is not necessary, as activating the multiMaterial correction
            // will make RhieChowCorrection false by default
            RhieChowCorrection false;    

            // The presence of this subdictionary activates the multiMaterial 
            // correction
            multiMaterialCorrection
            {
                type uniform;

                defaultWeights 1.0;

                // Not needed in all cases, by default equal 
                // to defaultWeights
                defaultWeightsGrad 1.0;
            }
        }
        
        // Rest of the solverDict etc...
    \endverbatim

SourceFiles
    uniformMultiMaterialInterface.C

\*---------------------------------------------------------------------------*/

#ifndef uniformMultiMaterialInterface_H
#define uniformMultiMaterialInterface_H

#include "multiMaterialInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class uniformMultiMaterialInterface Declaration
\*---------------------------------------------------------------------------*/

class uniformMultiMaterialInterface
:
    public multiMaterialInterface
{
    // Private data
    
        //- Default weight to apply uniformly for the
        //  additional component of the strain
        scalar defaultWeightEps_;
    
        //- Default weight to apply uniformly for the
        //  explicit components of the twoMuLambda*gradD term
        scalar defaultWeightGrad_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        uniformMultiMaterialInterface(const uniformMultiMaterialInterface&);

        //- Disallow default bitwise assignment
        void operator=(const uniformMultiMaterialInterface&);


public:

    //- Runtime type information
    TypeName("uniform");


    // Constructors

        //- Construct from Dictionary
        uniformMultiMaterialInterface
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    ~uniformMultiMaterialInterface();


    // Member Functions

        //- Return true if this interface applies a multi-material correction
        virtual bool active() const
        {
            return true;
        }
        
        //- Correct the face weights
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
