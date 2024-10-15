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
    Foam::misesPlasticity

Description
    Linear elastic-misesPlasticity law.

   \htmlinclude constitutiveLaw_misesPlasticity.html

SourceFiles
    misesPlasticity.C

\todo
    Add a correctAxialStrain?

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina, M.Reymond - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef misesPlasticity_H
#define misesPlasticity_H

#include "constitutiveLaw.H"
#include "yieldStressModel.H"
#include "Table.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class misesPlasticity Declaration
\*---------------------------------------------------------------------------*/

class misesPlasticity
:
    public constitutiveLaw
{

    //-autoPtr to yieldStressModel
    autoPtr<yieldStressModel> yieldStressModel_;

protected:

    // protected data    
        
        //- Plastic Strain Tensor
        volSymmTensorField& epsilonP_;        
        
        //- Incremental Plastic Strain tensor
        volSymmTensorField& DEpsilonP_;   

        //- yield stress field
        volScalarField& sigmaY_;

        //- Incrment of sigmaY field
        volScalarField& DSigmaY_;
       
        //- Equivalent plastic strain
        volScalarField& epsilonPEq_;

        //- Equivalent plastic strain increment
        volScalarField& DEpsilonPEq_;

        //- Plastic multiplier increment - plastric strain scaled by sqrt(2/3)
        volScalarField& DLambda_;

        //- Active yielding flag
        //     1.0 for active yielding
        //     0.0 otherwise 
        volScalarField& activeYield_;

        //- plasticN is the return direction to the yield surface
        volSymmTensorField& plasticN_;
        
        //- Linear plastic modulus. Defined by two points on the stress-plastic strain curve
        scalar Hp_;
        
    // protected Member Functions
        
        //- Disallow default bitwise copy construct
        misesPlasticity(const misesPlasticity&);

        //- Disallow default bitwise assignment
        void operator=(const misesPlasticity&);

public:
 
    //- Runtime type information
    TypeName("misesPlasticity");

    // Constructors

        //- Construct from mesh, global options, materials , dict and labelList
        misesPlasticity
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );


    //- Destructor
    ~misesPlasticity();


    // Member Functions

        //- Update hydrostatic and deviatoric stress    
        virtual void correct
        (
            volScalarField& sigmaHyd, 
            volSymmTensorField& sigmaDev, 
            volSymmTensorField& epsEl,
            const labelList& addr
        );  

        //- Correct the elastic strain tesnor removing oldTime copmonents of plastic
        // strains  
        virtual void correctEpsilonEl
        (
            volSymmTensorField& epsilonEl,
            const labelList& addr
        );   

        //- Update the yield stress: called at end of time-step
        virtual void updateTotalFields
        (
            const labelList& addr
        );   

        //- correct the additional strain field
        virtual void correctAdditionalStrain
        (
            volSymmTensorField& additionalStrain,
            const labelList& addr
        );       

        //- Access functions     
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
