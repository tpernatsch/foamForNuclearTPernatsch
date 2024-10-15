/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::neoHookeanMisesPlasticity

Description
    Rheology class for neo-Hookean hyperelasticity.
    
    \htmlinclude constitutiveLaw_neoHookeanElasticity.html


SourceFiles
    neoHookeanMisesPlasticity.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef neoHookeanMisesPlasticity_H
#define neoHookeanMisesPlasticity_H

#include "constitutiveLaw.H"
#include "Table.H"

#ifdef OPENFOAMESI
    #include "Function1.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class neoHookeanMisesPlasticity Declaration
\*---------------------------------------------------------------------------*/

class neoHookeanMisesPlasticity
:
    public constitutiveLaw
{

protected:

    // Protected data

        //- Cauchy yield stress
        volScalarField& sigmaY_;

        //- Incremental change in sigmaY
        volScalarField& DSigmaY_;

        //- Total plastic strain
        volSymmTensorField& epsilonP_;

        //- Incremental change of plastic strain
        volSymmTensorField& DEpsilonP_;

        //- Equivalent plastic strain
        volScalarField& epsilonPEq_;

        //- Equivalent plastic strain increment
        volScalarField& DEpsilonPEq_;

        //- Plastic multiplier increment - plastric strain scaled by sqrt(2/3)
        volScalarField& DLambda_;

        //- Elastic left Cauchy-Green trial strain tensor with volumetric term
        //  removed. Trial means that plasticity is neglected.
        volSymmTensorField& bEbarTrial_;  

        //- Elastic left Cauchy-Green strain tensor with volumetric term removed
        volSymmTensorField& bEbar_;

        //- plasticN is the return direction to the yield surface
        volSymmTensorField& plasticN_;

        //- Active yielding flag
        //     1.0 for active yielding
        //     0.0 otherwise
        volScalarField& activeYield_;

        //- Linear plastic modulus. It is only used if plasticity is linear,
        //  defined by two points on the stress-plastic strain curve
        scalar Hp_;

#ifdef OPENFOAMFOUNDATION
        //- Table of post-yield stress versus plastic strain
        autoPtr<Function1s::Table<scalar>> stressPlasticStrainSeries_;
#elif OPENFOAMESI
        typedef Foam::Function1Types::Table<scalar> Table;

        //- Table of post-yield stress versus plastic strain
        autoPtr<Table> stressPlasticStrainSeries_;
#endif  

    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        neoHookeanMisesPlasticity(const neoHookeanMisesPlasticity&);

        //- Disallow default bitwise assignment
        void operator=(const neoHookeanMisesPlasticity&);


public:
 
    //- Runtime type information
    TypeName("neoHookeanMisesPlasticity");

    // Constructors

        //- Construct from mesh, materials , dict and labelList
        neoHookeanMisesPlasticity
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );

    //- Destructor
    ~neoHookeanMisesPlasticity();


    // Member Functions

        //- Update mechanical law     
        virtual void correct
        (
            volScalarField& sigmaHyd, 
            volSymmTensorField& sigmaDev, 
            volSymmTensorField& epsEl,
            const labelList& addr
        );

        //- Update stress values according to rheology law     
        virtual void updateStress
        (
            volSymmTensorField& sigma, 
            volScalarField& sigmaHyd, 
            volSymmTensorField& sigmaDev,
            const labelList& addr
        );
        
        //- Update the yield stress: called at end of time-step
        virtual void updateTotalFields
        (
            const labelList& addr
        );   

        virtual void correctAdditionalStrain
        (
            volSymmTensorField& additionalStrain,
            const labelList& addr
        ); 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
