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
    Foam::misesPlasticCreep

Description
    Linear elastic-misesPlasticity-creep law.

    \htmlinclude constitutiveLaw_misesPlasticCreep.html

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef misesPlasticCreep_H
#define misesPlasticCreep_H

#include "misesPlasticity.H"
#include "creepModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class "misesPlasticCreep" Declaration
\*---------------------------------------------------------------------------*/

class misesPlasticCreep
:
    public misesPlasticity
{
    // Private data    

        //- autoPtr to creep model
        autoPtr<creepModel> creepModel_;


    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        misesPlasticCreep(const misesPlasticCreep&);

        //- Disallow default bitwise assignment
        void operator=(const misesPlasticCreep&);

public:
 
    //- Runtime type information
    TypeName("misesPlasticCreep");


    // Constructors

        //- Construct from mesh, materials , dict and labelList
        misesPlasticCreep
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    ~misesPlasticCreep();


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

        //- correctAxialStrain in case of plane stress
        virtual void correctAxialStrain
        (
            volSymmTensorField& epsilon,
            const labelList& addr
        ); 

        //- correct the additional strain field
        virtual void correctAdditionalStrain
        (
            volSymmTensorField& additionalStrain,
            const labelList& addr
        ); 

        //- Return next deltaT
        virtual scalar nextDeltaT
        (
            const labelList& addr
        );         
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
