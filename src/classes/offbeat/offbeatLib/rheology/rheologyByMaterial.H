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
    Foam::rheologyByMaterial

Description
    Rheology laws are specified on each material subdictionary.
    
    For more information see [here](rheology.html).

SourceFiles
    rheologyByMaterial.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef rheologyByMaterial_H
#define rheologyByMaterial_H

#include "rheology.H"
#include "sliceMapper.H"
#include "Table.H"

#ifdef OPENFOAMESI
    #include "Function1.H"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class rheologyByMaterial Declaration
\*---------------------------------------------------------------------------*/

class rheologyByMaterial
:
    public rheology
{
    // Private data
        
        //- Precision threshold for spring force balance 
        scalar precision_;

        //- Spring modulus
        scalar springModulus_;

        //- Spring initial loading (in m)
        scalar springPreCompression_;

        //- Pointer to slice mapper (necessary for modifiedPlaneStrain)
        const sliceMapper* mapper_;

#ifdef OPENFOAMFOUNDATION
        //- Fluid pressure time-dependent list (necessary for modifiedPlaneStrain)
        autoPtr<Function1s::Table<scalar>> fluidPressureList_;
#elif OPENFOAMESI
        //- Fluid pressure time-dependent list (necessary for modifiedPlaneStrain)
        autoPtr<Function1<scalar>> fluidPressureList_;
#endif

        //- SpringForce (only for modifiedPlaneStrain)
        autoPtr<scalarField> springForce_;

        //- Pressure force (only for modifiedPlaneStrain)
        autoPtr<scalarField> pressureForce_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        rheologyByMaterial(const rheologyByMaterial&);

        //- Disallow default bitwise assignment
        void operator=(const rheologyByMaterial&);

public:

    //- Runtime type information
    TypeName("byMaterial");
    

    // Constructors

        //- Construct from mesh, materials and dictionary
        rheologyByMaterial
        (
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& mechanicsDict
        ); 
    

    //- Destructor
    ~rheologyByMaterial();


    // Member Functions

        //- Correct/update the properties
        virtual void correct
        (
            volSymmTensorField& sigma, 
            volSymmTensorField& epsilon,
            volVectorField& D
        );     

        //- Correct out of plane strain for planeStress cases
        virtual void correctPlaneStress
        (
            volSymmTensorField& epsilon
        );      

        //- Correct out of plane strain for modifiedPlaneStress cases
        virtual void correctModifiedPlaneStrain
        (
            volSymmTensorField& epsilon,
            volVectorField& D
        );     

        //- Update total fields at the end of time step
        virtual void updateTotalFields();

        //- Collect additional strain field from zone constitutive laws
        volSymmTensorField getLawsAdditionalStrain();

        //- Return next deltaT according to rheology
        virtual scalar nextDeltaT(); 

        //- Return next macroTime according to rheology
        virtual scalar nextMacroTime(); 
          
        //- Return maximum time according to rheology
        virtual scalar lastMacroTime();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
