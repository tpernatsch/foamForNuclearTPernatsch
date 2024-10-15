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
    Foam::plenumSpringFixedDisplacementZeroShearFvPatchVectorField

Description
    Plenum spring patchfield for the top of the fuel column (fixedDisplacement-
    zeroShear version).

    The component of the displacement normal to the patch is derived from the
    force balance between internal stresses, spring force and plenum pressure.
    A plane strain (epsilonZ costant on the patch) or plane stress approximations
    (sigmaZ costant on the patch) are applied.
    The shear traction is set to zero.
    A pre-loading can be provided to the spring.

Usage
    \verbatim
    cladOuter
    {
        type            plenumSpringFixedDisplacementZeroShear;
        planeStrain     on;

        springModulus   1e4;
        initialSpringLoading 0;

        fuelTopPatches (fuelTop);
        topCapInnerPatches (cladOuter);

        value           $internalField;
    }    
    \endverbatim    

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)


\*---------------------------------------------------------------------------*/

#ifndef plenumSpringFixedDisplacementZeroShearFvPatchVectorField_H
#define plenumSpringFixedDisplacementZeroShearFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "plenumSpringBase.H"
#include "Switch.H"
#include "Table.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class plenumSpringFixedDisplacementZeroShearFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class plenumSpringFixedDisplacementZeroShearFvPatchVectorField
:
    public fixedDisplacementZeroShearFvPatchVectorField,
    public plenumSpringBase
{
    // Private data

protected:

    // Protected data

        //- Apply plane strain approximation to deduce axial displacement
        Switch planeStrain_;


public:

    //- Runtime type information
    TypeName("plenumSpringFixedDisplacementZeroShear");


    // Constructors

        //- Construct from patch and internal field
        plenumSpringFixedDisplacementZeroShearFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        plenumSpringFixedDisplacementZeroShearFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // plenumSpringFixedDisplacementZeroShearFvPatchVectorField
        //  onto a new patch
        plenumSpringFixedDisplacementZeroShearFvPatchVectorField
        (
            const plenumSpringFixedDisplacementZeroShearFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

    // Destructor

        virtual ~plenumSpringFixedDisplacementZeroShearFvPatchVectorField()
        {}


    // Member functions

        // Access Functions
        

        // Mapping functions


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
