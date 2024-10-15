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
    Foam::coolantPressureFixedDisplacementZeroShearFvPatchVectorField

Description
    Coolant pressure patch-field of fixed-displacement-zeroShear type.
    It is designed for the top cap outer or bottom cap outer patches.

    The component of the displacement normal to the patch is derived from the
    force balance between internal stresses and the fluid pressure.

    A plane strain (epsilonZ costant on the patch) or plane stress approximations
    (sigmaZ costant on the patch) are applied.

    The fluid pressure can be provided as a fixed-value or as a time-dependent
    list.

Usage
    \verbatim
    cladOuter
    {
        type            coolantPressureFixedDisplacementZeroShear;
        planeStrain     on;
        
        // coolantPressure 1e5;
        coolantPressureList
        {
            type                table;

            // Instead of "file", the user can insert "values"
            file                "$FOAM_CASE/constant/lists/topCladdingRingPressureCoolantList";
    
            // values
            // (
            //     (0.0 1e5)
            //     (100.0 1e6)
            // );

            format              foam;      // data format (optional)
            outOfBounds         clamp;     // optional out-of-bounds handling
            interpolationScheme linear;    // optional interpolation method
        }
        value           $internalField;
    }    
    \endverbatim  

\todo
    Link this patch to the plenum spring force. 
    Find a way to automatically detect if the patch is the cap or the top ring

SourceFiles
    coolantPressureFixedDisplacementZeroShearFvPatchVectorField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)


\*---------------------------------------------------------------------------*/

#ifndef coolantPressureFixedDisplacementZeroShearFvPatchVectorField_H
#define coolantPressureFixedDisplacementZeroShearFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "Switch.H"
#include "Table.H"

#ifdef OPENFOAMESI
    #include "Function1.H"
#endif       

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class coolantPressureFixedDisplacementZeroShearFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class coolantPressureFixedDisplacementZeroShearFvPatchVectorField
:
    public fixedDisplacementZeroShearFvPatchVectorField
{
    // Private Data

protected:

    // Protected Data

        //- Fluid pressure
        scalarField coolantPressure_;

#ifdef OPENFOAMFOUNDATION
        //- Fluid pressure time-dependent list
        autoPtr<Function1s::Table<scalar>> coolantPressureList_;
#elif OPENFOAMESI
        //- Fluid pressure time-dependent list
        autoPtr<Function1<scalar>> coolantPressureList_;
#endif   

        //- Apply plane strain approximation to deduce axial displacement
        Switch planeStrain_;

        //- Name of strain field
        word displacementName_;   

public:

    //- Runtime type information
    TypeName("coolantPressureFixedDisplacementZeroShear");


    // Constructors

        //- Construct from patch and internal field
        coolantPressureFixedDisplacementZeroShearFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        coolantPressureFixedDisplacementZeroShearFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // coolantPressureFixedDisplacementZeroShearFvPatchVectorField
        //  onto a new patch
        coolantPressureFixedDisplacementZeroShearFvPatchVectorField
        (
            const coolantPressureFixedDisplacementZeroShearFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        coolantPressureFixedDisplacementZeroShearFvPatchVectorField
        (
            const coolantPressureFixedDisplacementZeroShearFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new coolantPressureFixedDisplacementZeroShearFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~coolantPressureFixedDisplacementZeroShearFvPatchVectorField()
        {}


    // Member functions

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
