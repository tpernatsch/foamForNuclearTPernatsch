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
    Foam::coolantPressureFixedDisplacementFvPatchVectorField
 
Description
    Coolant pressure patch-field of fixed-displacement type.
    It is designed for the top cap outer or bottom cap outer patches.

    The component of the displacement normal to the patch is derived from the
    force balance between internal stresses and the coolant pressure.
    A plane strain (epsilonZ costant on the patch) or plane stress approximations
    (sigmaZ costant on the patch) are applied.
    A pre-loading can be provided to the spring.

    The fluid pressure can be provided as a fixed-value or as a time-dependent
    list.

Usage
    \verbatim
    cladTop
    {
        type            coolantPressureFixedDisplacement;
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
 
SourceFiles
     coolantPressureFixedDisplacementFvPatchVectorField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021
 
\*---------------------------------------------------------------------------*/
 
#ifndef coolantPressureFixedDisplacementFvPatchVectorField_H
#define coolantPressureFixedDisplacementFvPatchVectorField_H
 
#include "fixedDisplacementFvPatchVectorField.H"
#ifdef OPENFOAMESI
    #include "Function1.H"
#endif     
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
        Class coolantPressureFixedDisplacementFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/
 
class coolantPressureFixedDisplacementFvPatchVectorField
:
    public fixedDisplacementFvPatchVectorField
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

    // Protected Member Functions
 
public:
 
    //- Runtime type information
    TypeName("coolantPressureFixedDisplacement");
 
 
    // Constructors
 
        //- Construct from patch and internal field
        coolantPressureFixedDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );
 
        //- Construct from patch, internal field and dictionary
        coolantPressureFixedDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );
 
        //- Construct by mapping the given coolantPressureFixedDisplacementFvPatchVectorField
        //  onto a new patch
        coolantPressureFixedDisplacementFvPatchVectorField
        (
            const coolantPressureFixedDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );
 
 
        //- Copy constructor setting internal field reference
        coolantPressureFixedDisplacementFvPatchVectorField
        (
            const coolantPressureFixedDisplacementFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );
 
        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new coolantPressureFixedDisplacementFvPatchVectorField(*this, iF)
            );
        }
 
 
    // Member Functions

        //- Write
        virtual void write(Ostream&) const;
    
    // Evaluation functions
 
        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual void updateCoeffs ();

 };
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 } // End namespace Foam
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 #endif
 
 // ************************************************************************* //