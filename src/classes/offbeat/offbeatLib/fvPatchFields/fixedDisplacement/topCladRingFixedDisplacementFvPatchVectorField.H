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
     Foam::topCladRingFixedDisplacementFvPatchVectorField
 
Description
    Patchfield for the top cladding ring, to be used when the top cap is not 
    explicitly part of the model (i.e. the cladding tube in the model is 
    truncated). FixedDisplacement plane strain version.

    The component of the displacement normal to the patch is derived from the
    force balance between:
    - the fluid pressure
    - the gap pressure on the plenum side of the gap (acting on the opposite direction).
    - the fuel spring pressure (acting on the opposite direction).

    A plane strain (epsilonZ costant on the patch) or plane stress approximations
    (sigmaZ costant on the patch) are applied.
    
    The shear components of the displacement fields are set to zero.
    A pre-loading can be provided to the spring.

    The fluid pressure can be provided as a fixed-value or as a time-dependent
    list.

    NOTE: This patch assumes that the fuel plenum spring BC is applied to one 
    and only one patch

Usage
    \verbatim
    cladOuter
    {
        type            topCladRingSpringPlaneStrain;
        planeStrain     on;

        topCapRadius       0.00525;
        topCapPlenumRadius 0.0045;

        // fluidPressure 1e5;
        fluidPressureList
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
    topCladRingFixedDisplacementFvPatchVectorField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021
 
\*---------------------------------------------------------------------------*/
 
#ifndef topCladRingFixedDisplacementFvPatchVectorField_H
#define topCladRingFixedDisplacementFvPatchVectorField_H
 
#include "coolantPressureFixedDisplacementFvPatchVectorField.H"
#include "plenumSpringBase.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
        Class topCladRingFixedDisplacementFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/
 
class topCladRingFixedDisplacementFvPatchVectorField
:
    public coolantPressureFixedDisplacementFvPatchVectorField
{
    // Private Data

        //- Top cap radius on the fluid side
        scalar topCapRadius_;

        //- Top cap radius on the plenum side
        scalar topCapPlenumRadius_;

        //- Label of fuel top patch
        label fuelTopPatchID_;

        //- Apply plane strain approximation to deduce axial displacement
        Switch planeStrain_;

        //- Name of strain field
        word displacementName_;   

protected:

    // Protected Data

    // Protected Member Functions
 
public:
 
    //- Runtime type information
    TypeName("topCladRingFixedDisplacement");
 
 
    // Constructors
 
        //- Construct from patch and internal field
        topCladRingFixedDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );
 
        //- Construct from patch, internal field and dictionary
        topCladRingFixedDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );
 
        //- Construct by mapping the given topCladRingFixedDisplacementFvPatchVectorField
        //  onto a new patch
        topCladRingFixedDisplacementFvPatchVectorField
        (
            const topCladRingFixedDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );
 
 
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