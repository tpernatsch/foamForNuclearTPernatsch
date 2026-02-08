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
    Foam::resistiveLayerFvPatchScalarField

Description
    The `resistiveLayer` corrects the thermal conductivity boundary value
    to account for a thin resistive layer (oxide or crud) on the outer 
    surface.

________________________________________________________________________________

\par Options
    The `resistiveLayer` fvPatchField can be selected in the patch subdictionary
    inside the `boundaryField` subdictionary of the temperature field.

    Parameters in the patch subdictionary for `fixedTemperature`:

    - <b>`alpha`</b> - Outer layer heat transfer coefficient (W/m2.K). 
      If not present, zero is assumed.
    - <b>`value`</b> - The initial value, which will be overwritten when
      the boundary is updated.

________________________________________________________________________________

Usage
    To apply the `resistiveLayer` boundary condition to a given patch, the
    following example can be used as a template:
    
    \verbatim
        claddingOuterSurface
        {
            type            resistiveLayer;
            value           uniform 300;
            alpha           uniform 1e6;
        }
    \endverbatim        

SourceFiles
     resistiveLayerFvPatchScalarField.C

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    September, 2025
 
\*---------------------------------------------------------------------------*/
 
#ifndef resistiveLayerFvPatchScalarField_H
#define resistiveLayerFvPatchScalarField_H
 
#include "fixedValueFvPatchFields.H"

#include "primitiveFields.H"
#include "objectRegistry.H"
#include "fvMesh.H"
#include "volFields.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
        Class resistiveLayerFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/
 
class resistiveLayerFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Member functions
    
        //- Check that this BC is used for thermal conductivity
        void check() const;
        
        //- Correct the boundary conductivity
        void applyCorrection();

protected:

    // Protected Data

        //- Heat transfer coefficient of resistive layer
        mutable scalarField alpha_;
 
public:
 
    //- Runtime type information
    TypeName("resistiveLayer");
 
 
    // Constructors
 
        //- Construct from patch and internal field
        resistiveLayerFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );
 
        //- Construct from patch, internal field and dictionary
        resistiveLayerFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );
 
        //- Construct by mapping the given resistiveLayerFvPatchScalarField
        //  onto a new patch
        resistiveLayerFvPatchScalarField
        (
            const resistiveLayerFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );
 
 
        //- Copy constructor setting internal field reference
        resistiveLayerFvPatchScalarField
        (
            const resistiveLayerFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );
 
        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new resistiveLayerFvPatchScalarField(*this, iF)
            );
        }
              
    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchField<scalar>&, const labelList&);
    
        // Update
        
            //- Return the surface resistance
            inline const scalarField& alpha() const
            {
                return alpha_;
            }

            //- Return the surface resistance
            inline scalarField& alpha()
            {
                return alpha_;
            }
        
            //- Set the boundary value and apply correction
            virtual void operator=(const UList<scalar>&);
     
            //- Set the boundary value and apply correction
            virtual void operator=(const fvPatchScalarField&);
     
            //- Set the boundary value and apply correction
            virtual void operator=(const scalar&);

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual void updateCoeffs();

        // I-O

            //- Write
            virtual void write(Ostream&) const;
};
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 } // End namespace Foam
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 #endif
 
 // ************************************************************************* //