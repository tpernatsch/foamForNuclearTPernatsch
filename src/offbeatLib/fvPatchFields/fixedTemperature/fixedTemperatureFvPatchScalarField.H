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
    Foam::fixedTemperatureFvPatchScalarField

Description
    The `fixedTemperature` boundary condition allows to specify a constant fixed
    temperature on a patch.

    It is similar to the standard `fixedValue` OpenFOAM fvPatchField but, in case
    an oxide layer is present, an additional thermal resistance corresponding to
    the oxide is considered.

________________________________________________________________________________

\par Options
    The `fixedTemperature` fvPatchField can be selected in the patch subdictionary
    inside the `boundaryField` subdictionary of the temperature field.

    Parameters in the patch subdictionary for `fixedTemperature`:

    - <b>`outerOxideTemperature`</b> - If the corrosion model is present, this keyword
    specifies the initial outer temperature (i.e. the one outside the oxide layer).
    **If not present, the initial oxide outer temperature is set to the one specified
    in  `value`.**
    - <b>`value`</b> - The initial **temperature** value (in case of oxidation model, it is
    the outer temperature of the metallic portion of the body).

________________________________________________________________________________

Usage
    To apply the `fixedTemperature` boundary condition to a given patch, the
    following example can be used as a template:
    
    \verbatim
        claddingOuterSurface
        {
            type            fixedTemperature;

            // "value" is:
            // - the initial metal-oxide interface temperature, if the corrosion
            // model is present
            // - the initial outer temperature, if the corrosion is not present
            value           uniform 300;
                
            outerOxideTemperature           uniform 300;
        }
    \endverbatim	    

SourceFiles
     fixedTemperatureFvPatchScalarField.C

\mainauthor
    E. Brunetto, A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021
 
\*---------------------------------------------------------------------------*/
 
#ifndef fixedTemperatureFvPatchScalarField_H
#define fixedTemperatureFvPatchScalarField_H
 
#include "fixedValueFvPatchField.H"
#include "axialProfile.H"
#include "InterpolateTables.H"

#include "primitiveFields.H"
#include "objectRegistry.H"
#include "fvMesh.H"
#include "volFields.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
        Class fixedTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/
 
class fixedTemperatureFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
    // Private Data

protected:

    // Protected Data

        //- Temperature field imposed at the boundary
        //  NOTE: the distinction between "oxideTemperature" and "value" becomes
        //  necessary when an oxide layer is present. In this case, "value" is 
        //  the temperature of the oxide/material interface, while "oxideTemperature"
        //  is the outer oxide temperature. If the oxide is absent, the two are
        //  the same, i.e. the outer material temperature. This might change
        //  if we include the oxide layer in the model explicitly.
        mutable scalarField T_;
    
        //- Name of conductivity field
        word kappaName_;

    // Protected Member Functions

        //- Updates surface temperature (includes oxide treatement if present)
        virtual void updateSurfaceTemperature();
 
public:
 
    //- Runtime type information
    TypeName("fixedTemperature");
 
 
    // Constructors
 
        //- Construct from patch and internal field
        fixedTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );
 
        //- Construct from patch, internal field and dictionary
        fixedTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );
 
        //- Construct by mapping the given fixedTemperatureFvPatchScalarField
        //  onto a new patch
        fixedTemperatureFvPatchScalarField
        (
            const fixedTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );
 
 
        //- Copy constructor setting internal field reference
        fixedTemperatureFvPatchScalarField
        (
            const fixedTemperatureFvPatchScalarField&,
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
                new fixedTemperatureFvPatchScalarField(*this, iF)
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

        // I-O

            //- Write
            virtual void write(Ostream&) const;
    
        // Evaluation functions
     
            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual void updateCoeffs();
 };
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 } // End namespace Foam
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 #endif
 
 // ************************************************************************* //