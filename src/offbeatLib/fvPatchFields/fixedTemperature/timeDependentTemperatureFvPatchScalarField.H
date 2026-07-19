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
    Foam::timeDependentTemperatureFvPatchScalarField

Description
    The `timeDependentTemperature` boundary condition allows to impose a time
    dependent temperature on a patch.

    It is similar to the standard `uniformFixedValue` OpenFOAM fvPatchField but, 
    in case an oxide layer is present, an additional thermal resistance corresponding to
    the oxide is considered.

________________________________________________________________________________

\par Options
    The `timeDependentTemperature` fvPatchField can be selected in the patch subdictionary
    inside the `boundaryField` subdictionary of the temperature field.

    Parameters in the patch subdictionary for `timeDependentTemperature`:

    - <b>`temperatureSeries`</b> - Dictionary containing the information related to the
    time-dependent temperature.
    - <b>`outerOxideTemperature`</b> - If the corrosion model is present, this keyword
    specifies the initial outer temperature (i.e. the one outside the oxide layer).
    **If not present, the initial oxide outer temperature is set to the one specified
    in  `value`.**
    - <b>`value`</b> - The initial **temperature** value (in case of oxidation model, it is
    the outer temperature of the metallic portion of the body).

    Parameters in the subdictionary for `temperatureSeries`:

    - <b>`values`</b> - Inside the dictionary `temperatureSeries`, a list of tuples of
    (time, temperature) values.
    - <b>`interpolationScheme`</b> - Inside the dictionary `temperatureSeries`, a keyword
    for selecting the type of interpolation. **By default is `linear`**.
    
________________________________________________________________________________

Usage
    To apply the `timeDependentTemperature` boundary condition to a given patch, the
    following example can be used as a template:
    
    \verbatim
        claddingOuterSurface
        {
            type            timeDependentTemperature;

            // Time dependent temperature
            temperatureSeries
            {
                // Including an external file containing the time dependent table
                // for a more readable input
                values #include    "$FOAM_CASE/constant/lists/claddingTemperature";

                // optional out-of-bounds handling
                outOfBounds         clamp;

                // optional interpolation method
                interpolationScheme linear;
            }

            // If the corrosion model is present, the initial outer temperature
            // (i.e. the one outside the oxide layer) can be specified with the
            // following keyword. Otherwise, it is set equal to "values"
            outerOxideTemperature           uniform 300;
        }
    \endverbatim

    An example of the time dependent table included in the previous example under
    the keyword `values` is given hereafter. The first column represents time units,
    while the second columns represent temperatures (K).
    
    \verbatim
        (
            (0      300.0)
            (3600   500.0)
            (7200   500.0)
            (10800  300.0)
        );
    \endverbatim

SourceFiles
    timeDependentTemperatureFvPatchScalarField.C

\mainauthor
    E. Brunetto, A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/
 
#ifndef timeDependentTemperatureFvPatchScalarField_H
#define timeDependentTemperatureFvPatchScalarField_H
 
#include "fixedTemperatureFvPatchScalarField.H"
#include "axialProfile.H"
#include "InterpolateTables.H"

#include "globalOptions.H"
#include "Table.H"

#ifdef OPENFOAMESI
    #include "Function1.H"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
        Class timeDependentTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/
 
class timeDependentTemperatureFvPatchScalarField
:
    public fixedTemperatureFvPatchScalarField
{
    // Private Data

#ifdef OPENFOAMFOUNDATION
        //- Temperature time series
        autoPtr<Function1s::Table<scalar> > tempSeries_;
#elif OPENFOAMESI
        typedef Foam::Function1Types::Table<scalar> scalarTable;

        //- Temperature time series
        autoPtr<scalarTable> tempSeries_;
#endif 
 
public:
 
    //- Runtime type information
    TypeName("timeDependentTemperature");
 
 
    // Constructors
 
        //- Construct from patch and internal field
        timeDependentTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );
 
        //- Construct from patch, internal field and dictionary
        timeDependentTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );
 
        //- Construct by mapping the given timeDependentTemperatureFvPatchScalarField
        //  onto a new patch
        timeDependentTemperatureFvPatchScalarField
        (
            const timeDependentTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );
 
 
        //- Copy constructor setting internal field reference
        timeDependentTemperatureFvPatchScalarField
        (
            const timeDependentTemperatureFvPatchScalarField&,
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
                new timeDependentTemperatureFvPatchScalarField(*this, iF)
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