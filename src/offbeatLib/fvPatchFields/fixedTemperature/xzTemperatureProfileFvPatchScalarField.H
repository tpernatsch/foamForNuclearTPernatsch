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
    Foam::xzTemperatureProfileFvPatchScalarField
 
Description
    The `xzTemperatureProfile` boundary condition allows to impose a temperature
    profile constant in time but varying along x and z directions.

    \warning
    The axial locations provided to this boundary conditions are NOT normalized
    but absolute values expressed in \f$m\f$.

________________________________________________________________________________

\par Options
    The `xzTemperatureProfile` fvPatchField can be selected in the patch subdictionary
    inside the `boundaryField` subdictionary of the temperature field.

    Parameters in the patch subdictionary for `xzTemperatureProfile`:

    - <b>`axialProfileDict`</b> - Dictionary containing the information related to the
    temperature time-dependent axial profile.
    - <b>`outerOxideTemperature`</b> - If the corrosion model is present, this keyword
    specifies the initial outer temperature (i.e. the one outside the oxide layer).
    **If not present, the initial oxide outer temperature is set to the one specified
    in  `value`.**
    - <b>`value`</b> - The initial **temperature** value (in case of oxidation model, it is
    the outer temperature of the metallic portion of the body).

    Parameters in the subdictionary `axialProfileDict`:

    - <b>`xLocations`</b> - Inside the dictionary `axialProfileDict`, a list of x-location
    at witch the axial profile values are given. **The x-values are NOT
    normalized**.
    - <b>`data`</b> - Inside the dictionary `axialProfileDict`, a 2D table of the axial
    profile values in K (one row per x location; in each list one value per axial
    location).
    - <b>`zLocations`</b> - Inside the dictionary `axialProfileDict`, a list of axial
    z-location at witch the axial profile values are given. **The axial values are
    NOT normalized**.
    - <b>`axialInterpolationMethod`</b> - Inside the dictionary `axialProfileDict`, a
    keyword for selecting the type of axial interpolation (`linear` or `step`). 
    **By default is `linear`**.
    - <b>`xInterpolationMethod`</b> - Inside the dictionary `axialProfileDict`, a keyword
    for selecting the type of x-interpolation (`linear` or `step`). **By default is
    `linear`**.
    
________________________________________________________________________________

Usage
    To apply the `xzTemperatureProfile` boundary condition to a given patch, the
    following example can be used as a template:
    
    \verbatim
        claddingOuterSurface
        {
            type            xzTemperatureProfile;
            axialProfileDict
            {
                timePoints  (0 1e2 2e2 ... );
                timeInterpolationMethod  linear;

                // Note: axial locations are NOT normalized
                axialLocations ( 0 0.5 1 1.5 ... );
                axialInterpolationMethod linear;
                
                data (
                    (400 420 430 425 ... )
                    (450 470 480 475 ... )
                    (440 460 470 465 ... )
                    :
                    :
                    :
                );
            }
            value           uniform 300;

            // If the corrosion model is present, the initial outer temperature
            // (i.e. the one outside the oxide layer) can be specified with the
            // following keyword. Otherwise, it is set equal to "values"
            outerOxideTemperature           uniform 300;
        }
    \endverbatim	    

SourceFiles
    xzTemperatureProfileFvPatchScalarField.C

\todo
    current version of write(os) is not ideal since it writes the entire 
    lists for zValues, timeValues and profileData in the boundaryField section 
    of T

\mainauthor
    E. Brunetto, A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021
 
\*---------------------------------------------------------------------------*/
 
#ifndef xzTemperatureProfileFvPatchScalarField_H
#define xzTemperatureProfileFvPatchScalarField_H
 
#include "fixedTemperatureFvPatchScalarField.H"
#include "axialProfile.H"
#include "InterpolateTables.H"

#include "globalOptions.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
        Class xzTemperatureProfileFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/
 
class xzTemperatureProfileFvPatchScalarField
:
    public fixedTemperatureFvPatchScalarField
{
    // Private Data
    
        //- Tabulated z locations
        scalarField zValues_;

        //- Tabulated x locations
        scalarField xValues_;
        
        //- Tabulated Profile data
        FieldField<Field, scalar> profileData_;
    
        //- z location interpolation method
        interpolateTableBase::interpolationMethod zMethod_;
        
        //- x location interpolation method
        interpolateTableBase::interpolationMethod xMethod_;
        
        //- Interpolation table
        scalarFieldInterpolateTable table_;
        
        //- Axial profile dict
        dictionary dict_;
 
public:
 
    //- Runtime type information
    TypeName("xzTemperatureProfile");
 
 
    // Constructors
 
        //- Construct from patch and internal field
        xzTemperatureProfileFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );
 
        //- Construct from patch, internal field and dictionary
        xzTemperatureProfileFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );
 
        //- Construct by mapping the given xzTemperatureProfileFvPatchScalarField
        //  onto a new patch
        xzTemperatureProfileFvPatchScalarField
        (
            const xzTemperatureProfileFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );
 
 
        //- Copy constructor setting internal field reference
        xzTemperatureProfileFvPatchScalarField
        (
            const xzTemperatureProfileFvPatchScalarField&,
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
                new xzTemperatureProfileFvPatchScalarField(*this, iF)
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