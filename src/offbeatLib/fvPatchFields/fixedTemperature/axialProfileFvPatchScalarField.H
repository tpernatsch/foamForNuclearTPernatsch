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
    Foam::axialProfileFvPatchScalarField
 
Description
    The `axialProfileT` boundary condition allows to impose a temperature axial
    profile (constant in time) on a patch. In case an oxide layer is present, 
    an additional thermal resistance corresponding to the oxide is considered.

    \warning
    The axial locations provided to this boundary conditions are NOT normalized
    but absolute values expressed in \f$m\f$.

________________________________________________________________________________

\par Options
    The `axialProfileT` fvPatchField can be selected in the patch subdictionary
    inside the `boundaryField` subdictionary of the temperature field.

    Parameters in the patch subdictionary for `axialProfileT`:

    - <b>`axialProfileDict`</b> - Dictionary containing the information related to the
    temperature axial-profile.
    - <b>`outerOxideTemperature`</b> - If the corrosion model is present, this keyword
    specifies the initial outer temperature (i.e. the one outside the oxide layer).
    **If not present, the initial oxide outer temperature is set to the one specified
    in  `value`.**
    - <b>`value`</b> - The initial **temperature** value (in case of oxidation model, it is
    the outer temperature of the metallic portion of the body).

    Parameters in the subdictionary `axialProfileDict`:
    
    - <b>`axialLocations`</b> - Inside the dictionary `axialProfileDict`, a list of axial
    location at witch the axial profile values are given. **The axial values are NOT
    normalized**.
    - <b>`data`</b> - Inside the dictionary `axialProfileDict`, a list of the axial profile
    values in K (one value per axial location given in the previous list).
    - <b>`axialInterpolationMethod`</b> - Inside the dictionary `axialProfileDict`, a
    keyword for selecting the type of axial interpolation (`linear` or `step`). 
    **By default is `linear`**.
    
________________________________________________________________________________

Usage
    To apply the `axialProfileT` boundary condition to a given patch, the following
    example can be used as a template:
    
    \verbatim
        claddingOuterSurface
        {
            type            xzTemperatureProfile;
            axialProfileDict
            {
                xLocations  ( -0.005 0.005 );
                xInterpolationMethod linear;

                zLocations ( 0 4 );
                zInterpolationMethod linear;

                data
                (
                    //  T @ z=0   T @ z=4
                    (   300       400   ) // T @ x = -0.005
                    (   300       500   ) // T @ x =  0.005
                );
            }

            value           uniform 300;

            // If the corrosion model is present, the initial outer temperature
            // (i.e. the one outside the oxide layer) can be specified with the
            // following keyword. Otherwise, it is set equal to "value"
            outerOxideTemperature           uniform 300;
        }
    \endverbatim	    

SourceFiles
    axialProfileFvPatchScalarField.C

\mainauthor
    E. Brunetto, A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021
 
\*---------------------------------------------------------------------------*/
 
#ifndef axialProfileFvPatchScalarField_H
#define axialProfileFvPatchScalarField_H
 
#include "fixedTemperatureFvPatchScalarField.H"
#include "axialProfile.H"
#include "InterpolateTables.H"

#include "globalOptions.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
        Class axialProfileFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/
 
class axialProfileFvPatchScalarField
:
    public fixedTemperatureFvPatchScalarField
{
    // Private Data
    
        //- Tabulated axial locations
        scalarField zValues_;
        
        //- Tabulated Profile data
        scalarField profileData_;
    
        //- Axial location interpolation method
        interpolateTableBase::interpolationMethod zMethod_;
        
        //- Interpolation table
        scalarInterpolateTable table_;

        //- Reference axial locations
        scalarField zRef_;
        
        //- Axial profile dict
        dictionary dict_;
 
public:
 
    //- Runtime type information
    TypeName("axialProfileT");
 
 
    // Constructors
 
        //- Construct from patch and internal field
        axialProfileFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );
 
        //- Construct from patch, internal field and dictionary
        axialProfileFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );
 
        //- Construct by mapping the given axialProfileFvPatchScalarField
        //  onto a new patch
        axialProfileFvPatchScalarField
        (
            const axialProfileFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );
 
 
        //- Copy constructor setting internal field reference
        axialProfileFvPatchScalarField
        (
            const axialProfileFvPatchScalarField&,
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
                new axialProfileFvPatchScalarField(*this, iF)
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