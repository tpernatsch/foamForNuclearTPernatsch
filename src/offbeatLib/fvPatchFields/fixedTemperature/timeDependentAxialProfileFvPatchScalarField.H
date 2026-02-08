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
    Foam::timeDependentAxialProfileFvPatchScalarField
 
Description
    The `timeDependentAxialProfileT` boundary condition allows to impose a time
    dependent temperature axial profile on a patch. In case an oxide layer is present, 
    an additional thermal resistance corresponding to the oxide is considered.

    \warning
    The axial locations provided to this boundary conditions are NOT normalized
    but absolute values expressed in \f$m\f$.

________________________________________________________________________________

\par Options
    The `timeDependentAxialProfileT` fvPatchField can be selected in the patch subdictionary
    inside the `boundaryField` subdictionary of the temperature field.

    Parameters in the patch subdictionary for `timeDependentAxialProfileT`:

    - <b>`axialProfileDict`</b> - Dictionary containing the information related to the
    temperature time-dependent axial profile.
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
    - <b>`data`</b> - Inside the dictionary `axialProfileDict`, a 2D table of the axial
    profile values in K (one row per time point; in each list one value per axial
    location).
    - <b>`axialInterpolationMethod`</b> - Inside the dictionary `axialProfileDict`, a
    keyword for selecting the type of axial interpolation (`linear` or `step`). 
    **By default is `linear`**.
    - <b>`timeInterpolationMethod`</b> - Inside the dictionary `axialProfileDict`, a
    keyword for selecting the type of time interpolation (`linear` or `step`). **By
    default is `linear`**.
    
________________________________________________________________________________

Usage
    To apply the `timeDependentAxialProfileT` boundary condition to a given patch, the
    following example can be used as a template:
    
    \verbatim
        claddingOuterSurface
        {
            type            timeDependentAxialProfileT;
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
    timeDependentAxialProfileFvPatchScalarField.C

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
 
#ifndef timeDependentAxialProfileFvPatchScalarField_H
#define timeDependentAxialProfileFvPatchScalarField_H
 
#include "fixedTemperatureFvPatchScalarField.H"
#include "axialProfile.H"
#include "InterpolateTables.H"
#include "globalOptions.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
        Class timeDependentAxialProfileFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/
 
class timeDependentAxialProfileFvPatchScalarField
:
    public fixedTemperatureFvPatchScalarField
{
    // Private Data

    // Private Functions
    
        //- Calculate the axial location of each boundary face
        void calculateAxialLocations();
    
protected:
    // Protected Data
	
        //- Axial profile dict
        dictionary profileDict_;
		
		//- Reference mesh
		const fvPatch& referencePatch_;

        //- Tabulated axial locations
        scalarField zValues_;
        
        //- Tabulated time values
        const scalarField timeValues_;
        
        //- Tabulated Profile data
        FieldField<Field, scalar> profileData_;
    
        //- Axial location interpolation method
        interpolateTableBase::interpolationMethod zMethod_;
        
        //- Time interpolation method
        interpolateTableBase::interpolationMethod tMethod_;
        
        //- Interpolation table
        scalarFieldInterpolateTable table_;

		//- Z-axis
		vector axis_;

        //- Axial locations of boundary faces
        scalarField z_;
         
public:
 
    //- Runtime type information
    TypeName("timeDependentAxialProfileT");
 
 
    // Constructors
 
        //- Construct from patch and internal field
        timeDependentAxialProfileFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );
 
        //- Construct from patch, internal field and dictionary
        timeDependentAxialProfileFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );
 
        //- Construct by mapping the given timeDependentAxialProfileFvPatchScalarField
        //  onto a new patch
        timeDependentAxialProfileFvPatchScalarField
        (
            const timeDependentAxialProfileFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );
 
 
        //- Copy constructor setting internal field reference
        timeDependentAxialProfileFvPatchScalarField
        (
            const timeDependentAxialProfileFvPatchScalarField&,
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
                new timeDependentAxialProfileFvPatchScalarField(*this, iF)
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