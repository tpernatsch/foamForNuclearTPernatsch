/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::trisoGapFvPatchScalarField

Description
    Region coupled patchField for temperature fields in solid heat conduction in
    a spherical gap.

    The total gap conductance (\( h \)) is calculated as the sum of:
    - \( h_{\text{gas}} \): Gas conductance, which depends on the gap width and the
      properties of the gas in the spherical gap.
    - \( h_{\text{rad}} \): Radiative heat transfer, modeled based on the emissivity of the surfaces
      and the temperature difference.
    - \( h_{\text{contact}} \): Contact conductance, which is influenced by the surface roughness
      and deformation of the bodies.

________________________________________________________________________________

\par Options  
    The `trisoGap` fvPatchField can be selected in the patch subdictionary
    inside the `boundaryField` subdictionary of the temperature field.

    Parameters in the patch subdictionary for `trisoGap`:

    - <b>`patchType`</b> - Specifies the type of underlying fvPatch and
    must be set to `regionCoupledOFFBEAT`.
    - <b>`kappa`</b> - The name of the conductivity field. **By default, it is set
    to "k", which is the field name for conductivity used by the OFFBEAT material class.**
    - <b>`emissivity`</b> - The name of the emissivity field. **By default, it is set
    to "emissivity", which is the field name for emissivity used by the OFFBEAT material class.**
    - <b>`gapGasModel`</b> - The name of the gapGasModel class. **By default, it is set
    to "gapGas", which is the name for the OFFBEAT gap gas model class.**
    - <b>`alpha`</b> - Since this class is a gap-conductance model itself, it is not
    necessary to specify the gap conductance as it is calculated by the code.
    If `alpha` is specified, it will be used as the initial condition for the gap conductance.
    - <b>`roughness`</b> - Specifies the roughness field for the surface, used to calculate 
    the effective gap width.
    - <b>`jumpDistance`</b> - Jump distance between the surfaces in the gap. 
    If the jumpDistance is provided, it must be provided for both patches. If not provided,
    it is calculated internally by the code.
    - <b>`relax`</b> - Relaxation factor for the gap conductance. 
    **It is set to 1.0 by default.**
    - <b>`value`</b> - Specifies the initial value of the patch temperature field.

________________________________________________________________________________

Usage
    Here is a code snippet of a typical patch field dictionary:

    \verbatim
        bufferOuter
        {
            type            trisoGap;
            patchType       regionCoupledOFFBEAT;

            roughness       uniform 1e-6;
            jumpDistance    uniform 0.0;

            value           $internalField;
            relax           1;
        }

        ipycInner
        {
            type            trisoGap;
            patchType       regionCoupledOFFBEAT;
            
            roughness       uniform 1e-6;
            jumpDistance    uniform 0.0;
            
            value           $internalField;
            relax           1;
        }
    \endverbatim

SourceFiles
    trisoGapFvPatchScalarField.C

\mainauthor
    F. Xiang - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour) & XJTU (Xi'an JiaoTong
    University, China, Nuclear Thermal-hydraulic Laboratory)

\contribution
    A. Scolaro - EPFL
    E. Brunetto - EPFL

\date
    Feburary 2023


\*---------------------------------------------------------------------------*/

#ifndef fuelRodGapFvPatchScalarField_H
#define fuelRodGapFvPatchScalarField_H

#include "resistiveGapFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class trisoGapFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class trisoGapFvPatchScalarField
:
    public resistiveGapFvPatchScalarField
{
    // Private data

        //- Gap width [m]
        mutable scalarField gapWidth_;

        //- Gap conductance
        mutable scalarField hGap_;

        //- Surface roughness [m]
        scalarField R_;

        //- jumpDistance;
        scalarField jd_;

        //- Name of emissivity field
        word emissivityName_;

        //- Name of gap gas model
        word gapGasModelName_;

        //- Relaxation factor for hGap
        scalar relax_;

        //- Relaxation factor for gapWidth
        scalar relaxGapWidth_;

        //- The radii of gap patchs
        mutable scalarField r1_;
        mutable scalarField r2_;


    // Private functions

protected:

    // Protected functions

            //- Return the neighbour patchField
            inline const trisoGapFvPatchScalarField& neighbour() const
            {
                return refCast<const trisoGapFvPatchScalarField>
                (
                    regionCoupledPatch_.neighbFvPatch().lookupPatchField
                    <volScalarField, scalar>
                    (
                        internalField().name()
                    )
                );
            }

            //- Return the gap width
            tmp<scalarField> gapWidth() const;

            //- Return the interface conductance
            tmp<scalarField> hGap() const;

public:

    //- Runtime type information
    TypeName("trisoGap");


    // Constructors

        //- Construct from patch and internal field
        trisoGapFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        trisoGapFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given trisoGapFvPatchScalarField
        // onto a new patch
        trisoGapFvPatchScalarField
        (
            const trisoGapFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        trisoGapFvPatchScalarField
        (
            const trisoGapFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new trisoGapFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~trisoGapFvPatchScalarField()
    {}


    // Member functions

        // Access

            //- Return the surface roughness
            inline const scalarField& R() const
            {
                return R_;
            }

            //- Return the jumpDistance
            inline const scalarField& jd() const
            {
                return jd_;
            }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            //  Sets Updated to true
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
