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
    Foam::fuelRodGapFvPatchScalarField

Description
    Coupled boundary condition for modeling the presence of an evolving gap
    conductance between two bodies. The gap conductance is based primarily
    on the Frapcon 4.0 manual, making it most suitable for standard fuel rods
    (fuel pellet + cladding).

    The total gap conductance \f$h\f$ is calculated as the sum of:
    
    - \f$h_{\text{gas}}\f$: Gas conductance, which depends on the gap width, surface
    roughness, jump distance and the properties of the gas mixture in the gap.
    - \f$h_{\text{rad}}\f$: Radiative heat transfer, modeled based on the emissivity of 
    the surfaces and the temperature difference.
    - \f$h_{\text{contact}}\f$: Contact conductance, which depends on the surface roughness and
    contact pressure

    \warning  
    This boundary condition assumes that the same boundary condition is applied to a
    pair of `regionCoupledOFFBEAT` patches, which are coupled using the AMI
    (Arbitrary Mesh Interface) mapping algorithm.

    \warning  
    This boundary condition requires the presence of a `gapGasModel` in the
    simulation.

________________________________________________________________________________

\par Options
    The `fuelRodGap` fvPatchField can be selected in the patch subdictionary
    inside the `boundaryField` subdictionary of the temperature field.

    Parameters in the patch subdictionary for `fuelRodGap`:

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
    - <b>`relax`</b> - Relaxation factor for the gap conductance. **It is set to 1.0 by default,
    but values as low as 0.1 are often very useful to aid (and sometimes accelerate) convergence.**
    - <b>`value`</b> - Specifies the initial value of the patch temperature field.

________________________________________________________________________________

Usage
    Here is a code snippet of a typical patch field dictionary:

    \verbatim
        fuelOuter
        {
            type            fuelRodGap;
            patchType       regionCoupledOFFBEAT;

            roughness       uniform 1e-6;
            relax           0.1;

            value           $internalField;
        } 
    \endverbatim

SourceFiles
    fuelRodGapFvPatchScalarField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef fuelRodGapFvPatchScalarField_H
#define fuelRodGapFvPatchScalarField_H

#include "resistiveGapFvPatchScalarField.H"
#include "userParameters.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class fuelRodGapFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fuelRodGapFvPatchScalarField
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
        
        //- Name of emissivity field
        word emissivityName_;

        //- Name of gap gas model
        word gapGasModelName_;    

        //- Relaxation factor for hGap
        scalar relax_;    

    // Private functions
    

    // Scaling parameters for UQ and SA
        
        //- Group name for parameter definition
        declareNameWithVirtual(fuelRodGapFvPatchScalarField, group_, group);
    
        //- Densification strain scaling factor [-]
        declareParameterWithVirtual(fuelRodGapFvPatchScalarField, 
            F_hGas_, F_hGas);
        declareParameterWithVirtual(fuelRodGapFvPatchScalarField, 
            F_hRad_, F_hRad);
        declareParameterWithVirtual(fuelRodGapFvPatchScalarField, 
            F_hContact_, F_hContact);

        //- Densification strain offset [-]
        declareParameterWithVirtual(fuelRodGapFvPatchScalarField, 
            delta_hGas_, delta_hGas);
        declareParameterWithVirtual(fuelRodGapFvPatchScalarField, 
            delta_hRad_, delta_hRad);
        declareParameterWithVirtual(fuelRodGapFvPatchScalarField, 
            delta_hContact_, delta_hContact);

protected:

    // Protected functions
    
            //- Return the neighbour patchField
            inline const fuelRodGapFvPatchScalarField& neighbour() const
            {
                return refCast<const fuelRodGapFvPatchScalarField>
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
    TypeName("fuelRodGap");


    // Constructors

        //- Construct from patch and internal field
        fuelRodGapFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fuelRodGapFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fuelRodGapFvPatchScalarField
        // onto a new patch
        fuelRodGapFvPatchScalarField
        (
            const fuelRodGapFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        fuelRodGapFvPatchScalarField
        (
            const fuelRodGapFvPatchScalarField&,
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
                new fuelRodGapFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~fuelRodGapFvPatchScalarField()
    {}


    // Member functions

        // Access
        
            //- Return the surface roughness
            inline const scalarField& R() const
            {
                return R_;
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
