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
    Foam::temperatureCoupledFvPatchScalarField

Description
    Coupled boundary condition that enforces continuity for the temperature field
    between two regions with different material properties.

    \warning
    This boundary condition assumes that the same boundary condition is applied to a
    pair of `regionCoupledOFFBEAT` patches, which are  coupled using the AMI
    (Arbitrary Mesh Interface) mapping algorithm.

________________________________________________________________________________

\par Options
    The `temperatureCoupled` fvPatchField can be selected in the patch subdictionary
    inside the `boundaryField` subdictionary of the temperature field.
    
    Parameters in patch subdictionary for `temperatureCoupled`:

    - <b>`patchType`</b> - Specifies the type of underlying fvPatch and  
    must be set to `regionCoupledOFFBEAT`.
    - <b>`kappa`</b> - The name of the conductivity field. **By default, it is set
    to "k", which is the field name for conductivity used by the OFFBEAT material class.**
    - <b>`value`</b> - Specifies the initial value of the patch temperature field.

________________________________________________________________________________

Usage
    Here is a code snippet of a typical patch field dictionary:

    \verbatim
        "fuelOuter|cladInner"
        {
            type            temperatureCoupled;
            patchType       regionCoupledOFFBEAT;

            kappa           "k";

            value           $internalField;
        }
    \endverbatim

SourceFiles
    temperatureCoupledFvPatchScalarField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef temperatureCoupledFvPatchScalarField_H
#define temperatureCoupledFvPatchScalarField_H

#include "regionCoupledFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class temperatureCoupledFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class temperatureCoupledFvPatchScalarField
:
    public regionCoupledFvPatchScalarField
{
    // Private data
    

    // Private functions

protected:

    // Protected data
    
        //- Name of conductivity field
        word kappaName_;

    // Protected functions

        //- Local weight for this coupled field
        virtual tmp<scalarField> weights() const;

public:

    //- Runtime type information
    TypeName("temperatureCoupled");


    // Constructors

        //- Construct from patch and internal field
        temperatureCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        temperatureCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given temperatureCoupledFvPatchScalarField
        // onto a new patch
        temperatureCoupledFvPatchScalarField
        (
            const temperatureCoupledFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        temperatureCoupledFvPatchScalarField
        (
            const temperatureCoupledFvPatchScalarField&,
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
                new temperatureCoupledFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~temperatureCoupledFvPatchScalarField()
    {}


    // Member functions

        // Access
        
            //- Return kappa
            inline tmp<scalarField> kappa() const
            {
                return patch().lookupPatchField<volScalarField, scalar>(kappaName_);
            }

        // Evaluation functions

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
