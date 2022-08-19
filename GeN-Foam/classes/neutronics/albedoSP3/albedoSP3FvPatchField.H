/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------

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
    Foam::albedoSP3FvPatchField

Description
    Albedo boundary condition for SP3 or diffusion calculations. Please note that the
    boundary condition needs to be set both for first and second moments in SP3.

Usage
    \verbatim
    <patchName>
    {
        type            albedoSP3;
        gamma           0.5; // defined as (1-alpha)/(1+alpha)/2, alpha being the albedo coefficient
        diffCoeffName   Dalbedo;  //not to be changed
        fluxStarAlbedo  fluxStarAlbedo; //not to be changed
        forSecondMoment false;  //true in case it is a condition for a second moment flux (for SP3 calculations)
        value           uniform 1;
    }
    \endverbatim

Author of this file (and associated .C or included .H files)
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>, EPFL (Switzerland)    

SourceFiles
    albedoSP3FvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef albedoSP3FvPatchField_H
#define albedoSP3FvPatchField_H

#include "fvPatchField.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class albedoSP3FvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class albedoSP3FvPatchField
:
    public fvPatchField<Type>
{
    // Private data

    Field<Type> gradient_;
    
    //- Coeff gamma
    scalar gamma_;
    
    //scalarField diffCoeff_;
    
    //- Name of the diffusion coefficient field
    word diffCoeffName_;

    //- Name of the other moment
    word fluxStarAlbedo_;

    //- Bool to change BC based on moment (0 or 2)
    bool forSecondMoment_;

public:

    //- Runtime type information
    TypeName("albedoSP3");


    // Constructors

    //- Construct from patch and internal field
    albedoSP3FvPatchField
    (
        const fvPatch&,
        const DimensionedField<Type, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    albedoSP3FvPatchField
    (
        const fvPatch&,
        const DimensionedField<Type, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping the given albedoSP3FvPatchField
    //  onto a new patch
    albedoSP3FvPatchField
    (
        const albedoSP3FvPatchField<Type>&,
        const fvPatch&,
        const DimensionedField<Type, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    albedoSP3FvPatchField
    (
        const albedoSP3FvPatchField<Type>&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchField<Type> > clone() const
    {
        return tmp<fvPatchField<Type> >
        (
            new albedoSP3FvPatchField<Type>(*this)
        );
    }

    //- Construct as copy setting internal field reference
    albedoSP3FvPatchField
    (
        const albedoSP3FvPatchField<Type>&,
        const DimensionedField<Type, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchField<Type> > clone
    (
        const DimensionedField<Type, volMesh>& iF
    ) const
    {
        return tmp<fvPatchField<Type> >
        (
            new albedoSP3FvPatchField<Type>(*this, iF)
        );
    }


    // Member functions

    // Return defining fields

    //- Return gradient at boundary
    virtual Field<Type>& gradient()
    {
        return gradient_;
    }

    virtual const Field<Type>& gradient() const
    {
        return gradient_;
    }
    
    virtual scalar& gamma()
    {
        return gamma_;
    }

    virtual const scalar& gamma() const
    {
        return gamma_;
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
        const fvPatchField<Type>&,
        const labelList&
    );


    // Evaluation functions

    //- Return gradient at boundary
    virtual tmp<Field<Type> > snGrad() const;

    //- Evaluate the patch field
    virtual void evaluate
    (
        const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
    );

    //- Return the matrix diagonal coefficients corresponding to the
    //  evaluation of the value of this patchField with given weights
    virtual tmp<Field<Type> > valueInternalCoeffs
    (
        const tmp<scalarField>&
    ) const;

    //- Return the matrix source coefficients corresponding to the
    //  evaluation of the value of this patchField with given weights
    virtual tmp<Field<Type> > valueBoundaryCoeffs
    (
        const tmp<scalarField>&
    ) const;

    //- Return the matrix diagonal coefficients corresponding to the
    //  evaluation of the gradient of this patchField
    virtual tmp<Field<Type> > gradientInternalCoeffs() const;

    //- Return the matrix source coefficients corresponding to the
    //  evaluation of the gradient of this patchField
    virtual tmp<Field<Type> > gradientBoundaryCoeffs() const;


    //- Write
    virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
