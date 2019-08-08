/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::pUDropFvPatchScalarField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a jump condition for a scalar field, 
    using the cyclic condition as a base. It was specifically developed
    for pressure drops at orifices/openings.

    The pressure jump/drop at the patch is computed as 0.5*K*rho*mag(U)**2,
    when U is the velocity of the fluid through the patch and K is a user
    provided, constant pressure loss coefficient. The direction of the pressure jump
    depends on the flow direction through the baffle.

    By Stefan Radman, stefan.radman@epfl.ch

Usage
    Example of the boundary condition specification:
    \verbatim
    <masterPatchName>
    {
        type            pUDrop;
        patchType       cyclic;
        K               1.17; // loss coefficient 
        UName           realU; // name of the physical velocity field in the solver
        uRelax          0; // if true, underrelaxes the pressure drop with an internally computed factor (use if K is too large)
        pInfo           0; // print total flow through patch and pressure drop
        value           uniform 0; // starting condition
    }
    <slavePatchName>
    {
        type            pUDrop;
        patchType       cyclic;
        value           uniform 0; // starting condition
    }
    \endverbatim

Note
     The underlying \c patchType should be set to \c cyclic

See also
    Foam::jumpCyclicFvPatchField

SourceFiles
    pUDropFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef pUDropFvPatchScalarField_H
#define pUDropFvPatchScalarField_H

#include "jumpCyclicFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class pUDropFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class pUDropFvPatchScalarField
:
    public jumpCyclicFvPatchField<Foam::scalar>
{

protected:

    // Protected data

        //- resistance coefficient so that Dp = K_*rho*(U & U)
        scalar K_;

        //- name of the velocity field
        word UName_; 

        //- crude old time step based under-relaxation
        bool uRelax_;

        //- print info on pDrop
        bool pInfo_;

public:

    //- Runtime type information
    TypeName("pUDrop");

    // Constructors

        //- Construct from patch and internal field
        pUDropFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pUDropFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given pUDropFvPatchScalarField onto a
        //  new patch
        pUDropFvPatchScalarField
        (
            const pUDropFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pUDropFvPatchScalarField
        (
            const pUDropFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new pUDropFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pUDropFvPatchScalarField
        (
            const pUDropFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new pUDropFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the "jump" across the patch
            virtual tmp<Field<scalar>> jump() const;


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
