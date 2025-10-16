/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::fixedMassFlowRateFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    Exactly the same as OpenFOAM's flowRateInletVelocity, but renormalizes the
    velocity to account for the porosity at the patch.

    Author of this file (and associated .C or included .H files):
    Giovanni Nervi <giovanni.nervi@epfl.ch>, EPFL (Switzerland)


Usage
    \table
        Property     | Description             | Required    | Default value
        massFlowRate | mass flow rate [kg/s]   | no          |
        volumetricFlowRate | volumetric flow rate [m3/s]| no |
        rho          | density field name      | no          | thermo:rho
        rhoInlet     | inlet density           | no          |
        extrapolateProfile | Extrapolate velocity profile | no | false
    \endtable

    Example of the boundary condition specification for a volumetric flow rate:
    \verbatim
    <patchName>
    {
        type                fixedMassFlowRate;
        volumetricFlowRate  0.2;
        extrapolateProfile  yes;
        value               uniform (0 0 0);
    }
    \endverbatim

    Example of the boundary condition specification for a mass flow rate:
    \verbatim
    <patchName>
    {
        type                fixedMassFlowRate;
        massFlowRate        0.2;
        extrapolateProfile  yes;
        rho                 rho;
        rhoInlet            1.0;
        value               uniform (0 0 0);
    }
    \endverbatim

    The \c flowRate entry is a \c Function1 of time, see Foam::Function1Types.

Note
    - \c rhoInlet is required for the case of a mass flow rate, where the
      density field is not available at start-up
    - The value is positive into the domain (as an inlet)
    - May not work correctly for transonic inlets
    - Strange behaviour with potentialFoam since the U equation is not solved

SourceFiles
    fixedMassFlowRateFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedMassFlowRateFvPatchVectorField_H
#define fixedMassFlowRateFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class fixedMassFlowRateFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class fixedMassFlowRateFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Inlet integral flow rate
        autoPtr<Function1<scalar>> flowRate_;

        //- Name of the density field used to normalize the mass flux
        word rhoName_;

        //- Rho initialisation value (for start; if value not supplied)
        scalar rhoInlet_;

        //- Is volumetric?
        bool volumetric_;

        //- Set true to extrapolate the velocity profile from the interior
        Switch extrapolateProfile_;


    // Private Member Functions

        //- Update the patch values given the appropriate density type and value
        template<class RhoType>
        void updateValues(const RhoType& rho);


public:

   //- Runtime type information
   TypeName("fixedMassFlowRate");


   // Constructors

        //- Construct from patch and internal field
        fixedMassFlowRateFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedMassFlowRateFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fixedMassFlowRateFvPatchVectorField
        //  onto a new patch
        fixedMassFlowRateFvPatchVectorField
        (
            const fixedMassFlowRateFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedMassFlowRateFvPatchVectorField
        (
            const fixedMassFlowRateFvPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        fixedMassFlowRateFvPatchVectorField
        (
            const fixedMassFlowRateFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return fvPatchField<vector>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return fvPatchField<vector>::Clone(*this, iF);
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
