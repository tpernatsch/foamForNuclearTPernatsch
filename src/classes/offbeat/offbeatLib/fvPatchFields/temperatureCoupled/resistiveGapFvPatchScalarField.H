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
    Foam::resistiveGapFvPatchScalarField

Description
    Coupled boundary condition for modeling the presence of a fixed gap (or contact)
    resistance between two bodies.

Usage
    \verbatim
    fuelOuter
    {
        type            resistiveGap;
        patchType       regionCoupledOFFBEAT;
        coupled         true;
        alpha           uniform 5000;
        value           $internalField;
    } 
    \endverbatim

SourceFiles
    resistiveGapFvPatchScalarField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef resistiveGapFvPatchScalarField_H
#define resistiveGapFvPatchScalarField_H

#include "temperatureCoupledFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class resistiveGapFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class resistiveGapFvPatchScalarField
:
    public temperatureCoupledFvPatchScalarField
{
    // Private data

    // Private functions

protected:

    // Protected functions
    
        //- Thermal resistance across interface [W/m2.K]
        scalarField alpha_gap_;

        //- Local weight for this coupled field
        virtual tmp<scalarField> weights() const;

public:

    //- Runtime type information
    TypeName("resistiveGap");


    // Constructors

        //- Construct from patch and internal field
        resistiveGapFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        resistiveGapFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct from patch, internal field, dictionary and bool deciding 
        //- if gap resistance must be read from dictionary
        resistiveGapFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool
        );

        //- Construct by mapping given resistiveGapFvPatchScalarField
        // onto a new patch
        resistiveGapFvPatchScalarField
        (
            const resistiveGapFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        resistiveGapFvPatchScalarField
        (
            const resistiveGapFvPatchScalarField&,
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
                new resistiveGapFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~resistiveGapFvPatchScalarField()
    {}


    // Member functions

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
            
        // Access
        
            //- Return the interface heat transfer coefficient
            inline const scalarField& alpha() const
            {
                return alpha_gap_;
            }

            //- Return the interface heat transfer coefficient
            inline scalarField& alpha()
            {
                return alpha_gap_;
            }

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
