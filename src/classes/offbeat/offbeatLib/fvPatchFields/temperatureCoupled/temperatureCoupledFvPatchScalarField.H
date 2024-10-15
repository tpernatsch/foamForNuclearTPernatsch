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
    Coupled boundary condition that forces continuity (for the temperature field) 
    between two bodies with different material properties.

Usage
    \verbatim
    fuelOuter
    {
        type            temperatureCoupled;
        patchType       regionCoupledOFFBEAT;
        coupled         true;
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
