/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::hydrogenTransportFvPatchScalarField

Description
    Boundary condition to use with hydrogen transport. 

    TO USE ON THE Css FIELD!. Use zeroGradient for Ctot/Cpp fields.
    
    This boundary condition applies, by default, a 0 current boundary condition.

    If the HPickUp keyword is set to true, it also applies an inward H flux corresponding
    to the hydrogen pick up due to the waterside corrosion of the cladding.
    In this case, if the volFactor keyword is set to true, the user must enter 
    the clad inner and outer radius (m).
    This keyword allows the user to take into account the reduced volume of clad 
    vs the growing oxide layer and leads to a more precise H flux. 

Usage

    boundaryField
    {
    ...
        cladOuter
        {
            type hydrogenTransport;

            Q* 25500; //Optional, heat of transport. Must corresponds to the value used in the transport solver.           
            
            HPickUp true; //Optional. Default is false. Activate Hydrogen pick up and applies an ingress flux.
            
            volFactor true; //Optional. Default is false.
            rinner 0.004565; //Required if volFactor is set to true. Clad inner radius (m)
            router 0.005125; //Required if volFactor is set to true. Clad outer radius (m)      

            fH 0.15;// Optional. Default is 0.15. Hydrogen pick up fraction.
                    // (NB. Current HNGD implementation is for zyr4)
                    // Zircaloyr4 --> 0.15 -- 0.2
                    // ZIRLO --> 0.25
                    // Optimized ZIRLO --> 0.25
                    // M5 --> 0.15
        }
    }
SourceFiles
    hydrogenTransportFvPatchScalarField.C

\mainauthor
    M. Reymond - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    I. Clifford - PSI
    E. Brunnetto, A. Scolaro, C. Fiorina - EPFL

\date 
    April 2023

\*---------------------------------------------------------------------------*/

#ifndef hydrogenTransportFvPatchScalarField_H
#define hydrogenTransportFvPatchScalarField_H

#include "fvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class hydrogenTransportFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class hydrogenTransportFvPatchScalarField
:
    public fvPatchScalarField
{
    // Private Data
            
protected:

    // Protected Data
    
        //- Constant hydrogen ingress rate
        scalarField ingressRate_;
        
        //- Soret and/or stress-driven diffusion fluxes
        scalarField phi_;
        
        //- DH*deltaCoeffs
        scalarField deltaCoeffsD_;

public:

    //- Runtime type information
    TypeName("hydrogenTransport");


    // Constructors

        //- Construct from patch and internal field
        hydrogenTransportFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        hydrogenTransportFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  zeroCurrentHDiffusion onto a new patch
        hydrogenTransportFvPatchScalarField
        (
            const hydrogenTransportFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        hydrogenTransportFvPatchScalarField
        (
            const hydrogenTransportFvPatchScalarField&,
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
                new hydrogenTransportFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

    // Member functions

        // Access
        
        //- Return true if this patch field fixes a value.
        //  Needed to check if a level has to be specified while solving
        //  Poissons equations.
        virtual bool fixesValue() const
        {
            return true;
        }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatchScalarField&, const labelList&);


        //- Update
        virtual void updateCoeffs();

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        //- Return the matrix diagonal coefficients corresponding to the evaluation of the value of this patchField
        virtual tmp<Field<scalar> > valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix source coefficients corresponding to the evaluation of the value of this patchField,
        virtual tmp<Field<scalar> > valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix diagonal coefficients corresponding to the evaluation of the gradient of this patchField
        virtual tmp<Field<scalar> > gradientInternalCoeffs() const;

        //- Return the matrix source coefficients corresponding to the evaluation of the gradient of this patchField
        virtual tmp<Field<scalar> > gradientBoundaryCoeffs() const;
        
        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
