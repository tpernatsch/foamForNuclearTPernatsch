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
    Foam::timeDependentTemperatureFvPatchScalarField

Description
    Boundary condition to impose a time-dependent fixed temperature.

    \htmlinclude fvPatchField_timeDependentTemperatureFvPatchScalarField.html

SourceFiles
    timeDependentTemperatureFvPatchScalarField.C

\mainauthor
    E. Brunetto, A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/
 
#ifndef timeDependentTemperatureFvPatchScalarField_H
#define timeDependentTemperatureFvPatchScalarField_H
 
#include "fixedTemperatureFvPatchScalarField.H"
#include "axialProfile.H"
#include "InterpolateTables.H"

#include "globalOptions.H"
#include "Table.H"

#ifdef OPENFOAMESI
    #include "Function1.H"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
        Class timeDependentTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/
 
class timeDependentTemperatureFvPatchScalarField
:
    public fixedTemperatureFvPatchScalarField
{
    // Private Data

#ifdef OPENFOAMFOUNDATION
        //- Temperature time series
        autoPtr<Function1s::Table<scalar> > tempSeries_;
#elif OPENFOAMESI
        typedef Foam::Function1Types::Table<scalar> Table;

        //- Temperature time series
        autoPtr<Table> tempSeries_;
#endif 
 
public:
 
    //- Runtime type information
    TypeName("timeDependentTemperature");
 
 
    // Constructors
 
        //- Construct from patch and internal field
        timeDependentTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );
 
        //- Construct from patch, internal field and dictionary
        timeDependentTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );
 
        //- Construct by mapping the given timeDependentTemperatureFvPatchScalarField
        //  onto a new patch
        timeDependentTemperatureFvPatchScalarField
        (
            const timeDependentTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );
 
 
        //- Copy constructor setting internal field reference
        timeDependentTemperatureFvPatchScalarField
        (
            const timeDependentTemperatureFvPatchScalarField&,
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
                new timeDependentTemperatureFvPatchScalarField(*this, iF)
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