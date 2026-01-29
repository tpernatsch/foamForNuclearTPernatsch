/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::andersonMixingScheme

Description
    Solution acceleration using restarted Anderson mixing.

Usage
    Here is a code snippet of the dictionary to be used for activating
    `andersonMixing` acceleration:

    \verbatim
        accelerationScheme
        {
            type        andersonMixing;

            // Order of acceleration (number of snapshots minus 1)
            order       5;

            // Blending coefficient (default 1)
            alpha       1;

            // Increase the main diagonal of the projection matrix by this
            // factor to ensure diagonal dominance (default 1e-4)
            diagonalFactor 1e-4;
        }
    \endverbatim

SourceFiles
    andersonMixingScheme.C
    andersonMixingSchemes.H
    andersonMixingSchemes.C

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date
    January 2025

\*---------------------------------------------------------------------------*/

#ifndef andersonMixingScheme_H
#define andersonMixingScheme_H

#include "accelerationScheme.H"
#include "FieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class andersonMixingScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class andersonMixingScheme
:
    public accelerationScheme<Type>
{
    typedef GeometricField<Type, fvPatchField, volMesh> volTypeField;
    
    // Private data

        //- Multiply the diagonal of the projection matrix by this factor
        //- to ensure diagonal dominance
        const scalar diagonalFactor_;

        //- Blending factor
        const scalar alpha_;

        //- Acceleration order (number of snapshots minus 1)
        const label order_;

        //- Solution snapshots for acceleration
        PtrList<PtrList<volTypeField> > x_;

        //- Number of fields added
        label count_;

    // Private Member Functions

        //- Project a field onto another
        scalar projectFields(const Field<Type>& f1, const Field<Type>& f2) const;

        //- Disallow default bitwise copy construct
        andersonMixingScheme(const andersonMixingScheme&);

        //- Disallow default bitwise assignment
        void operator=(const andersonMixingScheme&);


public:

    //- Runtime type information
    TypeName("andersonMixing");

    // Constructors

        //- Construct from dictionary
        andersonMixingScheme(label size, const dictionary& dict);


    //- Destructor
    virtual ~andersonMixingScheme();


    // Member Functions

        //- Reset
        virtual void reset();

        //- Perform acceleration
        virtual void accelerate();

        //- Minumum number of iterations required for acceleration
        virtual int minIter() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "andersonMixingScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
