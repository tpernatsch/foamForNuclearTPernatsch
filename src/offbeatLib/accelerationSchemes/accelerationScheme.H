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
    Foam::accelerationScheme

Description
    Base class for solution acceleration/extrapolation based on a series of
    solution snapshots.
    The base class performs no acceleration.

Usage
    Here is a code snippet of the dictionary to be used for disabling
    acceleration:

    \verbatim
        accelerationScheme
        {
            type        none;
        }
    \endverbatim

SourceFiles
    accelerationScheme.C
    accelerationSchemes.H
    accelerationSchemes.C

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date
    January 2025

\*---------------------------------------------------------------------------*/

#ifndef accelerationScheme_H
#define accelerationScheme_H

#include "volFields.H"
#include "FieldField.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class accelerationScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class accelerationScheme
{
    typedef GeometricField<Type, fvPatchField, volMesh> volTypeField;
    typedef accelerationScheme<Type> typeAccelerationScheme;

protected:

    // Protected data
    
        //- Reference to volFields to be accelerated
        UPtrList<volTypeField> psi_;

        //- Number of fields added
        label nFields_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        accelerationScheme(const accelerationScheme&);

        //- Disallow default bitwise assignment
        void operator=(const accelerationScheme&);


public:

    //- Runtime type information
    TypeName("none");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            typeAccelerationScheme,
            dictionary,
            (
                label size,
                const dictionary& dict
            ),
            (size, dict)
        );

    // Selectors

        //- Return a reference to the selected thermal-fluid region
        static autoPtr<accelerationScheme<Type> > New
        (
            label size,
            const dictionary& dict
        );

    // Constructors

        //- Construct from dictionary
        accelerationScheme(label size, const dictionary& dict);


    //- Destructor
    virtual ~accelerationScheme();


    // Member Functions

        //- Add a new volField for acceleration
        void addField(volTypeField& psi);

        //- Reset
        virtual void reset();

        //- Perform acceleration
        virtual void accelerate();

        //- Minumum number of iterations required for acceleration
        virtual int minIter() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Macros for instantiating templates
#define CAT_(a, b) a ## b
#define CAT(a, b) CAT_(a, b)
#define CAPITALIZE(name) CAPITALIZE_##name
#define CAPITALIZE_accelerationScheme AccelerationScheme

#define declareAccelerationSchemeType(schemeType, type)                     \
    typedef schemeType<type> CAT(type,CAPITALIZE(schemeType))

//- Create the base acceleration scheme for a single type
#define makeAccelerationSchemeType(type)                                    \
    defineNamedTemplateTypeNameAndDebug                                     \
    (                                                                       \
        type##AccelerationScheme,                                           \
        0                                                                   \
    );                                                                      \
    defineTemplateRunTimeSelectionTable                                     \
    (                                                                       \
        type##AccelerationScheme,                                           \
        dictionary                                                          \
    );                                                                      \
    addToRunTimeSelectionTable                                              \
    (                                                                       \
        type##AccelerationScheme,                                           \
        type##AccelerationScheme,                                           \
        dictionary                                                          \
    )

//- Add a new runtime selectable acceleration scheme for a single type
#define addAccelerationSchemeType(schemeType, type)                         \
    defineNamedTemplateTypeNameAndDebug                                     \
    (                                                                       \
        schemeType<type>,                                                   \
        0                                                                   \
    );                                                                      \
    addTemplateToRunTimeSelectionTable                                      \
    (                                                                       \
        type##AccelerationScheme,                                           \
        schemeType,                                                         \
        type,                                                               \
        dictionary                                                          \
    )

//- Declare a new runtime selectable acceleration scheme
#define declareAccelerationSchemes(schemeType)                              \
    declareAccelerationSchemeType(schemeType, scalar);                      \
    declareAccelerationSchemeType(schemeType, vector);                      \
    declareAccelerationSchemeType(schemeType, tensor);                      \
    declareAccelerationSchemeType(schemeType, symmTensor)

//- Create the base acceleration scheme for all types
#define makeAccelerationSchemes()                                           \
    makeAccelerationSchemeType(scalar);                                     \
    makeAccelerationSchemeType(vector);                                     \
    makeAccelerationSchemeType(tensor);                                     \
    makeAccelerationSchemeType(symmTensor)

//- Add a new runtime selectable acceleration scheme for all types
#define addAccelerationSchemes(schemeType)                                  \
    addAccelerationSchemeType(schemeType, scalar);                          \
    addAccelerationSchemeType(schemeType, vector);                          \
    addAccelerationSchemeType(schemeType, tensor);                          \
    addAccelerationSchemeType(schemeType, symmTensor)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "accelerationScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
