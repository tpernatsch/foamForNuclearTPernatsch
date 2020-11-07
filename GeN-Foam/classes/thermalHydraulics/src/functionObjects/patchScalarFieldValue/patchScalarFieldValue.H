/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::functionObjects::patchScalarFieldValue

Group
    grpFieldFunctionObjects

Description
    Writes the values of an fvPatchScalarField 

Usage
    Example of function object specification:
    \verbatim
    patchScalarFieldValue1
    {
        type                patchScalarFieldValue;
        libs                ("libfieldFunctionObjects.so");
        ...
        writeToFile             true;
        log                     false;
        field               fieldName
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property      | Description              | Required   | Default value
        type          | type name: patchScalarFieldValue  | yes        |
        writeToFile   | write extents data to file | no       | yes
        log           | write extents data to standard output | no | yes
        field         | field name               | yes        | no
    \endtable

    Output data is written to the file \<timeDir\>/patchScalarFieldValue.dat

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    patchScalarFieldValue.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_patchScalarFieldValue_H
#define functionObjects_patchScalarFieldValue_H

#include "fvMeshFunctionObject.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class patchScalarFieldValue Declaration
\*---------------------------------------------------------------------------*/

class patchScalarFieldValue
:
    public fvMeshFunctionObject,
    public writeFile
{
public: 


protected:

    // Protected data

        //- Name of patch
        word patchName_;

        //- ID of patch
        label patchID_;

        //- Name of field
        word fieldName_;

        //- Field ptr
        const volScalarField* fieldPtr_;

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- No copy construct
        patchScalarFieldValue(const patchScalarFieldValue&) = delete;

        //- No copy assignment
        void operator=(const patchScalarFieldValue&) = delete;


public:

    //- Runtime type information
    TypeName("patchScalarFieldValue");


    // Constructors

        //- Construct from Time and dictionary
        patchScalarFieldValue
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~patchScalarFieldValue() = default;


    // Member Functions

        //- Read the field extents data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the patchScalarFieldValue
        virtual bool write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
