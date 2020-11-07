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
    Foam::functionObjects::stopIfMaxFieldDiff

Group
    grpFieldFunctionObjects

Description
    Stops the simulation when max(field1-field2) > 0 and prints out the 
    simulation time. Both field1 and field2 must be volSclaraFields

Usage
    Example of function object specification:
    \verbatim
    stopIfMaxFieldDiff1
    {
        type                stopIfMaxFieldDiff;
        libs                ("libfieldFunctionObjects.so");
        ...
        writeToFile             yes;
        log                     yes;
        field1                  fieldName;
        field2                  fieldName;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property      | Description              | Required   | Default value
        type          | type name: stopIfMaxFieldDiff  | yes        |
        writeToFile   | write extents data to file | no       | yes
        log           | write extents data to standard output | no | yes
        field1        |                          | yes 
        field2        |                          | yes
    \endtable

Note
    For non-scalar fields, the magnitude of the field is employed and compared
    to the threshold value.

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    stopIfMaxFieldDiff.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_stopIfMaxFieldDiff_H
#define functionObjects_stopIfMaxFieldDiff_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldSelection.H"
#include "rhoThermo.H"
#include "IOobjectList.H"
#include "faceSet.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class stopIfMaxFieldDiff Declaration
\*---------------------------------------------------------------------------*/

class stopIfMaxFieldDiff
:
    public fvMeshFunctionObject
{
public: 

protected:

    // Protected data
        
        //- Name of fields
        word field1Name_;
        word field2Name_;

        //- Field ptrs
        const volScalarField* field1Ptr_;
        const volScalarField* field2Ptr_;

    // Protected Member Functions

        //- No copy construct
        stopIfMaxFieldDiff(const stopIfMaxFieldDiff&) = delete;

        //- No copy assignment
        void operator=(const stopIfMaxFieldDiff&) = delete;


public:

    //- Runtime type information
    TypeName("stopIfMaxFieldDiff");


    // Constructors

        //- Construct from Time and dictionary
        stopIfMaxFieldDiff
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~stopIfMaxFieldDiff() = default;


    // Member Functions

        //- Read the field extents data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the stopIfMaxFieldDiff
        virtual bool write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
