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
    Foam::functionObjects::pressureDrop

Group
    grpFieldFunctionObjects

Description
    Calculates the spatial minimum and maximum extents of a field

    The extents are derived from the bound box limits after identifying the
    locations where field values exceed the user-supplied threshold value.

Usage
    Example of function object specification:
    \verbatim
    pressureDrop1
    {
        type                pressureDrop;
        libs                ("libfieldFunctionObjects.so");
        ...
        writeToFile             yes;
        log                     yes;
        upstreamRegionType      patch/faceSet;
        upstreamRegion          inlet;
        downstreamRegionType    faceSet;
        downstreamRegion        outlet;
        pressure                p_rgh
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property      | Description              | Required   | Default value
        type          | type name: pressureDrop  | yes        |
        writeToFile   | write extents data to file | no       | yes
        log           | write extents data to standard output | no | yes
        regionType    | Either cellZone or patch | Yes 
        patches       | list of patches to process | depends on regionType
        cellZones     | list of cellZones to process | depends on regionType
    \endtable

    Output data is written to the file \<timeDir\>/pressureDrop.dat

Note
    For non-scalar fields, the magnitude of the field is employed and compared
    to the threshold value.

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    pressureDrop.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_pressureDrop_H
#define functionObjects_pressureDrop_H

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
                        Class pressureDrop Declaration
\*---------------------------------------------------------------------------*/

class pressureDrop
:
    public fvMeshFunctionObject,
    public writeFile
{
public: 

    //-
    enum regionType
    {
        patch,
        faceSet_
    };

    static const Enum<regionType> regionTypeNames_;

protected:

    // Protected data

        //- Upstream region type
        regionType region1Type_;
        //- Downstream region type
        regionType region2Type_;
        
        //- Regions to assess
        word region1Name_;
        word region2Name_;
        
        //- 
        label patchID1_;
        label patchID2_;
        labelList faces1_;
        labelList faces2_;

        //- Surface areas of region 1,2
        scalar S1_;
        scalar S2_;
        
        //- Name of pressure field
        word pName_;

        //- Field ptr
        const volScalarField* pPtr_;

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- No copy construct
        pressureDrop(const pressureDrop&) = delete;

        //- No copy assignment
        void operator=(const pressureDrop&) = delete;


public:

    //- Runtime type information
    TypeName("pressureDrop");


    // Constructors

        //- Construct from Time and dictionary
        pressureDrop
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~pressureDrop() = default;


    // Member Functions

        //- Read the field extents data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the pressureDrop
        virtual bool write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
