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
    Foam::functionObjects::TBulk

Group
    grpFieldFunctionObjects

Description
    Calculates the spatial minimum and maximum extents of a field

    The extents are derived from the bound box limits after identifying the
    locations where field values exceed the user-supplied threshold value.

Usage
    Example of function object specification:
    \verbatim
    TBulk1
    {
        type                TBulk;
        libs                ("libfieldFunctionObjects.so");
        ...
        writeToFile         yes;
        log                 yes;
        temperatureField    T.liquid
        velocityField       U.liquid
        densityField        rho.liquid
        regionType          patch;
        patches             ();
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property      | Description              | Required   | Default value
        type          | type name: TBulk  | yes        |
        writeToFile   | write extents data to file | no       | yes
        log           | write extents data to standard output | no | yes
        temperatureField |                       | Yes
        velocityField |                          | Yes
        densityField  |                          | Yes
        regionType    | Either cellZone or patch | Yes 
        patches       | list of patches to process | depends on regionType
        cellZones     | list of cellZones to process | depends on regionType
    \endtable

    Output data is written to the file \<timeDir\>/TBulk.dat

Note
    For non-scalar fields, the magnitude of the field is employed and compared
    to the threshold value.

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    TBulk.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_TBulk_H
#define functionObjects_TBulk_H

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
                        Class TBulk Declaration
\*---------------------------------------------------------------------------*/

class TBulk
:
    public fvMeshFunctionObject,
    public writeFile
{
public: 

    //-
    enum regionType
    {
        patch,
        faceSet_,
        faceZone_
    };

    static const Enum<regionType> regionTypeNames_;

protected:

    // Protected data

        //- Region type
        regionType regionType_;
        
        //- Regions to assess
        wordReList regionNames_;
        
        //- List of patchIDs if regions are patches
        labelHashSet patchIDs_;

        //- List of cellFaces if regions are faceSets or faceZones
        labelList faces_;

        //- Field names
        word thermoName_;
        word alphaRhoPhiName_;
        
        //- Field ptrs
        const rhoThermo* thermoPtr_;
        const surfaceScalarField* alphaRhoPhiPtr_;

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- No copy construct
        TBulk(const TBulk&) = delete;

        //- No copy assignment
        void operator=(const TBulk&) = delete;


public:

    //- Runtime type information
    TypeName("TBulk");


    // Constructors

        //- Construct from Time and dictionary
        TBulk
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~TBulk() = default;


    // Member Functions

        //- Read the field extents data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the TBulk
        virtual bool write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
