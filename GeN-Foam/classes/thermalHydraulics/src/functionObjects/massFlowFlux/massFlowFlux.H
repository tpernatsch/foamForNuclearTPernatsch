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
    Foam::functionObjects::massFlowFlux

Group
    grpFieldFunctionObjects

Description
    Calculates the total mass flux through a set of faceSets, faceZones or
    patches.

Usage
    Example of function object specification:
    \verbatim
    massFlowFlux1
    {
        type                massFlowFlux;
        libs                ("libfieldFunctionObjects.so");
        ...
        writeToFile         yes;
        log                 yes;
        thermo              thermophysicalProperties.liquid;
        U                   U.liquid
        alpha               alpha.liquid
        regionType          patch;
        patches             (outlet);
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property      | Description              | Required   | Default value
        type          | type name: massFlowFlux      | yes        |
        writeToFile   | write extents data to file | no       | yes
        log           | write extents data to standard output | no | yes
        thermo        | thermophysicalProperties | yes
        U             | velocity field name      | yes
        alpha         | phase fraction field name| no | defaults to field of 1
        regionType    | Either cellZone or patch | Yes 
        patches       | list of patches to process | depends on regionType
        cellZones     | list of cellZones to process | depends on regionType
    \endtable

    Output data is written to the file \<timeDir\>/massFlowFlux.dat

Note
    For non-scalar fields, the magnitude of the field is employed and compared
    to the threshold value.

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    massFlowFlux.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_massFlowFlux_H
#define functionObjects_massFlowFlux_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldSelection.H"
#include "rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class massFlowFlux Declaration
\*---------------------------------------------------------------------------*/

class massFlowFlux
:
    public fvMeshFunctionObject,
    public writeFile
{
public: 

    //-
    enum regionType
    {
        patch,
        faceSet,
        faceZone
    };

    static const Enum<regionType> regionTypeNames_;

protected:

    // Protected data

        //- Patches to assess
        wordReList patchNames_;
        labelHashSet patchIDs_;
        word compoundPatchesName_;

        //- Region type
        regionType regionType_;

        //- Field names
        word thermoName_;
        word alphaName_;
        word UName_;

        //- Field ptrs
        const rhoThermo* thermoPtr_;
        const volScalarField* alphaPtr_;
        const volVectorField* UPtr_;

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- No copy construct
        massFlowFlux(const massFlowFlux&) = delete;

        //- No copy assignment
        void operator=(const massFlowFlux&) = delete;


public:

    //- Runtime type information
    TypeName("massFlowFlux");


    // Constructors

        //- Construct from Time and dictionary
        massFlowFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~massFlowFlux() = default;


    // Member Functions

        //- Read the field extents data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the massFlowFlux
        virtual bool write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
