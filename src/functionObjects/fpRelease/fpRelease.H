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
    fpRelease

Description
    FunctionObject that writes fission product release fractions vs time.

    Each species solved by fissionProductsDiffusionSolver registers a
    uniformDimensionedScalarField named "fpRelease_<species>" in the mesh
    registry.  This function object auto-detects all such fields and writes
    a single CSV file:

    \verbatim
    postProcessing/fpRelease/<startTime>/fpRelease.dat
    \endverbatim

Usage
    In controlDict:
    \verbatim
    functions
    {
        writeFpRelease
        {
            libs ("libFFNFunctionObjects.so");
            type fpRelease;
        }
    }
    \endverbatim

SourceFiles
    fpRelease.C

\*---------------------------------------------------------------------------*/

#ifndef fpRelease_H
#define fpRelease_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fpRelease Declaration
\*---------------------------------------------------------------------------*/

class fpRelease
:
    public functionObject
{
    // Private data

        const word name_;
        const Time& time_;
        word meshName_;
        autoPtr<OFstream> filePtr_;

        //- Species names, populated on the first execute() call
        wordList speciesNames_;


    // Private Member Functions

        //- Detect species from registry and write one line to file
        bool writeData();

        fpRelease(const fpRelease&);
        void operator=(const fpRelease&);


public:

    TypeName("fpRelease");


    // Constructors

        fpRelease
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        virtual bool start();
        virtual bool execute();
        virtual bool read(const dictionary&);
        virtual bool write();
        virtual void updateMesh(const mapPolyMesh&) {}
        virtual void movePoints(const pointField&) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
