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

\*----------------------------------------------------------------------------*/

#include "rodPressure.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "OSspecific.H"
#include "gapGasModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rodPressure, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        rodPressure,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::rodPressure::writeData()
{
    // Lookup mesh
    const fvMesh& mesh = time_.lookupObject<fvMesh>("region0");

    const word objName("gapGas");

    if(mesh.foundObject<gapGasModel>(objName))
    {
        const scalar& rodPressure = 
            mesh.lookupObject<gapGasModel>(objName).p();

        // Write to file
        if (Pstream::master())
        {
            historyFilePtr_() << mesh.time().value() << " " << rodPressure << endl;
        }
    }
    else
    {
        WarningIn("Foam::rodPressure::writeData()") << nl  
        << "    Cannot find gapGasModel object for postprocessing !" 
        <<  nl << endl;

        // Write to file
        if (Pstream::master())
        {
            historyFilePtr_() << mesh.time().value() << " " << "-1" << endl;
        }
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rodPressure::rodPressure
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    historyFilePtr_()
{
    Info<< "Creating " << this->name() << " function object" << endl;

    const fvMesh& mesh = time_.lookupObject<fvMesh>("region0");

    // Create history file if not already created
    if (historyFilePtr_.empty())
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
                time_.timeName(mesh.time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                historyDir = time_.path()/".."/"postProcessing/rodPressure"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"postProcessing/rodPressure"/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(historyDir);

            // Open new file at start up
            // Use the function object name in the file name to allow multiple
            // objects defined on the same patch
            historyFilePtr_.reset
            (
                new OFstream
                (
                    historyDir + "/" + name + ".dat"
                )
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_() << "# Time" << " " << "rodPressure(Pa)" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::rodPressure::start()
{
    return writeData();
}

bool Foam::rodPressure::execute()
{
    return writeData();
}


bool Foam::rodPressure::read(const dictionary& dict)
{
    return true;
}

bool Foam::rodPressure::write()
{
    // return writeData();
    return true;
}

// ************************************************************************* //