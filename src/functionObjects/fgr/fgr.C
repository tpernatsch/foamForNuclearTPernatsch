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

#include "fgr.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "OSspecific.H"
#include "fgrSCIANTIX.H"
#include "offbeatTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fgr, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        fgr,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fgr::writeData()
{
    const fvMesh& mesh = time_.lookupObject<fvMesh>("region0");
    const word objName("SCIANTIX_Fields");

    const offbeatTime& runTime = refCast<const offbeatTime>(time_);

    if(mesh.foundObject<fgrSCIANTIX>(objName))
    {
        const scalar fgr = (fgrM3_) 
        ? mesh.lookupObject<fgrSCIANTIX>(objName).fgrM3()
        : mesh.lookupObject<fgrSCIANTIX>(objName).fgrPercent();

        // Write to file
        if (Pstream::master())
        {
            historyFilePtr_() << runTime.userTime() << " " << fgr << endl;
        }
    }
    else
    {
        WarningIn("Foam::fgr::writeData()") << nl  
        << "    Cannot find fgrSCIANTIX object for postprocessing !" 
        <<  nl << endl;

        // Write to file
        if (Pstream::master())
        {
            historyFilePtr_() << runTime.userTime() << " " << "-1" << endl;
        }
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fgr::fgr
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    fgrM3_(dict.lookupOrDefault<bool>("fgrM3", false)),
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

            const offbeatTime& runTime = refCast<const offbeatTime>(time_);
            word startTimeName =
            runTime.timeName(runTime.timeToUserTime(time_.startTime().value()));

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                historyDir = time_.path()/".."/"postProcessing/fgr"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"postProcessing/fgr"/startTimeName;
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
                if (fgrM3_)
                {
                    historyFilePtr_() << "# Time" << " " << "fgr(m3)" << endl;
                }
                else
                {
                    historyFilePtr_() << "# Time" << " " << "fgr(%)" << endl;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fgr::start()
{
    return writeData();
}

bool Foam::fgr::execute()
{
    return writeData();
}


bool Foam::fgr::read(const dictionary& dict)
{
    return true;
}

bool Foam::fgr::write()
{
    // return writeData();
    return true;
}

// ************************************************************************* //