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

\*---------------------------------------------------------------------------*/

#include "fpRelease.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "offbeatTime.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fpRelease, 0);
    addToRunTimeSelectionTable(functionObject, fpRelease, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fpRelease::writeData()
{
    const fvMesh& mesh = time_.lookupObject<fvMesh>(meshName_);
    const offbeatTime& runTime = refCast<const offbeatTime>(time_);

    // Auto-detect species on first call
    if (speciesNames_.empty())
    {
        const wordList allNames =
            mesh.sortedNames(uniformDimensionedScalarField::typeName);

        static const word prefix("fpRelease_");
        for (const word& n : allNames)
        {
            if (n.size() > prefix.size() && n.substr(0, prefix.size()) == prefix)
            {
                speciesNames_.append(n.substr(prefix.size()));
            }
        }

        if (speciesNames_.empty())
        {
            WarningIn("Foam::fpRelease::writeData()")
                << "No fpRelease_* fields found in mesh registry. "
                << "Is fissionProductsDiffusionSolver active?" << nl;
            return true;
        }

        // Write CSV header
        if (Pstream::master())
        {
            filePtr_() << "# Time";
            for (const word& s : speciesNames_)
                filePtr_() << " " << s;
            filePtr_() << endl;
        }
    }

    if (Pstream::master())
    {
        filePtr_() << runTime.userTime();
        for (const word& s : speciesNames_)
        {
            const word fieldName("fpRelease_" + s);
            if (mesh.foundObject<uniformDimensionedScalarField>(fieldName))
            {
                filePtr_() << " "
                    << mesh.lookupObject<uniformDimensionedScalarField>
                       (fieldName).value();
            }
            else
            {
                filePtr_() << " -1";
            }
        }
        filePtr_() << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fpRelease::fpRelease
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    meshName_(dict.lookupOrDefault<word>("region", polyMesh::defaultRegion)),
    filePtr_(),
    speciesNames_()
{
    Info<< "Creating " << this->name() << " function object" << endl;

    if (Pstream::master())
    {
        const offbeatTime& runTime = refCast<const offbeatTime>(time_);
        const word startTimeName =
            runTime.timeName(runTime.timeToUserTime(time_.startTime().value()));

        fileName historyDir;
        if (Pstream::parRun())
            historyDir = time_.path()/".."/"postProcessing/fpRelease"/startTimeName;
        else
            historyDir = time_.path()/"postProcessing/fpRelease"/startTimeName;

        mkDir(historyDir);

        filePtr_.reset(new OFstream(historyDir/name + ".dat"));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fpRelease::start()
{
    return writeData();
}

bool Foam::fpRelease::execute()
{
    return writeData();
}

bool Foam::fpRelease::read(const dictionary&)
{
    return true;
}

bool Foam::fpRelease::write()
{
    return true;
}

// ************************************************************************* //
