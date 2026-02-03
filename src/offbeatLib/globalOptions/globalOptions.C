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

#include "globalOptions.H"
#include "wedgePolyPatch.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::globalOptions::calcAngularFraction
(
    const dictionary& dict
) const
{
    angularFraction_ = globalOptDict_.lookupOrDefault("angularFraction", 0.0);
    if(angularFraction_ > 0.0)
    {
        Info << "Angular fraction provided by user in globalOptions dict:" << nl
        << tab << angularFraction_<< nl << endl;
        return;
    }
    //- If the angular fraction is not provided, look for wedge patches.
    //- TODO Current algorithm assumes that wedge patches are symmetric wrt
    //- horizontal axis. Is this always the case?
    else
    {
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();
        forAll(mesh_.boundaryMesh(), i)
        {
            if
            (
                patches[i].size()
                && isType<wedgePolyPatch>(patches[i])
            )
            {
                 const wedgePolyPatch& pp =
                     refCast<const wedgePolyPatch>(patches[i]);
                scalar wedgeAngle = acos(pp.cosAngle());
                angularFraction_ = 2*radToDeg(wedgeAngle)/360;
                Info << "The model is axisymmetric (wedge patches are used)." << nl
                << "Angular fraction calculated by "
                << "timeDependentRadialAxialProfile::calcAngularFraction:" << nl
                << tab << angularFraction_<< nl << endl;
                return;
            }
        }
        FatalErrorIn("timeDependentRadialAxialProfile::calcAngularFraction")
            << "The model does not contain wedge patches and automatic " << nl
            << "calculation of model angularFraction is not possible." << nl
            << "Please provide angularFraction keyword in the globalOptions dict"
            << exit(FatalError);
    }

    if (angularFraction_ <= 0 or angularFraction_ > 1)
    {        
        FatalErrorIn("timeDependentRadialAxialProfile::calcAngularFraction")
        << "The angular fraction of " << angularFraction_
        << " is outside the allowed range (0,1]." << nl
        << "Check the angularFraction keyword in the globalOptions dict"
        << " or the model geometry."
        << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalOptions::globalOptions
(
    const fvMesh& mesh,
    const dictionary& solverDict
)
:
    regIOobject
    (
        IOobject
        (
            "globalOptions",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    globalOptDict_(solverDict.subOrEmptyDict("globalOptions")),
    reactorType_(globalOptDict_.lookupOrDefault<word>("reactorType", "LWR")),
    pinDirection_(globalOptDict_.lookupOrDefault<vector>("pinDirection", Foam::vector(0, 0, 1))),
    radialDirection_(globalOptDict_.lookupOrDefault<vector>("radialDirection", Foam::vector(1, 0, 0))),
    angularFraction_(-1)
{
    pinDirection_ /= mag(pinDirection_);
    radialDirection_ /= mag(radialDirection_);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
