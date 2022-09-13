/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::timeProfile

Description
    Base class for time profile functions.

SourceFiles
    timeProfile.C

Author of this file (and associated .C or included .H files)
    Thomas Guilbaud <cthomas.guilbaud@epfl.ch;>, EPFL (Swizerland)


\*---------------------------------------------------------------------------*/

#ifndef timeProfile_H
#define timeProfile_H

#include "volFields.H"
#include "rhoThermo.H"
#include "fvMatrices.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "coordinateSystem.H"
#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "volPointInterpolation.H"
#include "meshToMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class timeProfile Declaration
\*---------------------------------------------------------------------------*/

class timeProfile
{

protected:

    //- Dictionary name
    // word dictName_;

    //- Dictionary
    dictionary dict_;

    //- Type of the function used for calculation (e.g. table, ...)
    word type_;

    //- Time offset
    scalar startTime_;

    //- Function pointer
    autoPtr<Function1<scalar>> functionPtr_;

public:

    //- Constructors
        //- Disallow default bitwise copy construct
        timeProfile(const timeProfile&);

        //- Disallow default bitwise assignment
        void operator=(const timeProfile&);

        //- Constructor from dictionary
        timeProfile(const dictionary& dict);

        //- Constructor from IOobject and subDictname
        timeProfile(IOdictionary object, word timeProfileName);

    //- Destructor
        virtual ~timeProfile(){}


    //- Member functions

        //-
        bool valid();

        //- Get the value at a given time minus the start time offset.
        scalar value(scalar time);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
