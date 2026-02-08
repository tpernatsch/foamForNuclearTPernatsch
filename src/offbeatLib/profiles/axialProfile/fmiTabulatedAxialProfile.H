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
    Foam::fmiTabulatedAxialProfile

Description
    Axial profile defined from an FMI.
    Both "step" and "linear" interpolation methods are supported (see 
    Foam::InterpolateTable).
    The supplied profile must be normalised to 1.

Usage
    \verbatim
    axialProfile
    {
        type                        fmiTabulated;
        axialProfileNameFromFMU     offbeatHeatSrcAxial;
        axialLocations              ( -0.5 0 0.5 1 1.5 );
        initialData                 5{1};
        axialInterpolationMethod    linear;
    }
    // or
    axialProfile // Imposed a uniform axial profile
    {
        type                        fmiTabulated;
        axialLocationsNameFromFMU   offbeatAxialLocations;
        axialProfileNameFromFMU     offbeatHeatSrcAxial;
        axialInterpolationMethod    linear;
    }
    \endverbatim


SourceFiles
    fmiTabulatedAxialProfile.C

\mainauthor
    Thomas Guilbaud, EPFL/Transmutex SA, Switzerland

\date 
    October 2023

\*---------------------------------------------------------------------------*/

#ifndef fmiTabulatedAxialProfile_H
#define fmiTabulatedAxialProfile_H

#include "axialProfile.H"
#include "InterpolateTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fmiTabulatedAxialProfile Declaration
\*---------------------------------------------------------------------------*/

class fmiTabulatedAxialProfile
:
    public axialProfile
{
    // Private Data

        //- FMI port name for axial location
        word axialLocationsNameFromFMU_;

        //- FMI port name for axial profile
        word axialProfileNameFromFMU_;

        //- Axial location interpolation method
        interpolateTableBase::interpolationMethod zMethod_;

        //- Dictionary from which FMI variables are written
        IOdictionary fmiState_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fmiTabulatedAxialProfile(const fmiTabulatedAxialProfile&);

        //- Disallow default bitwise assignment
        void operator=(const fmiTabulatedAxialProfile&);

public:

    //- Runtime type information
    TypeName("fmiTabulated");

    
    // Constructors

        //- Construct from dictionary
        fmiTabulatedAxialProfile
        (
            const fvMesh& mesh,    
            const dictionary& dict,
            const labelList& addr,
            const vector& axialDirection,
            scalar zMin,
            scalar zMax 
        );

        
    //- Destructor
    virtual ~fmiTabulatedAxialProfile();

    // Member Functions
    
        //- Update the profile
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
