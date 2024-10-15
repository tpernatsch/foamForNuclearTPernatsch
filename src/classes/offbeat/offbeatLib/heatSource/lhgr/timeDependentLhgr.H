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
    Foam::timeDependentLhgr

Description
    Derived from the base constantLhgr class, this class considers a time-dependent
    linear heat generation rate (lhgr).

    \htmlinclude heatSource_timeDependentLhgr.html   
    
SourceFiles
    timeDependentLhgr.C

\todo
    - The conversion from lhgr to volumetric power might be wrong when fuel with 
    different dimension (e.g. hollow and full pellet in the same rod) are used

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef timeDependentLhgr_H
#define timeDependentLhgr_H

#include "constantLhgr.H"
#include "axialProfile.H"
#include "radialProfile.H"
#include "volFields.H"
#include "InterpolateTables.H"
#include "globalOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class timeDependentLhgr Declaration
\*---------------------------------------------------------------------------*/

class timeDependentLhgr
:
    public constantLhgr
{
    // Private typedefs
    
        typedef interpolateTableBase::interpolationMethod interpolationMethod;
 
    // Private data
        
        //- Tabulated time-dependent rod-average LHGR (W/m)
        scalarField lhgrData_;

        //- Tabulated time values
        scalarField timeData_;
        
        interpolationMethod method_;

        //- Time-dependent table for rod-average LHGR (W/m)
        scalarInterpolateTable lhgrTable_;

        //- Maximum allowed Lhgr increase for adjustable time step
        scalar maxAbsoluteLhgrIncrease_;
        
    // Private Member Functions

        //- Lookup the LHGR from the supplied physical time in seconds
        inline scalar lookupLHGR(scalar t) const
        {
            return lhgrTable_(mesh_.time().timeToUserTime(t));
        }
        
        //- Disallow default bitwise copy construct
        timeDependentLhgr(const timeDependentLhgr&);

        //- Disallow default bitwise assignment
        void operator=(const timeDependentLhgr&);

public:

    //- Runtime type information
    TypeName("timeDependentLhgr");

    // Constructors

        //- Construct from mesh, materials and dict
        timeDependentLhgr
        (
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& heatSourceOptDict
        );


    //- Destructor
    ~timeDependentLhgr();


    // Member Functions
    
        //- update the power distribution in the supplied list of cells
        virtual void correct();
        
        //- Return the next average power density (assuming the deltaT stays the
        // same)
        virtual scalar predictNextPowerDensity() const;    

        //- Provide next deltaT due to absolute change in heat source (not relative)
        virtual scalar nextDeltaTAbsoluteChange(scalar deltaT) const;   

        //- For time-dependent heat sources, return the next time marker,
        //- i.e. the next time point at which the rate of power change
        //- will be discontinuous or non-smooth.
        virtual scalar nextTimeMarker() const;

        //- Return the final time marker after which the power distribution
        //- is no longer defined.
        virtual scalar lastTimeMarker() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
