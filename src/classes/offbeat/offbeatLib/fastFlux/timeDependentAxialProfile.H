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
    Foam::timeDependentAxialProfile

Description
    Derived from the base constantFastFlux class, this class considers a time-dependent
    rod-average fast flux with possible axial profile.    
    The class is provided with a mesh and addressing array, and then deals with
    the mapping of the profile onto an addressed list, which can be accessed
    through the profile function.      

    \htmlinclude fastFlux_timeDependentAxialProfile.html

SourceFiles
    timeDependentAxialProfile.C

\todo
    - material zones provided or decided based on material.isFuel() or everywhere
    - do we need also thermalFlux classes?
    - do we need the profile to be normalized?
    - do we need a radial profile for the fast/thermal flux?    

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef timeDependentAxialProfile_H
#define timeDependentAxialProfile_H

#include "constantFastFlux.H"
#include "axialProfile.H"
#include "radialProfile.H"
#include "volFields.H"
#include "InterpolateTables.H"
#include "globalOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class timeDependentAxialProfile Declaration
\*---------------------------------------------------------------------------*/

class timeDependentAxialProfile
:
    public constantFastFlux
{
    // Private typedefs
    
        typedef interpolateTableBase::interpolationMethod interpolationMethod;
 
    // Private data
        
        //- Tabulated time-dependent fast flux (n/cm2/s)
        scalarField fluxData_;

        //- Tabulated time values
        scalarField timeData_;
        
        interpolationMethod method_;

        //- Time-dependent rod-average fast flux (n/cm2/s)
        scalarInterpolateTable fluxTable_;
        
        //- Pin axial direction
        vector pinDirection_;
        
        //- Minimum z
        scalar zMin_;
        
        //- Maximum z
        scalar zMax_;
        
        //- Time-dependent axial profile
        autoPtr<axialProfile> axialModel_;

        //- Addressing for fast flux cells
        labelList addr_;
    
    // Private Member Functions

        //- Get addressing for fast flux cells
        void calcAddressing(const dictionary& dict);
        
        //- Calculate the fuel-average cross-sectional area
        void calcMinMaxZ(const dictionary& dict);
        
        //- Disallow default bitwise copy construct
        timeDependentAxialProfile(const timeDependentAxialProfile&);

        //- Disallow default bitwise assignment
        void operator=(const timeDependentAxialProfile&);

public:

    //- Runtime type information
    TypeName("timeDependentAxialProfile");

    // Constructors

        //- Construct from mesh, materials and dict
        timeDependentAxialProfile
        (
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& fastFluxOptDict
        );


    //- Destructor
    ~timeDependentAxialProfile();


    // Member Functions
        
        //- Update the fastFlux distribution in the supplied list of cells
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
