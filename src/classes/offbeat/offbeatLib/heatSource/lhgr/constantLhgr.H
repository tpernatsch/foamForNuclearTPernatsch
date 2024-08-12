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
    Foam::constantLhgr

Description
    Derived from the base heatSource class, this class act as the parent class 
    for all heat source models that define the linear heat generation rate (lhgr)
    as an additional field. 

    \htmlinclude heatSource_constantLhgr.html
    
SourceFiles
    heatSource.C

\todo
    - The lhgr might be needed by some models such as the relocation, therefore
    it might be useful to have it as a field in the registry. Alternatively, 
    one could make the heatSource class an IOObject and ask directly for the lhgr 
    with a function. In this case having the lhgr as a vol field might not be necessary.

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef constantLhgr_H
#define constantLhgr_H

#include "heatSource.H"
#include "axialProfile.H"
#include "radialProfile.H"
#include "azimuthalProfile.H"
#include "volFields.H"
#include "InterpolateTables.H"
#include "globalOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantLhgr Declaration
\*---------------------------------------------------------------------------*/

class constantLhgr
:
    public heatSource
{
 
    // Private data
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantLhgr(const constantLhgr&);

        //- Disallow default bitwise assignment
        void operator=(const constantLhgr&);
        
protected:   

    // Protected data
        
        //- Linear heat generation rate (W/m)
        volScalarField lhgr_;

        //- Constant lhgr average value
        scalar constLhgr_;
        
        //- Pin axial direction
        vector pinDirection_;
        
        //- Radial direction
        vector radialDirection_;
        
        //- Minimum z
        scalar zMin_;
        
        //- Maximum z
        scalar zMax_;

        //- Maximum radius
        scalar rMax_;

        //- Time-dependent radial profile
        autoPtr<radialProfile> radialModel_;

        //- Time-dependent azimuthal profile
        autoPtr<azimuthalProfile> azimuthalModel_;
        
        //- Time-dependent axial profile
        autoPtr<axialProfile> axialModel_;

        //- Addressing for power-producing cells
        labelList addr_;

        //- Bool to see if update from PK is necessary
        const bool updateFromPK_;

        scalar powerOld_;

public:

    //- Runtime type information
    TypeName("constantLhgr");

    // Constructors

        //- Construct from mesh, materials and dict
        constantLhgr
        (
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& heatSourceOptDict,
            const bool timeDependent=false
        );


    //- Destructor
    ~constantLhgr();


    // Member Functions

        //- Get addressing for power-producing cells
        void calcAddressing(const dictionary& dict);
        
        //- Calculate the fuel height
        void calcReferenceDimensions(const dictionary& dict);
    
        //- update the power distribution in the supplied list of cells
        virtual void correct();
        
        //- Return the next average power density (assuming the deltaT stays the
        // same)
        virtual scalar predictNextPowerDensity() const;    

        //- Return the average power density 
        virtual scalar averagePowerDensity() const;

        //- Reshape Q profile according to porosity
        virtual void reshapePowerPorosity();

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

