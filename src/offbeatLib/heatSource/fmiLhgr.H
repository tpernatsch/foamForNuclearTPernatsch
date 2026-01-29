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
    Foam::fmiLhgr

Description
    Derived from the base heatSource class, this class act as the parent class 
    for all heat source models that define the linear heat generation rate (lhgr)
    as an additional field. 

    The class typeName is "fmiLhgr", and it allows to define a lhgr 
    which come from an FMI port.

Usage
    In solverDict file:
    \verbatim
    heatSource fmiLhgr;

    heatSourceOptions
    {
        lhgrNameFromFMU     offbeatLhgr;
        initialLhgr         230e2;
        // or
        initialQavg         0.0;

        // Optional
        alpha               1.0;
        ...
    }
    \endverbatim

    where the entries mean:
    \table
        Property        | Description                   | Type   | Req'd | Dflt
        lhgrNameFromFMU | FMI port name                 | word   |  yes  | -
        initialLhgr     | Initial LHGR                  | scalar |  no   | 0.0
        initialQavg     | Initial average power density | scalar |  no   | 0.0
        alpha           | Relation factor               | scalar |  no   | 1.0
    \endtable
    
SourceFiles
    heatSource.C

\todo
    - The lhgr might be needed by some models such as the relocation, therefore
    it might be useful to have it as a field in the registry. Alternatively, 
    one could make the heatSource class an IOObject and ask directly for the lhgr 
    with a function. In this case having the lhgr as a vol field might not be necessary.

\mainauthor
    Thomas Guilbaud, EPFL/Transmutex SA, Switzerland

\date 
    October 2023

\*---------------------------------------------------------------------------*/

#ifndef fmiLhgr_H
#define fmiLhgr_H

#include "heatSource.H"
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
                         Class fmiLhgr Declaration
\*---------------------------------------------------------------------------*/

class fmiLhgr
:
    public constantLhgr
{
 
    // Private data
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        fmiLhgr(const fmiLhgr&);

        //- Disallow default bitwise assignment
        void operator=(const fmiLhgr&);
        
protected:   

    // Protected data

        //- Dictionary from which FMI variables are written
        IOdictionary fmiState_;

        //- LHGR average value, used as initial value and previous value for 
        //  implicit
        scalar lhgrAvg_;

        //- Average power density, used as initial value and previous value for 
        //  implicit
        scalar Qavg_;

        //- Alpha relation factor for lhgr and Qavg
        scalar alpha_;

        //- FMI input name for LHGR
        word lhgrNameFromFMU_;

        //- FMI input name for average power density
        word QavgNameFromFMU_;

        //- 
        bool isUseLhgr_;

        //-
        bool isUseQavg_;

public:

    //- Runtime type information
    TypeName("fmiLhgr");

    // Constructors

        //- Construct from mesh, materials and dict
        fmiLhgr
        (
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& heatSourceOptDict,
            const bool timeDependent=false
        );


    //- Destructor
    ~fmiLhgr();


    // Member Functions

        //- update the power distribution in the supplied list of cells
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
