/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::FSDragCoefficientModels::Engel

Description
    Engel correlation for pressure drop

    @article{doi:10.13182/NSE79-A20618,
    author = {F. C. Engel and R. A. Markley and A. A. Bishop},
    title = {   Laminar, Transition, and Turbulent Parallel Flow Pressure Drop 
                Across Wire-Wrap-Spaced Rod Bundles},
    journal = {Nuclear Science and Engineering},
    volume = {69},
    number = {2},
    pages = {290-296},
    year  = {1979},
    publisher = {Taylor & Francis},
    doi = {10.13182/NSE79-A20618},

    URL = { 
            https://doi.org/10.13182/NSE79-A20618
        
    },
    eprint = { 
            https://doi.org/10.13182/NSE79-A20618
        
    }
    ,
        abstract = { 
            AbstractLaminar, transition, and turbulent parallel flow pressure 
            drop across wire-wrapped hexagonal rod bundles positioned inside a 
            duct were determined in tests using water, sodium, and air. A 
            smooth transition region from turbulent to laminar flow that 
            occurred over the Reynolds number range from 5000 to 400 
            characterized the resulting friction factor behavior. The 
            continuous transition region could be explained in terms of the 
            fraction of the flow area in turbulent flow.Laminar friction 
            factors calculated from individual subchannel measurements could be
            correlated by the same expression found for rod-bundle-averaged 
            conditions. In the laminar range, the friction factor was 
            correlated by the expression f = 110/Re, in the turbulent range by 
            f = 0.55/Re0.25, and in the transition range bywhere is the 
            intermittency factor. A general laminar flow friction factor 
            correlation was developed:This correlation agrees satisfactorily 
            with limited laminar flow data from rod bundles having different 
            wire-wrap lead pitch-to-diameter ratios. }
    }

Author
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    Engel.C

\*---------------------------------------------------------------------------*/

#ifndef EngelFSDragCoefficient_H
#define EngelFSDragCoefficient_H

#include "FSDragCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class FSPair;

namespace FSDragCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class Engel Declaration
\*---------------------------------------------------------------------------*/

class Engel
:
    public FSDragCoefficientModel
{
public:

    TypeName("Engel");

    //- Constructors

        //- Construct for a fluid-structure pair
        Engel
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~Engel(){}

    //- Member Functions

        //- Drag coefficient
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FSDragCoefficientModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //