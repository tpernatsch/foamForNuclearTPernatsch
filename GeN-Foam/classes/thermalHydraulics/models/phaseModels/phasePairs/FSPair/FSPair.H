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

Description
    Class to handle dimensionless numbers that characterise the system, e.g.
    Reynolds, etc. While some can (rightfully) think that some of these numbers
    should be directly members of other classes (e.g., Reynolds of fluids
    should be in the fluid class), in general these numbers depend on a
    multitude of parameters. You cannot define a Reynolds for a fluid by
    itself, it needs to be geometrically constratined in a certain domain 
    (so that a hydraulic diameter can be defined), so it is not properly
    "just" of the fluid. While this logic can of course be debate, I wanted
    to keep it simple and not start dispersing dimensionless numbers randomly
    around the code. 

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

Class
    Foam::FSPair

SourceFiles
    FSPair.C

\*---------------------------------------------------------------------------*/

#ifndef FSPair_H
#define FSPair_H

#include "phasePair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class fluid;
class structureModel;

/*---------------------------------------------------------------------------*\
                           Class FSPair Declaration
\*---------------------------------------------------------------------------*/

class FSPair
    : public phasePair
{
protected :

    //- Ref to fluid
    const fluid& fluid_;

    //- Ref to structure
    const structureModel& structure_; 

    const dimensionedScalar minRe_;

    //- Reynolds of the pair
    volScalarField Re_;

    //- Local components of the Reynolds number in the local frame of
    //  reference directions
    volVectorField lRe_;

    //- Fluid-structure contact fractions for two-phase scenarios 
    autoPtr<volScalarField> fPtr_;

public:

    TypeName("FSPair");

    //- Constructors

        FSPair
        (
            const fluid& fluid,
            const structureModel& structure
        );

        //- Destructor
        virtual ~FSPair(){}

    //- Member Functions

        //- Update dimensionless numbers
        void correct();

        //- Access

            const word& fluidName() const
            {
                return name1_;
            }

            const fluid& fluidRef() const
            {
                return fluid_;
            }

            const structureModel& structure() const
            {
                return structure_;
            }

            const volScalarField& Re() const
            {
                return Re_;
            }

            const volVectorField& lRe() const
            {
                return lRe_;
            }  

            autoPtr<volScalarField>& fPtr()
            {
                return fPtr_;
            }

            const autoPtr<volScalarField>& fPtr() const
            {
                return fPtr_;
            }  

            const volScalarField& f() const
            {
                return fPtr_();
            }  
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
