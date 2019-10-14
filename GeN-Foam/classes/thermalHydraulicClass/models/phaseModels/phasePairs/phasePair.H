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
    Base class to represent a phasePair with basing name access
    functionality.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

Class
    Foam::phasePair

SourceFiles
    phasePair.C

\*---------------------------------------------------------------------------*/

#ifndef phasePair_H
#define phasePair_H

#include "volFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class fluid;

/*---------------------------------------------------------------------------*\
                           Class phasePair Declaration
\*---------------------------------------------------------------------------*/

class phasePair
    : public regIOobject
{
protected :

    //- Ref to mesh
    const fvMesh& mesh_;

    //- Name of the pair
    const word name_;

    //- Name of first phase in pair
    const word name1_;

    //- Name of second phase in pair
    const word name2_;

public:

    TypeName("phasePair");

    //- Constructors

        phasePair
        (
            const fvMesh& mesh,
            const word name1,
            const word name2
        );

        //- Destructor
        virtual ~phasePair() {}

    //- Member Functions

        //- Update dimensionless numbers
        void correct();

        //- To remove abstraction from regIOobject
        virtual bool writeData(Ostream& os) const;

        //- Access

            const fvMesh& mesh() const
            {
                return mesh_;
            }

            const word& name() const
            {
                return name_;
            }

            const word& name1() const
            {
                return name1_;
            }

            const word& name2() const
            {
                return name2_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
