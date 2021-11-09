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
    Foam::regimeMapModel::twoParameters::regimeDomain2D

Description
    Regime domain representation in 2D parameter-space.

Autor
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    twoParameters.C

\*---------------------------------------------------------------------------*/

#ifndef regimeDomain2D_H
#define regimeDomain2D_H

#include "regimeBoundary2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace regimeMapModels
{

/*---------------------------------------------------------------------------*\
                           Class regimeDomain2D Declaration
\*---------------------------------------------------------------------------*/

class regimeDomain2D
{
private:

    //- Domain label
    label id_;

    //- List of points that make up the domain (relatively useless)
    List<const Vector2D<scalar>*> pointPtrs_;

    //- List of boundaries that make up the domain
    List<regimeBoundary2D> boundaries_;

    //- List of boundaries in boundaries_ that are shared by this and another
    //  regimeDomain. Used for interpolation purposes
    List<const regimeBoundary2D*> internalBoundaryPtrs_;

public:

    //- Runtime type information
    TypeName("regimeDomain2D");

    //- Constructors

        //- Default
        regimeDomain2D();

        //- Constructor
        regimeDomain2D
        (
            label id,
            List<const Vector2D<scalar>*> pointPtrs
        );

        //- Destructor
        virtual ~regimeDomain2D(){}

    //- Member Functions

        //-
        Tuple2<Vector2D<scalar>, Vector2D<scalar>> boundingBox() const;

        //- 
        void setInternalBoundaries();

        //- Check whether point p is contained within the domain
        bool containsPoint(const Vector2D<scalar>& p) const;

        //- Access

            const label& id() const
            {
                return id_;
            }

            //- Non const-access as well to allow the use of the setNbr 
            //  function of regimeBoundary2D
            List<regimeBoundary2D>& boundaries()
            {
                return boundaries_;
            }

            const List<regimeBoundary2D>& boundaries() const
            {
                return boundaries_;
            }

            const List<const regimeBoundary2D*>& internalBoundaryPtrs() const
            {
                return internalBoundaryPtrs_;
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regimeMapModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //