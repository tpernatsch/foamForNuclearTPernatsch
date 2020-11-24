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
    N/A

Description
    Given a word "name", the function returns:

        -   if NOT on the last PIMPLE iteration:
            --  the scalar found under 
                fvSolution.relaxationFactors.fields."name" if name exists;
            --  1.0 otherwise;
        -   if on the last PIMPLE iteration:
            --  the scalar found under 
                fvSolution.relaxationFactors.fields.("name"+"Final") if name 
                exists;
            --  1.0 otherwise;

    This is the standard behaviour of the relax function of a geometric field,
    however, I need to be able to read relaxation factors in this way even
    without them being tied to a geometric field

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    None


\*---------------------------------------------------------------------------*/

#ifndef MYOPS_H
#define MYOPS_H

#include "fvMesh.H"

namespace Foam
{
namespace myOps
{
    bool relax(const fvMesh& mesh, word name);


    scalar relaxationFactor(const fvMesh& mesh, word name);
}
}

#endif

/*---------------------------------------------------------------------------*/
