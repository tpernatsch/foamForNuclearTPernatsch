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
    Function to split a word along a char and return a word list. While this is
    definitely doable in standard OpenFOAM (with the functions in the stringOps
    namespace), for some reason, these do not appear to always work on the CSCS
    supercomputing cluster-compiled version of the solver. For sure it is a 
    compiler issue, but too deep to track, yet nonetheless related to the use 
    of stringOps.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    None


\*---------------------------------------------------------------------------*/

#ifndef MYSTRINGOPS_H
#define MYSTRINGOPS_H

#include "List.H"

namespace Foam
{
namespace myStringOps
{
    //- The deliminter is a char, so it must be provided in single-quotation
    //  marks. Type can be either word or string
    template<class Type>
    List<Type> split(const Type& compoundItem, char d)
    {
        List<Type> items(1, "");
        int i = 0;
        for (char const &c : compoundItem)
        {
            if (c == d)
            {
                i += 1;
                items.append("");
            }
            else items[i] += c;
        }
        return items;
    }
}
}

#endif

/*---------------------------------------------------------------------------*/
