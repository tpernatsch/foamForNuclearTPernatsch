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
    List of various typedefs that appear throughout various classes.
    Rather that re-type-deffing the same type over and over in each file,
    I preferred to include all of them in a single file.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    None



\*---------------------------------------------------------------------------*/

#include"HashPtrTable.H"

namespace Foam
{
    //- Forward declarations
    class regime;

    //- Typedefs
    typedef HashTable
    <
        labelList,
        word,
        word::hash
    >   labelListTable;

    typedef HashTable
    <
        scalarField,
        word,
        word::hash
    >   scalarFieldTable;

    typedef HashTable
    <
        volScalarField,
        word,
        word::hash
    >   volScalarFieldTable;

    typedef HashPtrTable
    <
        volScalarField, 
        word, 
        word::hash
    > volScalarFieldPtrTable;

    typedef HashTable
    <
        volTensorField,
        word,
        word::hash
    >   volTensorFieldTable;

    typedef HashPtrTable
    <
        volTensorField,
        word,
        word::hash
    >   volTensorFieldPtrTable;

    typedef HashTable
    <
        autoPtr<regime>,
        word,
        word::hash
    > regimeTable;
}

/*---------------------------------------------------------------------------*/
