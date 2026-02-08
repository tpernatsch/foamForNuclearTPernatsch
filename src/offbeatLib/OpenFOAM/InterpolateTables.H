/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    A. Scolaro, E. Brunetto, C. Fiorina - EPFL (ECOLE POLYTECHNIQUE FEDERALE 
    DE LAUSANNE, Switzerland, Laboratory for Reactor Physics and Systems 
    Behaviour)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef InterpolateTables_H
#define InterpolateTables_H

#include "InterpolateTable.H"
#include "primitiveFields.H"
#include "FieldFields.H"
#include "fieldTypes.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#define defineInterpolateTable(type, nullArg)                               \
typedef InterpolateTable<type##Field, type, type> type##InterpolateTable;   \
typedef InterpolateTable<FieldField<Field, type>, type##Field,              \
                         tmp<type##Field> > type##FieldInterpolateTable;    \
typedef InterpolateTable<PtrList<FieldField<Field, type>>,                  \
                         FieldField<Field, type>,                           \
                         tmp<FieldField<Field, type>>>                      \
                         type##FieldFieldInterpolateTable;
    

FOR_ALL_FIELD_TYPES(defineInterpolateTable, type);


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
