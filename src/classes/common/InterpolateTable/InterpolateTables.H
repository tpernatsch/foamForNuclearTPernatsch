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
#include "scalar.H"
#include "scalarField.H"
#include "scalarFieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- 1D interpolation table returning scalar values
typedef InterpolateTable<scalarField, scalar, scalar> scalarInterpolateTable;

//- 2D interpolation table, returning tmp<scalarField>
typedef InterpolateTable<FieldField<Field, scalar>, scalarField, tmp<scalarField> >
    scalarFieldInterpolateTable;

    //- 3D interpolation table, returning tmp<scalarField>
typedef InterpolateTable<PtrList<FieldField<Field, scalar>>, FieldField<Field, scalar>, tmp<FieldField<Field, scalar>> >
    scalarFieldFieldInterpolateTable;


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
