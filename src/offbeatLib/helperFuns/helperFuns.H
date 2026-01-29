/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Useful helper functions for debugging, etc.
   
SourceFiles
    helperFuns.C

\mainauthor
    I. Clifford - PSI

\contribution

\date 
    January 2024

\*---------------------------------------------------------------------------*/

#ifndef helperFuns_H
#define helperFuns_H

#include "fieldTypes.H"
#include "fvPatch.H"
#include "string.H"
#include "MeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace helperFuns
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Debugging helper function to gather values from each processor
//- returning a Field with one value per processor.
template<class Type>
tmp<Field<Type> > gather
(
    const Type& value
)
{
    tmp<Field<Type> > tAllValues(new Field<Type>(Pstream::nProcs()));
    Field<Type>& allValues = tAllValues.ref();
    
    allValues[Pstream::myProcNo()] = value;
    Pstream::gatherList(allValues);
    return tAllValues;
}


//- Write a field of values and associated fvPatch in VTK format
void writePatchFieldVTK
(
    const primitivePatch& patch,
    const scalarField& psi,
    const word& patchName,
    const word& psiName,
    const fileName& instance
);


//- Generate a MeshedSurface<face> of a unit sphere using subdivision of a
//- tetrahedron, see http://eugene-eeo.github.io/blog/sphere-triangles.html
autoPtr<MeshedSurface<face> > meshedUnitSphere(label nDivisions=3);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace helperFuns

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
