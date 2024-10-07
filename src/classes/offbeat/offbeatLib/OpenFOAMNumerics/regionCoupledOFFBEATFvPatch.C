/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "regionCoupledOFFBEATFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionCoupledOFFBEATFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, regionCoupledOFFBEATFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::regionCoupledOFFBEATFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}

#ifdef OPENFOAMESI
Foam::tmp<Foam::labelField> Foam::regionCoupledOFFBEATFvPatch::interfaceInternalField
(
    const labelUList &  internalData,
    const labelUList &  faceCells
) const
{
    return patchInternalField(internalData); // , faceCells);
}
#endif



Foam::tmp<Foam::labelField> Foam::regionCoupledOFFBEATFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    if (neighbFvPatch().sameRegion())
    {
        return neighbFvPatch().patchInternalField(iF);
    }
    else
    {
        return tmp<labelField>(new labelField(iF.size(), 0));

    }

    return tmp<labelField>(nullptr);
}
// ************************************************************************* //
