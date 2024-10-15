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

\*---------------------------------------------------------------------------*/

#include "InterpolateTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //



const Foam::Enum<Foam::interpolateTableBase::interpolationMethod>
    interpolateTableBase::interpolationMethodNames_
({
    {interpolationMethod::STEP, "step"},
    {interpolationMethod::LINEAR, "linear"},
});

const Foam::Enum<Foam::interpolateTableBase::outofBoundsMehtod>
    interpolateTableBase::outofBoundsMehtodNames_
({
    {outofBoundsMehtod::ERROR, "error"},
    {outofBoundsMehtod::EXTRAPOLATION, "extrapolation"},
    {outofBoundsMehtod::FIXED, "fixed"},
});


// * * * * * * * * * * * * * * Member Functions ** * * * * * * * * * * * * * //

template<>
scalar scalarInterpolateTable::integral(scalar k) const
{
    scalar res = 0;
    
    for (int i=1; i < xValues_.size(); i++)
    {
        scalar dx = pow(xValues_[i], k) - pow(xValues_[i-1], k);
        
        if (method_ ==  STEP)
        {
            res += data_[i-1]*dx;
        }
        else if (method_ == LINEAR)
        {
            // Trapezoidal integration
            res += 0.5*(data_[i-1] + data_[i])*dx;
        }
    }
    
    return res;
}


//- Explicit class instantiations
template class InterpolateTable<scalarField, scalar, scalar>;
template class InterpolateTable
            <FieldField<Field, scalar>, scalarField,tmp<scalarField> >;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
