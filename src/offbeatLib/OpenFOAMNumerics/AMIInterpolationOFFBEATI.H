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

\*---------------------------------------------------------------------------*/

inline Foam::label Foam::AMIInterpolationOFFBEAT::singlePatchProc() const
{
    return singlePatchProc_;
}


inline Foam::scalar Foam::AMIInterpolationOFFBEAT::lowWeightCorrection() const
{
    return lowWeightCorrection_;
}


inline bool Foam::AMIInterpolationOFFBEAT::applyLowWeightCorrection() const
{
    return lowWeightCorrection_ > 0;
}


inline const Foam::scalarField& Foam::AMIInterpolationOFFBEAT::srcMagSf() const
{
    return srcMagSf_;
}

inline const Foam::vectorField& Foam::AMIInterpolationOFFBEAT::srcNf() const
{
    return srcNf_;
}


inline const Foam::vectorField& Foam::AMIInterpolationOFFBEAT::srcCf() const
{
    return srcCf_;
}


inline const Foam::labelListList& Foam::AMIInterpolationOFFBEAT::srcAddress() const
{
    return srcAddress_;
}


inline const Foam::scalarListList& Foam::AMIInterpolationOFFBEAT::srcWeights() const
{
    return srcWeights_;
}


inline Foam::scalarListList& Foam::AMIInterpolationOFFBEAT::srcWeights()
{
    return srcWeights_;
}


inline const Foam::scalarField& Foam::AMIInterpolationOFFBEAT::srcWeightsSum() const
{
    return srcWeightsSum_;
}


inline Foam::scalarField& Foam::AMIInterpolationOFFBEAT::srcWeightsSum()
{
    return srcWeightsSum_;
}


inline const Foam::mapDistribute& Foam::AMIInterpolationOFFBEAT::srcMap() const
{
    return srcMapPtr_();
}


inline const Foam::scalarField& Foam::AMIInterpolationOFFBEAT::tgtMagSf() const
{
    return tgtMagSf_;
}

inline const Foam::vectorField& Foam::AMIInterpolationOFFBEAT::tgtNf() const
{
    return tgtNf_;
}

inline const Foam::vectorField& Foam::AMIInterpolationOFFBEAT::tgtCf() const
{
    return tgtCf_;
}


inline const Foam::labelListList& Foam::AMIInterpolationOFFBEAT::tgtAddress() const
{
    return tgtAddress_;
}


inline const Foam::scalarListList& Foam::AMIInterpolationOFFBEAT::tgtWeights() const
{
    return tgtWeights_;
}


inline Foam::scalarListList& Foam::AMIInterpolationOFFBEAT::tgtWeights()
{
    return tgtWeights_;
}


inline const Foam::scalarField& Foam::AMIInterpolationOFFBEAT::tgtWeightsSum() const
{
    return tgtWeightsSum_;
}


inline Foam::scalarField& Foam::AMIInterpolationOFFBEAT::tgtWeightsSum()
{
    return tgtWeightsSum_;
}


inline const Foam::mapDistribute& Foam::AMIInterpolationOFFBEAT::tgtMap() const
{
    return tgtMapPtr_();
}


// ************************************************************************* //
