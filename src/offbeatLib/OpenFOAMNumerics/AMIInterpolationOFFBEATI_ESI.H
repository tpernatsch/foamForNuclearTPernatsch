/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

inline const Foam::primitivePatch& Foam::AMIInterpolationOFFBEAT::srcPatch0() const
{
    if (!tsrcPatch0_.valid())
    {
        FatalErrorInFunction
            << "tsrcPatch0Ptr_ not set"
            << abort(FatalError);
    }

    return tsrcPatch0_();
}


inline const Foam::primitivePatch& Foam::AMIInterpolationOFFBEAT::tgtPatch0() const
{

    if (!ttgtPatch0_.valid())
    {
        FatalErrorInFunction
            << "ttgtPatch0Ptr_ not set"
            << abort(FatalError);
    }

    return ttgtPatch0_();
}


inline bool Foam::AMIInterpolationOFFBEAT::upToDate() const
{
    return upToDate_;
}


inline bool& Foam::AMIInterpolationOFFBEAT::upToDate()
{
    return upToDate_;
}


inline bool Foam::AMIInterpolationOFFBEAT::distributed() const
{
    return singlePatchProc_ == -1;
}


inline bool Foam::AMIInterpolationOFFBEAT::requireMatch() const
{
    return requireMatch_;
}


inline bool Foam::AMIInterpolationOFFBEAT::setRequireMatch(const bool flag)
{
    requireMatch_ = flag;
    return requireMatch_;
}


inline bool Foam::AMIInterpolationOFFBEAT::mustMatchFaces() const
{
    return requireMatch_ && !applyLowWeightCorrection();
}


inline bool Foam::AMIInterpolationOFFBEAT::reverseTarget() const
{
    return reverseTarget_;
}


inline Foam::scalar Foam::AMIInterpolationOFFBEAT::lowWeightCorrection() const
{
    return lowWeightCorrection_;
}


inline bool Foam::AMIInterpolationOFFBEAT::applyLowWeightCorrection() const
{
    return lowWeightCorrection_ > 0;
}


inline Foam::label Foam::AMIInterpolationOFFBEAT::singlePatchProc() const
{
    return singlePatchProc_;
}


inline const Foam::List<Foam::scalar>& Foam::AMIInterpolationOFFBEAT::srcMagSf() const
{
    return srcMagSf_;
}


inline Foam::List<Foam::scalar>& Foam::AMIInterpolationOFFBEAT::srcMagSf()
{
    return srcMagSf_;
}


inline const Foam::labelListList& Foam::AMIInterpolationOFFBEAT::srcAddress() const
{
    return srcAddress_;
}


inline Foam::labelListList& Foam::AMIInterpolationOFFBEAT::srcAddress()
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


inline const Foam::pointListList& Foam::AMIInterpolationOFFBEAT::srcCentroids() const
{
    return srcCentroids_;
}


inline Foam::pointListList& Foam::AMIInterpolationOFFBEAT::srcCentroids()
{
    return srcCentroids_;
}


inline const Foam::mapDistribute& Foam::AMIInterpolationOFFBEAT::srcMap() const
{
    return *srcMapPtr_;
}


inline const Foam::List<Foam::scalar>& Foam::AMIInterpolationOFFBEAT::tgtMagSf() const
{
    return tgtMagSf_;
}


inline Foam::List<Foam::scalar>& Foam::AMIInterpolationOFFBEAT::tgtMagSf()
{
    return tgtMagSf_;
}


inline const Foam::labelListList& Foam::AMIInterpolationOFFBEAT::tgtAddress() const
{
    return tgtAddress_;
}


inline Foam::labelListList& Foam::AMIInterpolationOFFBEAT::tgtAddress()
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
    return *tgtMapPtr_;
}

// ************************************************************************* //
