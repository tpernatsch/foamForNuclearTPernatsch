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

Class
    Foam::regionCoupledBaseOFFBEAT

Description
    Base class with common functinality for regionCoupled polyPatch.
    It includes AMI.


SourceFiles
    regionCoupledBaseOFFBEAT.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledBaseOFFBEAT_H
#define regionCoupledBaseOFFBEAT_H

#ifdef OPENFOAMFOUNDATION
#include "AMIInterpolation.H"
#elif OPENFOAMESI
#include "AMIInterpolation.H"
#include "meshTools.H"
#endif

#include "polyBoundaryMesh.H"
#include "volFields.H"
#include "PrimitivePatchInterpolation.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class regionCoupledBaseOFFBEAT Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledBaseOFFBEAT
{
private:

    // Private data

        //- Patch to sample
        const polyPatch& patch_;

        // Owner or not NOTE:: difference with standard rCoupledBase
        const bool owner_;

        //- Name of other half
        const word nbrPatchName_;

        //- Index of other half
        mutable label nbrPatchID_;

        //- Name of other region
        const word nbrRegionName_;

        //- Same region
        const bool sameRegion_;

        //- AMI interpolation class
        mutable autoPtr<AMIInterpolation> AMIPtr_;

        //- Flag to indicate that slave patch should be reversed for AMI
        const bool AMIReverse_;

         //- Projection surface
        mutable autoPtr<searchableSurface> surfPtr_;

        //- Dictionary used during projection surface construction
        const dictionary surfDict_;

        //- Update AMI weights or not
        const bool updateAMI_;

        //- AMI method
        const word methodName_;


    // Private Member Functions

        //- Reset the AMI interpolator
        void resetAMI() const;


protected:


        //- Clear geometry
        virtual void clearGeom();


public:

    //- Runtime type information
    TypeName("regionCoupledOFFBEAT");


    // Constructors

        //- Construct from patch
        regionCoupledBaseOFFBEAT(const polyPatch&);

        //- Construct from dictionary
        regionCoupledBaseOFFBEAT(const polyPatch&, const dictionary&);

         //- Construct as copy, resetting patch
        regionCoupledBaseOFFBEAT(const polyPatch&, const regionCoupledBaseOFFBEAT&);


    //- Destructor
    virtual ~regionCoupledBaseOFFBEAT();


    // Member Functions


        // Access

            //- Neighbour patch name
            const word& nbrPatchName() const
            {
                return nbrPatchName_;
            }

            //- Neighbour region name
            const word& nbrRegionName() const
            {
                return nbrRegionName_;
            }

            //- Cached sampleRegion != mesh.name()
            bool sameRegion() const
            {
                return sameRegion_;
            }

            //- Neighbour patch ID
            label nbrPatchID() const;

            //- Does this side own the patch?
            bool owner() const;

            //- Return a reference to the neighbour patch
            const regionCoupledBaseOFFBEAT& nbrPatch() const;

            //- Return a reference to the projection surface
            const autoPtr<searchableSurface>& surfPtr() const;

            //- Return a reference to the AMI interpolator
            const AMIInterpolation& AMI() const;


            // Interpolations

                //- Interpolate field
                template<class Type>
                tmp<Field<Type>> interpolate(const Field<Type>& fld) const;

                //- Interpolate tmp field
                template<class Type>
                tmp<Field<Type>> interpolate
                (
                    const tmp<Field<Type>>& tFld
                ) const;

                //- Low-level interpolate List
                template<class Type, class BinaryOp>
                void interpolate
                (
                    const UList<Type>& fld,
                    const BinaryOp& bop,
                    List<Type>& result
                ) const;


        //- Return new ordering for primitivePatch.
        //  Ordering is -faceMap: for every face
        //  index of the new face -rotation:for every new face the clockwise
        //  shift of the original face. Return false if nothing changes
        //  (faceMap is identity, rotation is 0), true otherwise.
        virtual bool order
        (
            PstreamBuffers&,
            const primitivePatch&,
            labelList& faceMap,
            labelList& rotation
        ) const;


        //- Return the type
        virtual const word& regionCoupleType() const = 0;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;

        //- Update AMI weights
        bool updateAMI() const
        {
            return updateAMI_;
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "regionCoupledBaseOFFBEATTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
