/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::regionCoupledBaseOFFBEATFvPatch

Description
    Base class of regionCoupledFvPatch with common functionality for
    regionCoupledFvPatch and regionCoupledWallFvPatch

    OFFBEAT-version of this class contains the possibility to update the AMI
    addressing during the simulation.

SourceFiles
    regionCoupledBaseOFFBEATFvPatch.C

\mainauthorMod
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contributionMod
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledBaseOFFBEATFvPatch_H
#define regionCoupledBaseOFFBEATFvPatch_H

#include "regionCoupledBaseOFFBEAT.H"
#include "fvMesh.H"
#include "Time.H"
#include "polyPatch.H"
#include "regionCoupledOFFBEATLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class regionCoupledBaseOFFBEATFvPatch Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledBaseOFFBEATFvPatch
:
    public regionCoupledOFFBEATLduInterface
{
protected:

    // Protected data

        //- Reference to regionCoupledBase
        const regionCoupledBaseOFFBEAT& regionCoupledBaseOFFBEAT_;

        //- Reference to fvPatch
        const fvPatch& patch_;


    // Protected members

        //- Returns fvMesh
        const fvMesh& nbrFvMesh() const
        {
            return
            (
                patch_.boundaryMesh().mesh().time().lookupObject<fvMesh>
                (
                    regionCoupledBaseOFFBEAT_.nbrRegionName()
                )
            );
        }


public:

    //- Runtime type information
    TypeName("regionCoupledOFFBEAT");


    // Constructors

        //- Construct from polyPatch
        regionCoupledBaseOFFBEATFvPatch
        (
            const polyPatch& pp,
            const fvPatch& patch
        )
        :
            regionCoupledOFFBEATLduInterface(),
            regionCoupledBaseOFFBEAT_
            (
                refCast<const regionCoupledBaseOFFBEAT>(pp)
            ),
            patch_(patch)
        {}


    //- Destructor
    virtual ~regionCoupledBaseOFFBEATFvPatch()
    {}


    // Member Functions

        // Access

            //- Return neighbour
            virtual label nbrPatchID() const
            {
                return regionCoupledBaseOFFBEAT_.nbrPatchID();
            }

            //- Is it the owner?
            virtual bool owner() const
            {
                return regionCoupledBaseOFFBEAT_.owner();
            }

            //- Return regionCoupledBase neighb Patch
            virtual const regionCoupledBaseOFFBEATFvPatch& nbrPatch() const
            {
                return refCast<const regionCoupledBaseOFFBEATFvPatch>
                (
                    nbrFvMesh().boundary()
                    [
                        regionCoupledBaseOFFBEAT_.nbrPatchID()
                    ]
                );
            }

            //- Return a reference to the AMI interpolator
            virtual const AMIInterpolation& AMI() const
            {
                return regionCoupledBaseOFFBEAT_.AMI();
            }

            //- Returns neighbour polyMesh
            virtual const polyMesh& nbrMesh() const
            {
                return
                (
                    patch_.boundaryMesh().mesh().time().lookupObject<polyMesh>
                    (
                        regionCoupledBaseOFFBEAT_.nbrRegionName()
                    )
                );
            }

            //- Return fvPatch
            const fvPatch& patch() const
            {
                return patch_;
            }

            //- Returns if it is the same Region
            bool sameRegion() const
            {
                return regionCoupledBaseOFFBEAT_.sameRegion();
            }

            //- Return regionCoupledPolyPatch
            const regionCoupledBaseOFFBEAT& regionCoupledPatch() const
            {
                return regionCoupledBaseOFFBEAT_;
            }

            //- Return neighbor fvPatch
            const fvPatch& neighbFvPatch() const
            {
                return refCast<const fvPatch>
                (
                    nbrFvMesh().boundary()
                    [
                        regionCoupledBaseOFFBEAT_.nbrPatchID()
                    ]
                );
            }

            //- Return the interface type
            const word& regionCoupleType() const
            {
                return regionCoupledBaseOFFBEAT_.regionCoupleType();
            }


             //- Return faceCell addressing
            virtual const labelUList& faceCells() const = 0;

            //- Update AMI
            void updateAMI() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
