/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::regionCoupledOFFBEATFvPatch

Description
    Common functionality for regionCoupleFvPatch and regionCoupledWallFvPatch

SourceFiles
    regionCoupledOFFBEATFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledOFFBEATFvPatch_H
#define regionCoupledOFFBEATFvPatch_H

#include "fvPatch.H"
#include "fvMesh.H"
#include "Time.H"
#include "regionCoupledOFFBEATPolyPatch.H"
#include "regionCoupledBaseOFFBEATFvPatch.H"

#include "fvPatchFields.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class regionCoupledOFFBEATFvPatch Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledOFFBEATFvPatch
:
    public lduInterface,
    public fvPatch,
    public regionCoupledBaseOFFBEATFvPatch
{
    // Private data

        const regionCoupledOFFBEATPolyPatch& regionCoupledOFFBEATPolyPatch_;

    // Private members

        //- Return regionCoupledOFFBEATFvPatch nbr
        const regionCoupledOFFBEATFvPatch& neighbFvPatch() const
        {
            return refCast<const regionCoupledOFFBEATFvPatch>
            (
                nbrFvMesh().boundary()
                [
                    nbrPatchID()
                ]
            );
        }

public:

    //- Runtime type information
    TypeName(regionCoupledOFFBEATPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        regionCoupledOFFBEATFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm),
            regionCoupledBaseOFFBEATFvPatch
            (
                patch,
                *this
            ),
            regionCoupledOFFBEATPolyPatch_
            (
                refCast<const regionCoupledOFFBEATPolyPatch>(patch)
            )
        {}


    //- Destructor
    ~regionCoupledOFFBEATFvPatch()
    {}


    // Member Functions


        // Access

            //- Return faceCell addressing
            virtual const labelUList& faceCells() const
            {
                return fvPatch::faceCells();
            }

            //- Return true because this patch is coupled
            virtual bool coupled() const
            {
                return regionCoupledOFFBEATPolyPatch_.coupled();
            }

            const scalarField& weights() const;


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

#ifdef OPENFOAMESI
            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList &  internalData,
                const labelUList &  faceCells
            ) const;
#endif            

            //- Inherit initInternalFieldTransfer from lduInterface
            using lduInterface::initInternalFieldTransfer;

            //- Initialise neighbour field transfer
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                labelUList& iF
            ) const
            {}

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
