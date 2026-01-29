/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::regionCoupledBaseOFFBEATGAMGInterface

Description
    Base class for GAMG agglomerated coupled region interface.

SourceFiles
    regionCoupledBaseOFFBEATGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledBaseOFFBEATGAMGInterface_H
#define regionCoupledBaseOFFBEATGAMGInterface_H

#include "GAMGInterface.H"
#include "regionCoupledOFFBEATLduInterface.H"
#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class regionCoupledBaseOFFBEATGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledBaseOFFBEATGAMGInterface
:
    public GAMGInterface,
    virtual public regionCoupledOFFBEATLduInterface
{
    // Private data

        //- Reference for the regionCoupledBaseFvPatch from which this is
        //  agglomerated
        const regionCoupledOFFBEATLduInterface& fineregionCoupledOFFBEATLduInterface_;

        //- AMI interface
        //autoPtr<AMIPatchToPatchInterpolation> amiPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionCoupledBaseOFFBEATGAMGInterface(const regionCoupledBaseOFFBEATGAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const regionCoupledBaseOFFBEATGAMGInterface&);


public:


    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        regionCoupledBaseOFFBEATGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& restrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );


    //- Destructor
    virtual ~regionCoupledBaseOFFBEATGAMGInterface();


    // Member Functions

        // Interface transfer functions

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;


        //- Cyclic interface functions

            //- Return neigbour processor number
            virtual label nbrPatchID() const
            {
                return fineregionCoupledOFFBEATLduInterface_.nbrPatchID();
            }

            virtual bool owner() const
            {
                return fineregionCoupledOFFBEATLduInterface_.owner();
            }

            virtual const regionCoupledBaseOFFBEATGAMGInterface& nbrPatch() const
            {
                return dynamic_cast<const regionCoupledBaseOFFBEATGAMGInterface&>
                (
                    fineregionCoupledOFFBEATLduInterface_.nbrPatch()
                );
            }

            //virtual const AMIPatchToPatchInterpolation& AMI() const
            //{
            //    return amiPtr_();
            //}

            virtual const polyMesh& nbrMesh() const
            {
                return fineregionCoupledOFFBEATLduInterface_.nbrMesh();
            }

            const lduInterface& nbrLduInterface
            (
                const label fineLevelIndex,
                const label index
            ) const
            {
                const GAMGAgglomeration& agg = nbrMesh().thisDb().lookupObject
                <
                    GAMGAgglomeration
                >
                (
                    GAMGAgglomeration::typeName
                );

                const lduMesh& mesh = agg.meshLevel(fineLevelIndex);
                return mesh.interfaces()[index];
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
