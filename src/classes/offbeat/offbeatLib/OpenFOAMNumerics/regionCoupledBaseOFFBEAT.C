/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "regionCoupledBaseOFFBEAT.H"
#include "SubField.H"
#include "polyMesh.H"
#include "Time.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionCoupledBaseOFFBEAT, 0);
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::regionCoupledBaseOFFBEAT::resetAMI() const
{
    if (owner())
    {
        AMIPtr_.clear();

        const polyPatch& nbr = refCast<const polyPatch>(nbrPatch());
        pointField nbrPoints = nbr.localPoints();

        if (debug)
        {
            const Time& t = patch_.boundaryMesh().mesh().time();
            OFstream os(t.path()/patch_.name() + "_neighbourPatch-org.obj");
            meshTools::writeOBJ(os, nbr.localFaces(), nbrPoints);
        }

        // transform neighbour patch to local system
        //transformPosition(nbrPoints);
        primitivePatch nbrPatch0
        (
            SubList<face>
            (
                nbr.localFaces(),
                nbr.size()
            ),
            nbrPoints
        );

        if (debug)
        {
            const Time& t = patch_.boundaryMesh().mesh().time();
            OFstream osN(t.path()/patch_.name() + "_neighbourPatch-trans.obj");
            meshTools::writeOBJ(osN, nbrPatch0.localFaces(), nbrPoints);

            OFstream osO(t.path()/patch_.name() + "_ownerPatch.obj");
            meshTools::writeOBJ
            (
                osO,
                patch_.localFaces(),
                patch_.localPoints()
            );
        }
#ifdef OPENFOAMFOUNDATION        
        // Construct/apply AMI interpolation to determine addressing and weights
        AMIPtr_.reset
        (
         new AMIInterpolation
         (
             patch_,
             nbrPatch0,
             surfPtr(),
             faceAreaIntersect::tmMesh,
             true,
             AMIInterpolation::wordTointerpolationMethod(methodName_),
             -1,
             AMIReverse_
         )
        );
#elif OPENFOAMESI        
        // Construct/apply AMI interpolation to determine addressing and weights
        AMIPtr_.reset
        (
            AMIInterpolation::New
            (
                methodName_,
                true,
                AMIReverse_,
                -1        
            )
        );  

        AMIPtr_->calculate(patch_, nbrPatch0, surfPtr());
        AMIPtr_->tgtMagSf() = mag(nbrPatch0.faceAreas());        
#endif

        if (debug)
        {
            Pout<< "regionCoupledBaseOFFBEAT : " << patch_.name()
                << " constructed AMI with " << nl
                << "    " << ":srcAddress:" << AMIPtr_().srcAddress().size()
                << nl
                << "    " << " tgAddress :" << AMIPtr_().tgtAddress().size()
                << nl << endl;
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::regionCoupledBaseOFFBEAT::clearGeom()
{
    AMIPtr_.clear();
    surfPtr_.clear();
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::regionCoupledBaseOFFBEAT::regionCoupledBaseOFFBEAT
(
    const polyPatch& pp
)
:
    patch_(pp),
    owner_(false),
    nbrPatchName_(word::null),
    nbrPatchID_(-1),
    nbrRegionName_(word::null),
    sameRegion_(false),
    AMIPtr_(nullptr),
    AMIReverse_(false),
    surfPtr_(nullptr),
    surfDict_(fileName("surface")),
    updateAMI_(false),
    methodName_("faceAreaWeightAMI")
{}


Foam::regionCoupledBaseOFFBEAT::regionCoupledBaseOFFBEAT
(
    const polyPatch& pp,
    const dictionary& dict
)
:
    patch_(pp),
    owner_(dict.lookupOrDefault<bool>("owner", false)),
    nbrPatchName_(dict.lookup("neighbourPatch")),
    nbrPatchID_(-1),
    nbrRegionName_(dict.lookup("neighbourRegion")),
    sameRegion_(nbrRegionName_ == patch_.boundaryMesh().mesh().name()),
    AMIPtr_(nullptr),
    AMIReverse_(dict.lookupOrDefault<bool>("flipNormals", false)),
    surfPtr_(nullptr),
    surfDict_(dict.subOrEmptyDict("surface")),
    updateAMI_(dict.lookupOrDefault<bool>("updateAMI", false)),
    methodName_(dict.lookupOrDefault<word>("AMIMethod", "faceAreaWeightAMI"))
{}


Foam::regionCoupledBaseOFFBEAT::regionCoupledBaseOFFBEAT
(
    const polyPatch& pp,
    const regionCoupledBaseOFFBEAT& mpb
)
:
    patch_(pp),
    owner_(mpb.owner_),
    nbrPatchName_(mpb.nbrPatchName_),
    nbrPatchID_(mpb.nbrPatchID_),
    nbrRegionName_(mpb.nbrRegionName_),
    sameRegion_(mpb.sameRegion_),
    AMIPtr_(nullptr),
    AMIReverse_(mpb.AMIReverse_),
    surfPtr_(mpb.surfPtr_),
    surfDict_(mpb.surfDict_),
    updateAMI_(mpb.updateAMI_),
    methodName_("faceAreaWeightAMI")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionCoupledBaseOFFBEAT::~regionCoupledBaseOFFBEAT()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::regionCoupledBaseOFFBEAT::nbrPatchID() const
{
    if (nbrPatchID_ == -1)
    {
        if
        (
            patch_.boundaryMesh().mesh().time().foundObject<polyMesh>
            (
                nbrRegionName_
            )
        )
        {
            const polyMesh& mesh =
                patch_.boundaryMesh().mesh().time().lookupObject<polyMesh>
                (
                    nbrRegionName_
                );

            nbrPatchID_ = mesh.boundaryMesh().findPatchID(nbrPatchName_);

            if (nbrPatchID_ == -1)
            {
                FatalErrorInFunction
                    << "Illegal neighbourPatch name " << nbrPatchName_
                    << nl << "Valid patch names are "
                    << mesh.boundaryMesh().names()
                    << exit(FatalError);
            }

            // Check that it is a cyclic AMI patch
            const regionCoupledBaseOFFBEAT& nbrPatch =
                refCast<const regionCoupledBaseOFFBEAT>
                (
                    mesh.boundaryMesh()[nbrPatchID_]
                );

            if (nbrPatch.nbrPatchName() != patch_.name())
            {
                WarningInFunction
                    << "Patch " << patch_.name()
                    << " specifies neighbour patch " << nbrPatchName()
                    << nl << " but that in return specifies "
                    << nbrPatch.nbrPatchName() << endl;
            }
        }
    }

    return nbrPatchID_;
}


bool Foam::regionCoupledBaseOFFBEAT::owner() const
{
    if(!owner_ and !nbrPatch().owner_)
    {
        if (nbrRegionName_ == patch_.boundaryMesh().mesh().name())
        {
            return patch_.index() < nbrPatchID();
        }
        else
        {
            return patch_.boundaryMesh().mesh().name() < nbrRegionName_;
        }
    }

    return owner_;
}


const Foam::autoPtr<Foam::searchableSurface>&Foam::regionCoupledBaseOFFBEAT::
surfPtr() const
{
    const word surfType(surfDict_.lookupOrDefault<word>("type", "none"));

    if (!surfPtr_.valid() && owner() && surfType != "none")
    {
        word surfName(surfDict_.lookupOrDefault("name", patch_.name()));

        const polyMesh& mesh = patch_.boundaryMesh().mesh();

        surfPtr_ =
            searchableSurface::New
            (
                surfType,
                IOobject
                (
                    surfName,
                    mesh.time().constant(),
                    "triSurface",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                surfDict_
            );
    }

    return surfPtr_;
}


const Foam::AMIInterpolation& Foam::regionCoupledBaseOFFBEAT::AMI() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "AMI interpolator only available to owner patch"
            << abort(FatalError);
    }

    if (!AMIPtr_.valid())
    {
        resetAMI();
    }

    return AMIPtr_();
}


const Foam::regionCoupledBaseOFFBEAT&
Foam::regionCoupledBaseOFFBEAT::nbrPatch() const
{
    const polyMesh& mesh =
        patch_.boundaryMesh().mesh().time().lookupObject<polyMesh>
        (
            nbrRegionName_
        );

    const polyPatch& pp = mesh.boundaryMesh()[nbrPatchID()];
    return refCast<const regionCoupledBaseOFFBEAT>(pp);
}


bool Foam::regionCoupledBaseOFFBEAT::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    return false;
}


void Foam::regionCoupledBaseOFFBEAT::write(Ostream& os) const
{
    os.writeKeyword("neighbourPatch") << nbrPatchName_
    << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourRegion") << nbrRegionName_
    << token::END_STATEMENT << nl;

    if (AMIReverse_)
    {
        os.writeKeyword("flipNormals") << AMIReverse_
            << token::END_STATEMENT << nl;
    }

    if (!surfDict_.empty())
    {
        os.writeKeyword(surfDict_.dictName());
        os  << surfDict_;
    }
}
// ************************************************************************* //
