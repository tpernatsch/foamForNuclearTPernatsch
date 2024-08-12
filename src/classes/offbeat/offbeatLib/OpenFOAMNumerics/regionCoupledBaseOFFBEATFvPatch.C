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

Description


\*---------------------------------------------------------------------------*/

#include "regionCoupledBaseOFFBEATFvPatch.H"
#include "volFields.H"

#ifdef OPENFOAMFOUNDATION
#include "AMIInterpolation.H"
#elif OPENFOAMESI
#include "AMIInterpolation.H"
#endif

#include "PrimitivePatchInterpolation.H"
#include "volPointInterpolation.H"
#include "twoDPointCorrector.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionCoupledBaseOFFBEATFvPatch, 0);
}


void Foam::regionCoupledBaseOFFBEATFvPatch::updateAMI() const
{
  // - Access next section only if owner
  // - Next section will update the AMI addressing
    if (owner() and regionCoupledPatch().updateAMI()) 
    {
        //- Take a reference to source and target polypatch
        const polyPatch& srcPatch
            = refCast<const polyPatch>(regionCoupledPatch());

        const polyPatch& tgtPatch
            = refCast<const polyPatch>(nbrPatch().regionCoupledPatch());

        //- Take a reference to displacement field and interpolate to mesh points
        const volVectorField disp = 
        patch_.boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh") ?
        patch_.boundaryMesh().mesh().lookupObject<volVectorField>("DD")
        :
        patch_.boundaryMesh().mesh().lookupObject<volVectorField>("D");

        const volPointInterpolation& meshPointInterpolation = volPointInterpolation::New( patch_.boundaryMesh().mesh());

        pointVectorField pointsDisplacement = meshPointInterpolation.interpolate(disp);

        pointField srcPoints = srcPatch.points() + pointsDisplacement.internalField();

        twoDPointCorrector twoDCorrectorSrc(patch_.boundaryMesh().mesh());
        twoDCorrectorSrc.correctPoints(srcPoints);

        //- Build source and target patch
        //- The two patches require a list of faces and the (displaced) mesh points.
        primitivePatch srcPatch0
        (
          SubList<face>
          (
              srcPatch,
              srcPatch.size(),
              0
          ),
          srcPoints
        );

        pointField tgtPoints = tgtPatch.points() + pointsDisplacement.internalField();

        twoDPointCorrector twoDCorrectorTgt(patch_.boundaryMesh().mesh());
        twoDCorrectorTgt.correctPoints(tgtPoints);

        primitivePatch tgtPatch0
        (
          SubList<face>
          (
              tgtPatch,
              tgtPatch.size(),
              0
          ),
          tgtPoints
        );

        //- Take a reference to AMI. 
        //- const_cast is necessary to access update() function
        AMIInterpolation& ami(const_cast<AMIInterpolation&>(regionCoupledPatch().AMI()));
        Info << "update AMI addressing" << endl;
#ifdef OPENFOAMFOUNDATION
        ami.update(srcPatch0, tgtPatch0, true);
#elif OPENFOAMESI  
        ami.upToDate(false);
        ami.calculate(srcPatch0, tgtPatch0, regionCoupledBaseOFFBEAT_.surfPtr());
        ami.tgtMagSf() = mag(tgtPatch0.faceAreas());
#endif        
    }
}

// ************************************************************************* //
