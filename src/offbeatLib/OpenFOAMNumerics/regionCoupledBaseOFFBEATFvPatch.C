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
    // Only owner updates AMI
    if (!owner())
    {
        return;
    }

    // External switch (e.g. solver-level control)
    if (!regionCoupledPatch().updateAMI())
    {
        return;
    }

    const fvMesh& mesh = patch_.boundaryMesh().mesh();

    // Moving frame → incremental field
    const bool movingFrame = mesh.foundObject<fvMesh>("referenceMesh");
    const word dispName = movingFrame ? "DD" : "D";
    const scalar timeValue = mesh.time().value();

    // Access patches
    const polyPatch& srcPatch =
        refCast<const polyPatch>(regionCoupledPatch());

    const polyPatch& tgtPatch =
        refCast<const polyPatch>(nbrPatch().regionCoupledPatch());

    // Displacement field
    const volVectorField& disp =
        mesh.lookupObject<volVectorField>(dispName);

    const fvPatchField<vector>& dispSrc =
        patch_.lookupPatchField<volVectorField, vector>(dispName);

    const vectorField dispTgt =
        regionCoupledPatch().interpolate
        (
            nbrPatch().patch().lookupPatchField<volVectorField, vector>(dispName)
        );

    // Relative motion across interface
    const vectorField slip = dispSrc - dispTgt;
    const scalarField slipNormal = slip & patch_.nf();
    const vectorField slipTangential = slip - patch_.nf()*slipNormal;

    // --- Initialize reference (first call)
    if (!slipReferenceInitialized_)
    {
        slipNormal0_ = slipNormal;
        slipTangential0_ = slipTangential;
        slipReferenceTimeValue_ = timeValue;
        slipReferenceInitialized_ = true;
        return;
    }

    // --- New time step in moving-frame mode → reset baseline
    if (movingFrame && timeValue != slipReferenceTimeValue_)
    {
        slipNormal0_ = slipNormal;
        slipTangential0_ = slipTangential;
        slipReferenceTimeValue_ = timeValue;
        return;
    }

    // --- Check drift since last AMI build
    const bool updateNormal =
        gMax(mag(slipNormal - slipNormal0_)) > regionCoupledPatch().updateAMINormalTol();

    const bool updateTangential =
        gMax(mag(slipTangential - slipTangential0_)) > regionCoupledPatch().updateAMITangentialTol();

    if (!(updateNormal || updateTangential))
    {
        return;
    }

    // --- Build displaced patches
    const volPointInterpolation& pointInterpolator =
        volPointInterpolation::New(mesh);

    const pointVectorField pointsDisp =
        pointInterpolator.interpolate(disp);

    pointField srcPoints(srcPatch.points() + pointsDisp.internalField());
    pointField tgtPoints(tgtPatch.points() + pointsDisp.internalField());

    twoDPointCorrector twoDCorrectorSrc(mesh);
    twoDCorrectorSrc.correctPoints(srcPoints);

    twoDPointCorrector twoDCorrectorTgt(mesh);
    twoDCorrectorTgt.correctPoints(tgtPoints);

    primitivePatch srcPatch0
    (
        SubList<face>(srcPatch, srcPatch.size(), 0),
        srcPoints
    );

    primitivePatch tgtPatch0
    (
        SubList<face>(tgtPatch, tgtPatch.size(), 0),
        tgtPoints
    );

    // --- Update AMI addressing
    AMIInterpolation& ami =
        const_cast<AMIInterpolation&>(regionCoupledPatch().AMI());

    Info<< "Updating AMI addressing on patch " << patch_.name() << nl;

#ifdef OPENFOAMFOUNDATION
    ami.update(srcPatch0, tgtPatch0, true);
#elif OPENFOAMESI
    ami.upToDate(false);
    ami.calculate(srcPatch0, tgtPatch0, regionCoupledBaseOFFBEAT_.surfPtr());
#endif

    // --- Refresh reference state
    slipNormal0_ = slipNormal;
    slipTangential0_ = slipTangential;
    slipReferenceTimeValue_ = timeValue;
}

// ************************************************************************* //
