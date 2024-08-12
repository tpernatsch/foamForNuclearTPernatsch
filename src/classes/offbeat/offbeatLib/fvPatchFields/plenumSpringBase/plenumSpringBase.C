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

#include "plenumSpringBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(plenumSpringBase, 0);
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::plenumSpringBase::plenumSpringBase
(
    const fvPatch& p
)
:
    patch_(p),
    springModulus_(3.5e3),
    springPreCompression_(0.0),
    fuelTopPatchNames_(),
    fuelTopPatchIDs_(),
    topCapInnerPatchNames_(),
    topCapInnerPatchIDs_(),
    displacementName_("D")
{}


Foam::plenumSpringBase::plenumSpringBase
(
    const fvPatch& p,
    const dictionary& dict
)
:
    patch_(p),
    springModulus_(readScalar(dict.lookup("springModulus"))),
    springPreCompression_(readScalar(dict.lookup("springPreCompression"))),
    fuelTopPatchNames_(),
    fuelTopPatchIDs_(),
    topCapInnerPatchNames_(),
    topCapInnerPatchIDs_(),
    displacementName_(dict.lookupOrDefault<word>("displacementName", "D"))
{
    const fvMesh& mesh(patch_.boundaryMesh().mesh());

    // Find IDs of patches connected to bottom of the spring
    // TODO: add selection by patchType and allow the user to use all three 
    // selection modes (nameList, patchType and groups) by converting to set/HashTables
    if(dict.found("fuelTopPatches"))
    {
        wordList patchNames(dict.lookup("fuelTopPatches"));

        forAll(patchNames, i)
        {
            word patchName = patchNames[i];

            // Lookup ID of patch 
            label patchID = 
            mesh.boundaryMesh().findPatchID(patchName);

            if (patchID == -1)
            {
                FatalIOErrorInFunction(dict.lookup("fuelTopPatches"))
                    << "Specified patch " << patchName << " not found."
                    << abort(FatalIOError);
            }

            // Add patch name and ID to lists
            fuelTopPatchNames_.append(patchName);
            fuelTopPatchIDs_.append(patchID);
        }
    }
    else if(mesh.boundaryMesh().groupPatchIDs().find("fuelTop") 
        != mesh.boundaryMesh().groupPatchIDs().end())
    {
        labelList patchIDs(mesh.boundaryMesh().groupPatchIDs()["fuelTop"]);

        forAll(patchIDs, i)
        {
            label patchID = patchIDs[i];
            word patchName = mesh.boundaryMesh()[patchID].name();

            // Add patch name and ID to lists
            fuelTopPatchNames_.append(patchName);
            fuelTopPatchIDs_.append(patchID);
        }
    }

    if(!fuelTopPatchIDs_.size())
    {                       
        FatalErrorInFunction()
        <<"fuelTop patches not found when costructing " 
        << "plenumSpring fvPatchVectorField for patch \"" 
        << patch_.name() << "\"." << nl << nl 
        << "Either: " << nl
        << " - group the patches using the group name \"fuelTop\" or" << nl 
        << " - list the patch names under \"fuelTopPatches\" " 
        << "in the patch subDict"  << abort(FatalError);
    }

    // Find IDs of patches connected to top of the spring
    // TODO: add selection by patchType and allow the user to use all three 
    // selection modes (nameList, patchType and groups) by converting to set/HashTables
    if(dict.found("topCapInnerPatches"))
    {
        wordList patchNames(dict.lookup("topCapInnerPatches"));

        forAll(patchNames, i)
        {
            word patchName = patchNames[i];

            // Lookup ID of patch 
            label patchID = 
            mesh.boundaryMesh().findPatchID(patchName);

            if (patchID == -1)
            {
                FatalIOErrorInFunction(dict.lookup("topCapInnerPatches"))
                    << "Specified patch " << patchName << " not found."
                    << abort(FatalIOError);
            }

            // Add patch name and ID to lists
            topCapInnerPatchNames_.append(patchName);
            topCapInnerPatchIDs_.append(patchID);
        }
    }
    else if(mesh.boundaryMesh().groupPatchIDs().find("topCapInner") 
        != mesh.boundaryMesh().groupPatchIDs().end())
    {
        labelList patchIDs(mesh.boundaryMesh().groupPatchIDs()["topCapInner"]);

        forAll(patchIDs, i)
        {
            label patchID = patchIDs[i];
            word patchName = mesh.boundaryMesh()[patchID].name();

            // Add patch name and ID to lists
            topCapInnerPatchNames_.append(patchName);
            topCapInnerPatchIDs_.append(patchID);
        }
    }

    if(!topCapInnerPatchIDs_.size())
    {                       
        FatalErrorInFunction()
        <<"topCapInner patches not found when costructing " 
        << "plenumSpring fvPatchVectorField for patch \"" 
        << patch_.name() << "\"." << nl << nl 
        << "Either: " << nl
        << " - group the patches using the group name \"topCapInner\" or" << nl 
        << " - list the patch names under \"topCapInnerPatches\" " 
        << "in the patch subDict"  << abort(FatalError);
    }
}

Foam::plenumSpringBase::plenumSpringBase
(
    const fvPatch& p,
    const plenumSpringBase& psb
)
:
    patch_(p),
    springModulus_(psb.springModulus_),
    springPreCompression_(psb.springPreCompression_),
    fuelTopPatchNames_(psb.fuelTopPatchNames_),
    fuelTopPatchIDs_(psb.fuelTopPatchIDs_),
    topCapInnerPatchNames_(psb.topCapInnerPatchNames_),
    topCapInnerPatchIDs_(psb.topCapInnerPatchIDs_),
    displacementName_(psb.displacementName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::plenumSpringBase::~plenumSpringBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::plenumSpringBase::springElongation() const
{
    const fvMesh& mesh(patch_.boundaryMesh().mesh());

    const globalOptions& globalOpt
    (mesh.lookupObject<globalOptions>("globalOptions"));

    vector pinDirection = globalOpt.pinDirection();

    // Reference to total displacement field
    const volVectorField& D =
        mesh.lookupObject<volVectorField>(displacementName_); 

    // Displacement at the bottom of the spring
    scalar DavgBottom(0.0);
    scalar Abottom(0.0);

    forAll(fuelTopPatchIDs_, i)
    {
        label patchID(fuelTopPatchIDs_[i]);

        // Displacement for current patch
        const fvPatchField<vector>& springBottomD =
            D.boundaryField()[patchID];

        // Patch connected to the top of the spring
        const fvPatch& springBottomPatch =
            mesh.boundary()[patchID];

        DavgBottom += 
        gSum((springBottomD & pinDirection)*springBottomPatch.magSf());

        Abottom += gSum(springBottomPatch.magSf());
    }

    DavgBottom = DavgBottom/max(Abottom, VSMALL);

    // Displacement at the top of the spring
    scalar DavgTop(0.0);
    scalar Atop(0.0);

    forAll(topCapInnerPatchIDs_, i)
    {
        label patchID(topCapInnerPatchIDs_[i]);

        // Displacement for current patch
        const fvPatchField<vector>& springTopD =
            D.boundaryField()[patchID];

        // Patch connected to the top of the spring
        const fvPatch& springTopPatch =
            mesh.boundary()[patchID];

        DavgTop += 
        gSum((springTopD & pinDirection)*springTopPatch.magSf());

        Atop += gSum(springTopPatch.magSf());
    }

    DavgTop = DavgTop/max(Atop, VSMALL);

    // Total spring elongation
    scalar Davg = DavgBottom - DavgTop;

    // Calculate pressure exterded by spring
    // Pressure cannot be negative (or spring force cannot be positive)
    scalar Dtot = 
    max(Davg + springPreCompression_, 0.0);

	return Dtot;
}

void Foam::plenumSpringBase::write(Ostream& os) const
{
#ifdef OPENFOAMFOUNDATION    
    writeEntry<scalar>(os, "springModulus", springModulus_);
    writeEntry<scalar>(os, "springPreCompression", springPreCompression_);
    writeEntry<wordList>(os, "fuelTopPatches", fuelTopPatchNames_);
    writeEntry<wordList>(os, "topCapInnerPatches", topCapInnerPatchNames_);
    writeEntryIfDifferent<word>(os, "displacementName", "D", displacementName_);
#elif OPENFOAMESI
    os.writeEntry<scalar>("springModulus", springModulus_);
    os.writeEntry<scalar>("springPreCompression", springPreCompression_);
    os.writeEntry<wordList>("fuelTopPatches", fuelTopPatchNames_);
    os.writeEntry<wordList>("topCapInnerPatches", topCapInnerPatchNames_);
    os.writeEntryIfDifferent<word>("displacementName", "D", displacementName_);
#endif
}

// ************************************************************************* //
