/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "gapFRAPCON.H"
#include "addToRunTimeSelectionTable.H"
#include "physicoChemicalConstants.H"

#include "regionCoupledOFFBEATFvPatch.H"

#include "PrimitivePatchInterpolation.H"
#include "volPointInterpolation.H"
#include "cuttingPlane.H"
#include "twoDPointCorrector.H"
#include "fuelMaterial.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gapFRAPCON, 0);
    addToRunTimeSelectionTable
    (
        gapGasModel, 
        gapFRAPCON, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::gapFRAPCON::calcInitialMass()
{
    // Total volume available to the gas
    const scalar Vtot = 
    (gapV_ + gapOffsetV_) + (holeV_) + (dishV_) + (topPlenumV_ + reserveV_) +
    (bottomPlenumV_) + crackV_;

    // Volume weighted gas temperature in the rod
    const scalar gasT = 
    (
        gapV_*gapT_ + gapOffsetV_*gapT_ + holeV_*holeT_ + dishV_*dishT_ +
        topPlenumV_*topPlenumT_ + reserveV_*reserveT_ + bottomPlenumV_*bottomPlenumT_ 
        + crackV_*crackT_
    )/Vtot;

    // Gas density using ideal gas law
    scalar rho = rhoMixture(gasP_, gasT);
    
    // Gas initial mass, in kg
    gasM0_ = rho*Vtot ;

    // Update gas mass
    gasM_ = gasM0_;
}


Foam::scalar Foam::gapFRAPCON::rhoMixture
(
    const scalar& gasP, 
    const scalar& gasT
) const
{
    scalar rhoMixture(0.0);
    scalar R_mixture(0.0);
    scalar R = Foam::constant::physicoChemical::R.value();

    forAll(species_, i)
    {
        R_mixture += R/speciesW_[i]*1000*Y_[i];
    }
    
    rhoMixture = gasP/(R_mixture*gasT);

    return rhoMixture;
}


void Foam::gapFRAPCON::findPatchIDs(const word groupName)
{
    if(!patchIDsPtr_.valid())
    {
        patchIDsPtr_.reset(new HashTable<std::set<label>, word>(10));
    }

    HashTable<std::set<label>, word>& patchIDs = patchIDsPtr_();

    // Find associated patches on the mesh that bound the gap gas
    if
    (
#ifdef OPENFOAMFOUNDATION        
        !gapGasOptDict_.isNull() 
#elif OPENFOAMESI
        !gapGasOptDict_.isNullDict() 
#endif        
        and gapGasOptDict_.found(groupName+"Patches")
    )
    {
        wordList gapPatches(gapGasOptDict_.lookup(groupName+"Patches"));

        forAll(gapPatches, patchI)
        {
            word gapPatchName(gapPatches[patchI]);

            label gapPatchID
            (mesh_.boundaryMesh().findPatchID(gapPatchName));

            if(gapPatchID == -1)
            {
                FatalErrorIn("Foam::gapFRAPCON::findPatchIDs(const word groupName)")
                << gapPatchName << " patch from list " << groupName << "Patches"
                << " not found in mesh." << exit(FatalError);
            }

            HashTable<std::set<label>, word>::iterator iter = patchIDs.find
            (
                groupName
            );

            if (iter != patchIDs.end())
            {
                iter().insert(gapPatchID);
            }
            else
            {
                std::set<label> s;
                s.insert(gapPatchID);
                patchIDs.insert(groupName, s);
            }
        }

        if(!gapPatches.size())
        {
            // Insert an empty set --> allow an empty patch list
            std::set<label> s;
            patchIDs.insert(groupName, s);
        }
    }    
    
    const Foam::HashTable<Foam::labelList, Foam::word>& groupPatchIDs =
    mesh_.boundaryMesh().groupPatchIDs();

    if(groupPatchIDs.find(groupName) != groupPatchIDs.end())
    {  
        labelList gapPatchIDs(mesh_.boundaryMesh().groupPatchIDs()[groupName]);

        forAll(gapPatchIDs, patchI)
        {
            label gapPatchID(gapPatchIDs[patchI]);

            HashTable<std::set<label>, word>::iterator iter = patchIDs.find
            (
                groupName
            );

            if (iter != patchIDs.end())
            {
                iter().insert(gapPatchID);
            }
            else
            {
                std::set<label> s;
                s.insert(gapPatchID);
                patchIDs.insert(groupName, s);
            }
        }
    }  

    if
    (
        patchIDs.find(groupName) == patchIDs.end()
    )
    {                       
        FatalErrorIn("Foam::gapFRAPCON::findPatchIDs(const word groupName)")
        << groupName << " patches not found." << nl 
        << "Either: " << nl
        << " - group the patches using the group name \""
        << groupName << "\" or" << nl 
        << " - list the patch names in \"" << groupName << "Patches\" " 
        << "in the gapGasOptions dict in solverDict file"  << nl 
        << " (the list can be empty)"  << nl ;

        if(groupName == "hole")
        {
            FatalError() << " - switch off the keyword \"includeCentralHole\" "
            << "in gapGasOptions dict"<< exit(FatalError);
        }

        Info << exit(FatalError);
    }
}


void Foam::gapFRAPCON::correctScalingFactors()
{
    // Create two cutting planes corresonding to zMin and zMax (the normal is 
    //  the pinDirection_). Find cut faces on the cladding surface. Intersect
    //  each edge of the cut face with the cutting plane. With the intersection 
    //  points build a face and calculate the area. The ratio between this area
    //  and the total face area is used to obtain the scaling factor for that 
    //  face. Above and below the cutting plane, set the scaling factor to 0.

    // TODO check if pinDirection as normal is correct (what if top patch
    //  rotates? Probably it is negligible.)

    // Move pointField to the latest position. Necessary to correctly intersect
    //  planes and edges

    // Constant reference to class for mesh point interpolation
    const volPointInterpolation& meshPointInterpolation = 
        volPointInterpolation::New(mesh_);

    // Compute the displacement of the points (due to deltaD)
    pointVectorField pointsDisplacement = 
        meshPointInterpolation.interpolate(DorDD_);

    // Calculate the new position of the displaced mesh points
    vectorField pointsDisplaced = 
        mesh_.points()
    +   pointsDisplacement.internalField();

    twoDPointCorrector twoDCorrector(mesh_);
    twoDCorrector.correctPoints(pointsDisplaced);

    // Take reference to gap patch IDs list
    HashTable<std::set<label>, word>& patchIDs = patchIDsPtr_();
    std::set<label>& gapPatchIDs(patchIDs["gap"]);

    // Loop over all gap patches
    for(const label& patchID : gapPatchIDs) 
    {
        // const label patchID(gapPatchIDs[patchI]);
        const fvPatch& p = mesh_.boundary()[patchID];
        
        const regionCoupledOFFBEATFvPatch& patch
            = refCast<const regionCoupledOFFBEATFvPatch>(p);    
        
        // NOTE: fuel is supposed to be the owner
        if(!patch.owner())
        {
            const fvPatchVectorField& Df = 
            DorDD_.boundaryField()[patchID];

            scalarField& factors_gap(scalingFactors_[patchID]);
            scalarField& factors_top(scalingFactorsTopPlenum_[patchID]);
            scalarField& factors_bottom(scalingFactorsBottomPlenum_[patchID]);

            // Above and below zMax and zMin, set scaling factors to 0
            factors_gap = 1;
            factors_top = 0;
            factors_bottom = 0;

            // Create top and bottom cutting planes
            vector planePointTop = zMax_*pinDirection_;
            vector planePointBottom = zMin_*pinDirection_;

            vector planeNormalTop = pinDirection_/mag(pinDirection_);
            vector planeNormalBottom = -pinDirection_/mag(pinDirection_);

            plane plTop(planePointTop, planeNormalTop);
            plane plBottom(planePointBottom, planeNormalBottom);

            // Loop over every face
            forAll(patch, faceID)
            {
                scalar axialLocation
                (
                    (patch.Cf()[faceID] + Df[faceID]) & pinDirection_
                );

                if(axialLocation > zMax_)
                {
                    factors_gap[faceID] = 0;
                    factors_top[faceID] = 1;
                }
                if(axialLocation < zMin_)
                {
                    factors_gap[faceID] = 0;
                    factors_bottom[faceID] = 1;
                }

                scalar maxAxialLocation = 
                max(p.patch()[faceID].points(pointsDisplaced) & pinDirection_);
                scalar minAxialLocation = 
                min(p.patch()[faceID].points(pointsDisplaced) & pinDirection_);

                // Find if face is cut by either of the two planes
                bool cutByTopPlane = false;
                bool cutByBottomPlane = false;

                if(maxAxialLocation > zMax_ and minAxialLocation < zMax_)
                {
                    cutByTopPlane = true;
                }

                if(maxAxialLocation > zMin_ and minAxialLocation < zMin_)
                {
                    cutByBottomPlane = true;
                }
    
                if(cutByTopPlane || cutByBottomPlane)
                {
                    const edgeList edges(p.patch()[faceID].edges());
                    
                    // Prepare pointfield and label list to construct face
                    List<point> topPlenumFacePts;
                    List<point> bottomPlenumFacePts;
                    List<point> gapFacePts;
                    List<point> totalFacePts;

                    labelList topPlenumFaceLabels;
                    labelList bottomPlenumFaceLabels;
                    labelList gapFaceLabels;
                    labelList totalFaceLabels;

                    // Loop over all edges
                    forAll(edges, edgeID)
                    {
                        // Select edge
                        const edge edgei(edges[edgeID]);

                        // Edge base point (from the displaced point list)
                        const point base(pointsDisplaced[edgei.start()]);

                        // Add base point to totalFace point list
                        totalFacePts.append(base);
                        totalFaceLabels.append
                        (totalFaceLabels.size());

                        // Edge vector (not normalized, actually it is the 
                        // vectorial distance between the two points of the edge)
                        const vector dir(edgei.vec(pointsDisplaced));

                        // Distance between base and intersection with the two planes (if any)
                        // NOTE: the distance is a fraction of the vector dir
                        scalar distToTop(plTop.normalIntersect(base, dir));
                        scalar distToBottom(plBottom.normalIntersect(base, dir));

                        if
                        (
                            // If point is below bottom plane 
                            (plBottom.sideOfPlane(base) != 
                            Foam::plane::side::FLIP) 
                        )
                        {
                            // Append to bottom list
                            bottomPlenumFacePts.append(base);
                            bottomPlenumFaceLabels.append
                            (bottomPlenumFaceLabels.size());

                            // If intersection with bottom plane is within edge lenght
                            if (distToBottom <= 1.0 & distToBottom >= 0)
                            {
                                // Add point to bottom list
                                bottomPlenumFacePts.append(base + distToBottom*dir);
                                bottomPlenumFaceLabels.append
                                (bottomPlenumFaceLabels.size());
                                
                                // Add point to gap list
                                gapFacePts.append(base + distToBottom*dir);
                                gapFaceLabels.append
                                (gapFaceLabels.size());
                            } 

                            // If intersection with top plane is withing edge lenght
                            if (distToTop <= 1.0 & distToTop >= 0)
                            {
                                // Append to gap list
                                gapFacePts.append(base + distToTop*dir);
                                gapFaceLabels.append
                                (gapFaceLabels.size());

                                // Append to top list 
                                topPlenumFacePts.append(base + distToTop*dir);
                                topPlenumFaceLabels.append
                                (topPlenumFaceLabels.size());
                            }
                        }
                        else if
                        (
                            // If point is above top plane 
                            (plTop.sideOfPlane(base) != 
                            Foam::plane::side::FLIP)
                        )
                        {

                            // Append to top list 
                            topPlenumFacePts.append(base);
                            topPlenumFaceLabels.append
                            (topPlenumFaceLabels.size());

                            // If intersection with top plane is within edge lenght
                            if (distToTop <= 1.0 & distToTop > VSMALL)
                            {  
                                // Append to top list
                                topPlenumFacePts.append(base + distToTop*dir);
                                topPlenumFaceLabels.append
                                (topPlenumFaceLabels.size());
                                
                                // Append to gap list   
                                gapFacePts.append(base + distToTop*dir);
                                gapFaceLabels.append
                                (gapFaceLabels.size());
                            } 

                            // If intersection with bottom plane is within edge lenght
                            if (distToBottom <= 1.0 & distToBottom > VSMALL)
                            {          
                                // Append to gap list
                                gapFacePts.append(base + distToBottom*dir);
                                gapFaceLabels.append
                                (gapFaceLabels.size());

                                // Append to bottom list  
                                bottomPlenumFacePts.append(base + distToBottom*dir);
                                bottomPlenumFaceLabels.append
                                (bottomPlenumFaceLabels.size());
                            }
                        }
                        else
                        {

                            // Append to gap list 
                            gapFacePts.append(base);
                            gapFaceLabels.append
                            (gapFaceLabels.size());

                            // If intersection with top plane is within edge lenght
                            if (distToTop <= 1.0 & distToTop > VSMALL)
                            {     
                                // Append to gap list
                                gapFacePts.append(base + distToTop*dir);
                                gapFaceLabels.append
                                (gapFaceLabels.size());
                                
                                // Append to top list 
                                topPlenumFacePts.append(base + distToTop*dir);
                                topPlenumFaceLabels.append
                                (topPlenumFaceLabels.size());
                            }

                            // If intersection with bottom plane is within edge lenght
                            if (distToBottom <= 1.0 & distToBottom > VSMALL)
                            {        
                                // Append to gap list
                                gapFacePts.append(base + distToBottom*dir);
                                gapFaceLabels.append
                                (gapFaceLabels.size());

                                bottomPlenumFacePts.append(base + distToBottom*dir);
                                bottomPlenumFaceLabels.append
                                (bottomPlenumFaceLabels.size());                                    
                            }
                        }
                    }

                    pointField overlapPtsField(gapFacePts);
                    face overlapFace(gapFaceLabels);
                    scalar overlapArea
#ifdef OPENFOAMFOUNDATION                        
                    (mag(overlapFace.area(overlapPtsField)));
#elif OPENFOAMESI
                    (mag(overlapFace.areaNormal(overlapPtsField)));
#endif

                    pointField facePtsField(totalFacePts);
                    face totalFace(totalFaceLabels);
                    scalar totalArea
#ifdef OPENFOAMFOUNDATION                        
                    (mag(totalFace.area(facePtsField)));
#elif OPENFOAMESI
                    (mag(totalFace.areaNormal(facePtsField)));
#endif

                    pointField topPlenumPtsField(topPlenumFacePts);
                    face topPlenumFace(topPlenumFaceLabels);
                    scalar topPlenumArea
#ifdef OPENFOAMFOUNDATION                        
                    (mag(topPlenumFace.area(topPlenumPtsField)));
#elif OPENFOAMESI
                    (mag(topPlenumFace.areaNormal(topPlenumPtsField)));
#endif

                    pointField bottomPlenumPtsField(bottomPlenumFacePts);
                    face bottomPlenumFace(bottomPlenumFaceLabels);
                    scalar bottomPlenumArea
#ifdef OPENFOAMFOUNDATION                        
                    (mag(bottomPlenumFace.area(bottomPlenumPtsField)));
#elif OPENFOAMESI
                    (mag(bottomPlenumFace.areaNormal(bottomPlenumPtsField)));
#endif

                    factors_gap[faceID] = overlapArea/totalArea;
                    factors_top[faceID] = topPlenumArea/totalArea;
                    factors_bottom[faceID] = bottomPlenumArea/totalArea;
                }
            }
        }
    }
}


Foam::scalar Foam::gapFRAPCON::patchListAvgAxialLocation
(const std::set<label>& patchIDs)
{
    // Initialize position vectors and areas
    vector posVector(Foam::vector(0,0,0));
    scalar areaTot(0.0);

    // Loop over and update posVector and areaTot
    for(const label& patchID : patchIDs)
    {
        // Necessary to access polyPatch as the fvPatch might be null in case 
        //  of empty constraint
        const polyPatch& patch = 
        mesh_.boundary()[patchID].patch(); 

        const fvPatchVectorField& Dp = 
        DorDD_.boundaryField()[patchID];

        const vectorField& Cf = 
        patch.faceCentres();

        const scalarField& magSf = 
        patch.magFaceAreas();

        if(Dp.size())
        {
            // If the field is of non-empty type
            forAll(patch, faceI)
            {
                posVector += 
                (Cf[faceI] + Dp[faceI])*magSf[faceI];
                areaTot += magSf[faceI];
            }  
        }
        else
        { 
            // If the field is of empty type take the internal field values
            labelList fCells = patch.faceCells();

            forAll(patch, faceI)
            {
                label cellI(fCells[faceI]);

                posVector += 
                (Cf[faceI] + Dp.internalField()[cellI])*magSf[faceI];
                areaTot += magSf[faceI];
            }         
        }
    }

    // Handle parallelisation
    reduce(posVector, sumOp<vector>());
    reduce(areaTot, sumOp<scalar>());
    posVector /= max(areaTot,VSMALL);

    return posVector & pinDirection_;
}


void Foam::gapFRAPCON::correctGap()
{  
    // Initialize useful quantities
    scalar volume(0);
    scalar area(0);
    scalar temperature(0.0);
    scalar VoverT(0.0);
        
    // Loop over specified patches to calculate total gas volume from the
    // geometry. If the cladding surface is larger than the fuel surface, the 
    // calculated volume will include the plena. 
    HashTable<std::set<label>, word>& patchIDs = patchIDsPtr_();
    std::set<label>& gapPatchIDs(patchIDs["gap"]);

    for(const label& patchID : gapPatchIDs) 
    {   
        const fvPatch& p = mesh_.boundary()[patchID];  

        const regionCoupledOFFBEATFvPatch& patch
            = refCast<const regionCoupledOFFBEATFvPatch>(p);
        const regionCoupledOFFBEATFvPatch& nbrPatch
            = refCast<const regionCoupledOFFBEATFvPatch>(patch.nbrPatch());  
 
        // Current patch quantities

        // Patch center face vectors, displacement and temperature fields
        const vectorField& cf = patch.Cf();
        const fvPatchVectorField& Df = DorDD_.boundaryField()[patchID];
        const fvPatchScalarField& T = T_.boundaryField()[patchID];

        // Patch normal vectors
        // TODO: apply F to rotate normal in case of large strain
        const vectorField nf = patch.nf();        
        
        // Take face areas from AMI as it is more precise than the patch
        // (only if updateAMI is true, otherwise it gives the same value)
        const AMIInterpolation& ami = patch.owner() ?
        patch.regionCoupledPatch().AMI()
        :
        patch.nbrPatch().regionCoupledPatch().AMI();       

        List<scalar> magSf = patch.owner() ? 
        ami.srcMagSf() : ami.tgtMagSf();


        // Scale the area with surface overlap scaling factors
        scalarField factors = scalingFactors_[patchID];
        magSf = magSf*factors;

        // Neihgbor patch quantities

        const word nbrPatchName = nbrPatch.name();
        const label nbrPatchID = 
        mesh_.boundaryMesh().findPatchID(nbrPatchName);

        const vectorField nbrCf = 
        patch.regionCoupledPatch().interpolate(nbrPatch.Cf());

        const vectorField nbrDf =  
        patch.regionCoupledPatch().interpolate
        (
            nbrPatch.lookupPatchField<volVectorField, vector>(DorDD_.name())
        );

        const vectorField nbrNf = 
        patch.regionCoupledPatch().interpolate(nbrPatch.nf());

        // The nbr face areas are obtained by scaling the owner areas by the 
        //  fraction of target face covered by source faces 

        // This is due to the fact that there is no simple way to obtain the 
        // corresponding nbrMagSf from the nbr side. Siply nbrMagSf*factors
        // would not work when interpolated back to the owner face in case of
        // different meshes 
        scalarField nbrScaleFactor = patch.owner() ?
        patch.regionCoupledPatch().interpolate
        (
            ami.tgtWeightsSum()/max(scalingFactors_[nbrPatchID], VSMALL)
        )
        :
        scalingFactors_[patchID]/max(ami.tgtWeightsSum(), VSMALL);

        scalarField nbrMagSf = magSf/max(nbrScaleFactor, VSMALL);

        // The pos(gapwidth) is used to scale the volumes so that 
        // there is no gap volume in case of gap closure
        scalarField gapWidth( ((nbrCf + nbrDf) - (cf + Df)) & nf);

        // A closed volume can be calculated as:
        // 1/3*surfIntegral(S*n dS)
        // A correction factor of 1.5 is needed as the gap volume is not 
        // bounded axially
        scalarField ownVolHole = 
        -(1/3.0)*(1.5*((cf + Df) & nf * magSf))*pos(gapWidth);

        scalarField nbrVolHole = 
        -(1/3.0)*(1.5*((nbrCf + nbrDf) & nbrNf * nbrMagSf))*pos(gapWidth);

        // The volume is calculated as the difference between the volume inside
        //  the current patch minus the volume inside the nbr patch.
        scalarField volHole =  max(nbrVolHole + ownVolHole, 0.0);

        // As the previous volume is calculated on both sides of the gap
        //  we need to divide by 2 when summing
        volume += 0.5*gSum(volHole); 
        area += gSum(magSf);

        temperature += gSum(T*magSf);
        VoverT += 0.5*gSum(1/T*volHole);
    }

    // Remove plena volumes 
    gapV_ = volume;
    gapT_ = temperature/max(area, VSMALL);
    gapVoverT_ = VoverT;
}


void Foam::gapFRAPCON::correctHole()
{  
    // Initialize useful quantities
    scalar volume(0);
    scalar temperature(0.0);
    scalar VoverT(0.0);
        
    // Loop over hole patches to calculate hole volume, temperature and VoverT.
    HashTable<std::set<label>, word>& patchIDs = patchIDsPtr_();
    std::set<label>& holePatchIDs(patchIDs["hole"]);

    for(const label& patchID : holePatchIDs) 
    {   
        const fvPatch& p = mesh_.boundary()[patchID];   

        const fvPatchVectorField& dispP = 
        DorDD_.boundaryField()[patchID];

        const fvPatchScalarField& T =
        T_.boundaryField()[patchID];

        // TODO apply F to rotate face normal vectors

        // A closed volume can be calculated as:
        // 1/3*surfIntegral(S*n dS)
        // A correction factor of 1.5 is needed as the hole volume is not 
        // bounded axially
        scalarField volHole = -(1/3.0)*(1.5*((p.Cf() + dispP) & p.Sf()));

        volume += gSum(volHole); 
        temperature += gSum(T*volHole);
        VoverT += gSum(1/T*volHole);
    }

    // Add offset 
    holeV_ = volume;
    holeT_ = temperature/max(volume, VSMALL);
    holeVoverT_ = VoverT;
}


void Foam::gapFRAPCON::correctPlena()
{        
    // Initialize plena volume, temperature and area
    scalar volumeTop(0);
    scalar areaTop(0);
    scalar temperatureTop(0);   

    scalar volumeBottom(0);
    scalar areaBottom(0);
    scalar temperatureBottom(0);    

    // Loop over gap patches to calculate plena temperature, volume and area
    HashTable<std::set<label>, word>& patchIDs = patchIDsPtr_();
    std::set<label>& gapPatchIDs(patchIDs["gap"]);

    for(const label& patchID : gapPatchIDs) 
    {
        const fvPatch& p = mesh_.boundary()[patchID];      

        const regionCoupledOFFBEATFvPatch& patch
            = refCast<const regionCoupledOFFBEATFvPatch>(p); 

        // Patch center face vectors, displacement and temperature fields
        const vectorField& cf = patch.Cf();
        const fvPatchVectorField& Df = DorDD_.boundaryField()[patchID];
        const fvPatchScalarField& T = T_.boundaryField()[patchID];

        // Patch normal vectors
        // TODO: apply F to rotate normal in case of large strain
        vectorField nf = patch.nf();        
        
        // Take face areas from AMI as it is more precise than the patch
        // (only if updateAMI is true, otherwise it gives the same value)
        const AMIInterpolation& ami = patch.owner() ?
        patch.regionCoupledPatch().AMI()
        :
        patch.nbrPatch().regionCoupledPatch().AMI();          

        List<scalar> magSf = patch.owner() ? 
        ami.srcMagSf() : ami.tgtMagSf();


        // // Pout<<"magSf" <<magSf<<endl;
        // // Info<<"tot magSf" <<magSf<<endl;

        // if(patch.owner())
        //     ami.srcMap().distribute(magSf);
        // else    
        //     ami.tgtMap().distribute(magSf);


        const scalarField factors = scalingFactors_[patchID];
        const scalarField factorsTop = scalingFactorsTopPlenum_[patchID];
        const scalarField factorsBottom = scalingFactorsBottomPlenum_[patchID];

        // Axial locations
        scalarField axialLocation = cf & pinDirection_;

        // NOTE: the fuel surface should have factors equal to 1.0, so that
        //  it is not taken into account

        // Top plenum quantities
        scalarField volTop = 
        -1/3.0*(1.5*((cf + Df) & nf * magSf*factorsTop));

        scalarField tempTop = T*magSf*factorsTop;
        scalarField SfTop = magSf*factorsTop;

        // Bottom plenum quantities

        scalarField volBottom = 
        -1/3.0*(1.5*((cf + Df) & nf * magSf*factorsBottom));

        scalarField tempBottom = T*magSf*factorsBottom;
        scalarField SfBottom= magSf*factorsBottom;

        // Reduce to scalar
        volumeTop += gSum(volTop);
        areaTop += gSum(SfTop);
        temperatureTop += gSum(tempTop);  

        volumeBottom += gSum(volBottom);
        areaBottom += gSum(SfBottom);
        temperatureBottom += gSum(tempBottom);
    }

    // For the average plenum temperatures, take into account also the faces
    //  on the fuel top/bottom patches
    std::set<label>& topFuelPatchIDs(patchIDs["topFuel"]);
    std::set<label>& bottomFuelPatchIDs(patchIDs["bottomFuel"]);

    for(const label& patchID : topFuelPatchIDs) 
    {
        // Calculate top plenum temperature on the top fuel surface side
        const polyPatch& patch = mesh_.boundary()[patchID].patch();                
        const fvPatchScalarField& T = T_.boundaryField()[patchID];

        if(!isType<emptyPolyPatch>(patch))
        {
            temperatureTop += gSum(T*patch.magFaceAreas());
            areaTop += gSum(patch.magFaceAreas());
        }
        else
        {
            labelList fCells = patch.faceCells();
            scalarField internalT(T.internalField(), fCells);

            temperatureTop += gSum(internalT*patch.magFaceAreas());
            areaTop += gSum(patch.magFaceAreas());
        }
    }


    for(const label& patchID : bottomFuelPatchIDs) 
    {
        // Calculate top plenum temperature on the bottom fuel surface side
        const polyPatch& patch = mesh_.boundary()[patchID].patch();              
        const fvPatchScalarField& T = T_.boundaryField()[patchID];

        if(!isType<emptyPolyPatch>(patch))
        {
            temperatureBottom += gSum(T*patch.magFaceAreas());
            areaBottom += gSum(patch.magFaceAreas());
        }
        else
        {
            labelList fCells = patch.faceCells();
            scalarField internalT(T.internalField(), fCells);

            temperatureBottom += gSum(internalT*patch.magFaceAreas());
            areaBottom += gSum(patch.magFaceAreas());
        }
    }
    
    // Normalize temperature
    temperatureTop /= max(areaTop,VSMALL);
    temperatureBottom /= max(areaBottom,VSMALL);

    // Assigne plena temperature and volumes
    topPlenumT_ = temperatureTop;  
    topPlenumV_ = volumeTop;

    bottomPlenumT_ = temperatureBottom;  
    bottomPlenumV_ = volumeBottom;
}


void Foam::gapFRAPCON::correctDish() 
{       
    const scalarField& Vi = mesh_.V();
    const scalarField& Ti = T_.internalField();

    // Prepare dishFraction scalar field
    scalarField dishFraction(mesh_.nCells(), 0.0);

    const labelListList& matAddrList(mat_.matAddrList());

    forAll(mat_.materialsList(), matI)
    {
        if(isA<fuelMaterial>(mat_.materialsList()[matI]))
        {
            const fuelMaterial& fuelMat = 
            refCast<const fuelMaterial>(mat_.materialsList()[matI]);

            const scalar dishFractionMat
            (fuelMat.dishFraction());

            const labelList& addr(matAddrList[matI]);

            forAll(addr, addrI)
            {
                const label& cellID(addr[addrI]);
                dishFraction[cellID] = dishFractionMat;
            }
        }
    }

    // Total dish volume
    dishV_ = gSum(Vi*dishFraction);

    // Average dish temperature
    dishT_ = gSum(Vi*dishFraction*Ti)/max(dishV_, VSMALL);

    // Total dish V/T
    dishVoverT_ = gSum(Vi*dishFraction/Ti);
}


void Foam::gapFRAPCON::correctCrack() 
{    
    if(mesh_.foundObject<volScalarField>(relocationName_))
    {       
        // Average temperature and total crack volume
        // Necessary to calculate average crack gass temperature
        scalar volume = 0;  
        scalar temperature = 0;
        scalar VoverT = 0;

        // Take a reference to relocation field and cell volumes    
        const scalarField& relocation = 
        mesh_.lookupObject<volScalarField>(relocationName_).internalField();
        const scalarField& cellVol = mesh_.V();  
        
        // Temperature of the cell assumed as temperature of gaas in cracks
        const scalarField& T = T_.internalField();
        
        // The following result comes from supposing relocation as a 2D phenomenon
        const scalarField volRelocation = 2*relocation*cellVol;

        volume = gSum(volRelocation);
        temperature = gSum(T*volRelocation);
        VoverT = gSum(1/T*volRelocation);

        crackV_ = volume;
        crackT_ = temperature/max(volume, VSMALL);
        crackVoverT_ = VoverT;
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gapFRAPCON::gapFRAPCON
(
    const fvMesh& mesh,
    const materials& mat,
    const volScalarField& T,
    const volVectorField& DorDD,
    const fissionGasRelease& fgr,
    const dictionary& gapGasOptDict
)
:
    gapGasModel(mesh, mat, T, DorDD, fgr, gapGasOptDict),
    dishV_(0.0),
    holeV_(0.0),
    topPlenumV_(0.0),
    bottomPlenumV_(0.0),
    crackV_(0.0),
    gapOffsetV_(readScalar(gapGasOptDict.lookup("gapVolumeOffset"))),
    reserveV_(readScalar(gapGasOptDict.lookup("gasReserveVolume"))),
    dishT_(290),
    holeT_(290),
    topPlenumT_(290),
    bottomPlenumT_(290), 
    crackT_(290),
    reserveT_(290),
    reserveTlist_(),
    gapVoverT_(0),
    dishVoverT_(0),
    holeVoverT_(0),
    crackVoverT_(0),
    zMax_(GREAT),
    zMin_(-GREAT),  
    includeCentralHole_
    (gapGasOptDict.lookupOrDefault<bool>("includeCentralHole", true)),
    relocationName_
    (gapGasOptDict.lookupOrDefault<word>("relocationName", "epsilonRelocation"))
{
    // Don't read reserveT if not necessary
    if(reserveV_ > 0)
    {
        if(gapGasOptDict.found("gasReserveTemperatureList"))
        {
    #ifdef OPENFOAMFOUNDATION        
            reserveTlist_.set
            ( 
                new Function1s::Table<scalar>
                (
                    "gasReserveTemperatureList", 
                    gapGasOptDict.subDict("gasReserveTemperatureList")
                )
            );
    #elif OPENFOAMESI        
            reserveTlist_.reset
            ( 
                Function1<scalar>::New
                (
                    "gasReserveTemperatureList", 
                    gapGasOptDict.subDict("gasReserveTemperatureList")
                )
            );
    #endif
            
            scalar currentTime(mesh_.time().value());

            reserveT_ = 
            reserveTlist_->value(mesh_.time().timeToUserTime(currentTime));  
        }
        else
        {
            reserveT_ = readScalar(gapGasOptDict.lookup("gasReserveTemperature"));
        }
    }

    findPatchIDs("gap");
    findPatchIDs("bottomFuel");
    findPatchIDs("topFuel");

    if(includeCentralHole_)
    {   
        findPatchIDs("hole");
    }
    else
    {
        // Insert an empty set
        HashTable<std::set<label>, word>& patchIDs = patchIDsPtr_();
        std::set<label> s;
        patchIDs.insert("hole", s);
    }

    // Set the correction factor list. One list per patch (initial value is 1)
    scalingFactors_.setSize(mesh_.boundaryMesh().size());
    scalingFactorsTopPlenum_.setSize(mesh_.boundaryMesh().size());
    scalingFactorsBottomPlenum_.setSize(mesh_.boundaryMesh().size());
    forAll(scalingFactors_, patchI)
    {
        scalingFactors_[patchI].setSize(mesh_.boundaryMesh()[patchI].size());
        scalingFactorsTopPlenum_[patchI].setSize(mesh_.boundaryMesh()[patchI].size());
        scalingFactorsBottomPlenum_[patchI].setSize(mesh_.boundaryMesh()[patchI].size());

        scalingFactors_[patchI] = 1;
        scalingFactorsTopPlenum_[patchI] = 0;
        scalingFactorsBottomPlenum_[patchI] = 0;
    }

    // Derive species names from gapGas dictionary
    dictionary gapGasDict(readStream(this->typeName));
    close();

    const dictionary& dictMassFr = gapGasDict.subDict("massFractions");
    forAll(dictMassFr.toc(), entryI)
    {
        species_.append(dictMassFr.toc()[entryI]);
    }

    // Read (or build default) dictionary for molar masses
    dictionary dictW("dictW");
    if(gapGasOptDict.found("speciesW"))
    {
        dictW = gapGasOptDict.subDict("speciesW");
    }
    else
    {
        dictW.add("Xe", 131.3);
        dictW.add("Ne", 20.183);
        dictW.add("Ar", 39.948);
        dictW.add("Kr", 83.8);
        dictW.add("Rn", 222);
        dictW.add("He", 4.0026);
    }

    // Set species molar weight
    speciesW_.setSize(species_.size());
    
    forAll(species_, i)
    {
        if(dictW.found(species_[i]))
        {
            speciesW_[i] = readScalar(dictW.lookup(species_[i]));
        }
        else
        {                        
            FatalErrorIn("gapGasModel")
            << "Gas component "
            << species_[i] << nl 
            << " not found in speciesW in gapGasProperties" << exit(FatalError);
        }
    }   

    // Read (or build default) dictionary for conductiviy coefficient A
    dictionary dictAs("coefficientA");
    if(gapGasOptDict.found("conductivity_A"))
    {
        dictAs = gapGasOptDict.subDict("conductivity_A");
    }
    else
    {
        dictAs.add("Xe", 9.825e-5);
        dictAs.add("Ne", 1.0);
        dictAs.add("Ar", 4.092e-4);
        dictAs.add("Kr", 1.966e-4);
        dictAs.add("Rn", 1.0);
        dictAs.add("He", 2.531e-3);
    } 

    // Set list of A coefficients for conductivity formula
    conductivityAs_.setSize(species_.size());
    forAll(species_, i)
    {
        if(dictAs.found(species_[i]))
        {
            conductivityAs_[i] = readScalar(dictAs.lookup(species_[i]));
        }
        else
        {                        
            FatalErrorIn("gapGasModel")
            << "Gas component "
            << species_[i] << nl 
            << " not found in conductivity_A in gapGasProperties" << exit(FatalError);
        }
    }   
    
    // Read (or build default) dictionary for conductiviy coefficient B 
    dictionary dictBs("conductivity_B");
    if(gapGasOptDict.found("conductivity_B"))
    {
        dictBs = gapGasOptDict.subDict("conductivity_B");
    }
    else
    {
        dictBs.add("Xe", 0.7334);
        dictBs.add("Ne", 1.0);
        dictBs.add("Ar", 0.6748);
        dictBs.add("Kr", 0.7006);
        dictBs.add("Rn", 1.0);
        dictBs.add("He", 0.7146);
    } 

    // Set list of B coefficients for conductivity formula
    conductivityBs_.setSize(species_.size()); 
    forAll(species_, i)
    {
        if(dictAs.found(species_[i]))
        {
            conductivityBs_[i] = readScalar(dictBs.lookup(species_[i]));
        }
        else
        {                        
            FatalErrorIn("gapGasModel")
            << "Gas component "
            << species_[i] << nl 
            << " not found in conductivity_B in gapGasProperties" << exit(FatalError);
        }
    }   

    // Set main pin direction
    const globalOptions& globalOpt
    (
        mesh_.lookupObject<globalOptions>("globalOptions")
    );

    pinDirection_ = (globalOpt.pinDirection());
    
    // Read current mass, mass fractions, pressure
    // and calculate molar fractions and gap gas volumes
    readData(readStream(this->typeName));
    close();

    // TODO
    // Temporary fix to change the write path of the file gapGas in time/uniform
    // folders even though it is initialized in the 0/ folder.
    // In the future it will be inizialized in 0/uniform.
    fileName& localPath = const_cast<fileName&>(this->local());
    localPath = "uniform";
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::gapFRAPCON::~gapFRAPCON()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::gapFRAPCON::kappa(const scalar T) const
{        
    scalarList ks_(species_.size(), 0.0);
    
    ks_ = conductivityAs_*pow(T, conductivityBs_);
    
    scalar kMix(0.0);
    
    forAll(speciesW_, specieI)
    {
            scalarList deltas_(species_.size(), 0.0);
            deltas_[specieI] = 1.0;
            scalar summTerm(0.0);
            
            scalar kI = ks_[specieI];
            scalar nI = M_[specieI];
            
            forAll(species_, specieJ)
            {
                    scalar phiIJ(0.0);
                    scalar psiIJ(0.0);
                    scalar kJ = ks_[specieJ];
                    scalar MI = speciesW_[specieI]/1000;
                    scalar MJ = speciesW_[specieJ]/1000;
                    
                    phiIJ = 
                    (
                        pow( (1 + pow(kI/kJ, 0.5)*pow(MI/MJ, 0.25)), 2.0)
                    )/(pow(2, 1.5)*pow(1 + MI/MJ , 0.5));
                    
                    psiIJ = phiIJ*
                    (
                        1 + 2.41*
                        (
                            (MI - MJ)*(MI - 0.142*MJ)
                        )/(pow(MI + MJ, 2.0))
                    );
                    
                    summTerm += (1-deltas_[specieJ])*psiIJ*M_[specieJ];
                    
            }
            
            kMix += ( (kI*nI)/(nI + summTerm) ) ; 
    } 
    
    return kMix;
}


Foam::scalar Foam::gapFRAPCON::a(const scalar T) const
{
    scalar a_He = 0.425 - (2.3e-4*min(T, 1300));
    scalar a_Xe = 0.749 - (2.5e-4*min(T, 1300));
    const scalar W_He = 4.0026;
    const scalar W_Xe = 131.3;
    
    scalar a_mix = 0;
    
    forAll(species_, i)
    {
        scalar W_i = speciesW_[i];
        scalar a_i = a_He + (W_i - W_He)/(W_Xe - W_He)*(a_Xe - a_He);
        
        a_mix += M_[i]*a_i/sqrt(W_i);
    }

    return a_mix;
}


void Foam::gapFRAPCON::correct()
{  
    // Check fgr model for additional gas release and update gas mass, 
    //  mass fractions and molar fractions
    correctMass();

    // Calculate fuel axial coordinates zMin_ and zMax_
    HashTable<std::set<label>, word>& patchIDs = patchIDsPtr_();
    std::set<label>& topFuelPatchIDs(patchIDs["topFuel"]);
    std::set<label>& bottomFuelPatchIDs(patchIDs["bottomFuel"]);

    zMax_ = patchListAvgAxialLocation(topFuelPatchIDs);
    zMin_ = patchListAvgAxialLocation(bottomFuelPatchIDs);

    // Update correction factors
    correctScalingFactors();

    // Update volumes, temperatures and V/T
    correctPlena();
    correctDish();

    correctHole();
    correctGap();
    correctCrack();

    // Update reserve temperature (if time list was given)
    if(reserveTlist_.valid())
    {   
        scalar currentTime(mesh_.time().value());

        reserveT_ = 
        reserveTlist_->value(mesh_.time().timeToUserTime(currentTime));
    }
    
    // Update pressure
    if(gasPType_ == "fromModel")
    {     
        // Total number of moles
        scalar molN(0);
   
        forAll(Y_, nuclideI)
        {
            molN += Y_[nuclideI]*gasM_/(speciesW_[nuclideI]/1000);
        }
    
        // Universal gas constant
        scalar R = Foam::constant::physicoChemical::R.value();

        // Update pressure (ideal gas)
        gasP_ = molN*R/(gapVoverT_ + gapOffsetV_/max(gapT_,VSMALL) +
                        dishVoverT_ +
                        holeVoverT_ +
                        topPlenumV_/max(topPlenumT_,VSMALL) +
                        reserveV_/reserveT_ +
                        bottomPlenumV_/max(bottomPlenumT_,VSMALL) + 
                        crackVoverT_);
    }
    else if(gasPType_ == "fromList")
    {
        scalar currentTime(mesh_.time().value());

        gasP_ = 
        gasPressureList_->value(mesh_.time().timeToUserTime(currentTime));
    }
    else if(gasPType_ == "fixed")
    {
        //gapP does not change
    }
    else
    {    
        FatalErrorIn("gapFRAPCON:correct()")
        << "gasPressureType: "
        << gasPType_ 
        << " not known. Use either:" << nl 
        << "- \"fromModel\"" << nl 
        << "- \"fromList\"" << nl 
        << "- \"fixed\"" << exit(FatalError);
    }

    //*******Output relevant information *****/   
  
    Info<< "Gap conditions:" << endl
        << tab << "Gap volume: " << gapV_ << endl
        << tab << "Gap temperature: " << gapT_ << endl
        << tab << "Reserve volume: " << reserveV_ << endl
        << tab << "Reserve temperature: " << reserveT_ << endl
        << tab << "Dish volume: " << dishV_ << endl
        << tab << "Dish temperature: " << dishT_ << endl  
        << tab << "Hole volume: " << holeV_ << endl
        << tab << "Hole temperature: " << holeT_ << endl
        << tab << "Top plenum volume: " << topPlenumV_ << endl
        << tab << "Top plenum temperature: " << topPlenumT_ << endl
        << tab << "Bottom plenum volume: " << bottomPlenumV_ << endl
        << tab << "Bottom plenum temperature: " << bottomPlenumT_ << endl
        << tab << "Crack gas volume: " << crackV_ << endl
        << tab << "Crack gas average temperature: " << crackT_ << endl
        << tab << "Gas pressure: " << gasP_ << endl    
        << tab << "Gas mass: " << gasM_ << endl    
        << tab << "Y: " << Y_ << endl    
        << tab << "M: " << M_ << endl

        << tab << "pinDirection_: " << pinDirection_ << endl  
        << tab << "zMax: " << zMax_ << endl
        << tab << "zMin: " << zMin_ << endl  

        << endl;        
}


bool Foam::gapFRAPCON::readData(Istream& is)
{    
    dictionary dict(is);

    // The gas mass is either calculated or given in gapGas dict.
    // V, T and p initial conditions are necessary for calculating
    // initial gass mass and number of mols.

    //-Calculate zMin_ and zMax_
    HashTable<std::set<label>, word>& patchIDs = patchIDsPtr_();
    std::set<label>& topFuelPatchIDs(patchIDs["topFuel"]);
    std::set<label>& bottomFuelPatchIDs(patchIDs["bottomFuel"]);

    zMin_ = patchListAvgAxialLocation(bottomFuelPatchIDs);
    zMax_ = patchListAvgAxialLocation(topFuelPatchIDs);

    //-Update scaling factors
    correctScalingFactors();

    // Read pressure
    gasPType_ = dict.lookupOrDefault<word>("gasPressureType", "fromModel");
    if(gasPType_ == "fromModel" or gasPType_ == "fixed")
    {
        gasP_ = readScalar(dict.lookup("gasPressure"));
    }
    else if(gasPType_ == "fromList")
    {
#ifdef OPENFOAMFOUNDATION        
        gasPressureList_.set
        ( 
            new Function1s::Table<scalar>
            (
                "gapPressureList", 
                dict.subDict("gapPressureList")
            )
        );
#elif OPENFOAMESI        
        gasPressureList_.reset
        ( 
            Function1<scalar>::New
            (
                "gapPressureList", 
                dict.subDict("gapPressureList")
            )
        );
#endif        

        scalar currentTime(mesh_.time().value());
        
        gasP_ = 
        gasPressureList_->value(mesh_.time().timeToUserTime(currentTime));  
    }
    else
    {    
        FatalErrorIn("gapFRAPCON:readData")
        << "Gas pressure type "
        << gasPType_ << nl 
        << " not known. Use fromModel, fromList or fixed" << exit(FatalError);
    }

    // Plena volume, temperature and V/T contribution. 
    correctPlena();

    // Hole volume, temperature and V/T contribution. 
    correctHole();

    // Dish volume, temperature and V/T contribution. 
    correctDish();

    // Gap volume, temperature and V/T contribution
    correctGap();

    // Calculate initial crack T, V and V/T contribution
    correctCrack();
    
// Read mass fraction list     
    const dictionary& dictMassFr = dict.subDict("massFractions");

    Y_.resize(species_.size());
    
    forAll(species_, i)
    {
        if(dictMassFr.found(species_[i]))
        {
            Y_[i] = readScalar(dictMassFr.lookup(species_[i]));
        }
        else
        {                        
            FatalErrorIn("gapGasModel")
            << "Gas component "
            << species_[i] << nl 
            << " not found in massFractions in gapGas" << exit(FatalError);
        }
    }   

    correctMassFractions();
    
    Y0_ = Y_; 
    
    // Read mass
    gasM0_ = dict.lookupOrDefault<scalar>("gasMass", 0.0) ;
    
    // if initial mass is not present, calculate it
    if( gasM0_ <= 0.0 )
    {
            calcInitialMass();                                     
    }

    gasM_ = gasM0_;
    
    // Set molar fraction list
    M_.resize(species_.size());

    correctMolFractions();

    M0_ = M_;

    Info<< "Gap conditions:" << endl
        << tab << "Gap volume: " << gapV_ << endl
        << tab << "Gap temperature: " << gapT_ << endl  
        << tab << "Reserve volume: " << reserveV_ << endl
        << tab << "Reserve temperature: " << reserveT_ << endl
        << tab << "Dish volume: " << dishV_ << endl
        << tab << "Dish temperature: " << dishT_ << endl  
        << tab << "Hole volume: " << holeV_ << endl 
        << tab << "Hole temperature: " << holeT_ << endl
        << tab << "Top plenum volume: " << topPlenumV_ << endl
        << tab << "Top plenum temperature: " << topPlenumT_ << endl
        << tab << "Bottom plenum volume: " << bottomPlenumV_ << endl
        << tab << "Bottom plenum temperature: " << bottomPlenumT_ << endl
        << tab << "Crack gas volume: " << crackV_ << endl
        << tab << "Crack gas average temperature: " << crackT_ << endl
        << tab << "Gas pressure: " << gasP_ << endl    
        << tab << "Gas mass: " << gasM_ << endl     
        << tab << "Y: " << Y_ << endl    
        << tab << "M: " << M_ << endl
        << tab << "pinDirection_: " << pinDirection_ << endl  
        << tab << "zMax_: " << zMax_ << endl
        << tab << "zMin_: " << zMin_ << endl  
        << endl;        

    return !is.bad();
}


bool Foam::gapFRAPCON::writeData(Ostream& os) const
{    
    const globalOptions& globalOpt
    (mesh_.lookupObject<globalOptions>("globalOptions"));

    // Obtain model angular fraction to calculate the gap info for the 
    // corresponding 360 degrees 3-D rod
    scalar angularFraction = globalOpt.angularFraction();

    dictionary data;
    
    forAll(species_, i)
    {
        data.add(species_[i], Y_[i]);
    }
    
    os.writeKeyword("massFractions") << data << token::END_STATEMENT 
    << nl << nl; 

    os.writeKeyword("gasPressureType") << gasPType_ << token::END_STATEMENT 
    << nl << nl; 

    os.writeKeyword("gapVolume (in the OFFBEAT model)") 
    << gapV_ << token::END_STATEMENT << nl;
    os.writeKeyword("gapVolume (in the corresponding 360 degree rod)") 
    << gapV_/angularFraction << token::END_STATEMENT << nl;
    os.writeKeyword("gapTemperature") << gapT_ << token::END_STATEMENT 
    << nl << nl;

    os.writeKeyword("reserveVolume (in the OFFBEAT model)") 
    << reserveV_ << token::END_STATEMENT << nl;
    os.writeKeyword("reserveVolume (in the corresponding 360 degree rod)") 
    << reserveV_/angularFraction << token::END_STATEMENT << nl;
    os.writeKeyword("reserveT_") << reserveT_ << token::END_STATEMENT 
    << nl << nl;

    os.writeKeyword("dishVolume (in the OFFBEAT model)") 
    << dishV_ << token::END_STATEMENT << nl;
    os.writeKeyword("dishVolume (in the corresponding 360 degree rod)") 
    << dishV_/angularFraction << token::END_STATEMENT << nl;
    os.writeKeyword("dishTemperature") << dishT_ << token::END_STATEMENT 
    << nl << nl;
    
    os.writeKeyword("holeVolume (in the OFFBEAT model)") 
    << holeV_ << token::END_STATEMENT << nl;
    os.writeKeyword("holeVolume (in the corresponding 360 degree rod)") 
    << holeV_/angularFraction << token::END_STATEMENT << nl;
    os.writeKeyword("holeTemperature") << holeT_ << token::END_STATEMENT 
    << nl << nl;
    
    os.writeKeyword("topPlenumVolume (in the OFFBEAT model)") 
    << topPlenumV_ << token::END_STATEMENT << nl;
    os.writeKeyword("topPlenumVolume (in the corresponding 360 degree rod)") 
    << topPlenumV_/angularFraction << token::END_STATEMENT << nl;
    os.writeKeyword("topPlenumTemperature") << topPlenumT_ 
    << token::END_STATEMENT << nl << nl;

    os.writeKeyword("bottomPlenumVolume (in the OFFBEAT model)") 
    << bottomPlenumV_ << token::END_STATEMENT << nl;
    os.writeKeyword("bottomPlenumVolume (in the corresponding 360 degree rod)") 
    << bottomPlenumV_/angularFraction << token::END_STATEMENT << nl;
    os.writeKeyword("bottomPlenumTemperature") << bottomPlenumT_ 
    << token::END_STATEMENT << nl << nl;

    os.writeKeyword("crackGasVolume (in the OFFBEAT model)") 
    << crackV_ << token::END_STATEMENT << nl;
    os.writeKeyword("crackGasVolume (in the corresponding 360 degree rod)") 
    << crackV_/angularFraction << token::END_STATEMENT << nl;
    os.writeKeyword("crackGasTemperature") << crackT_ << token::END_STATEMENT 
    << nl << nl;   

    os.writeKeyword("gasPressure") << gasP_ << token::END_STATEMENT << nl;
    os.writeKeyword("gasMass") << gasM_ << token::END_STATEMENT << nl;
    
    return os.good();
}

// ************************************************************************* //
