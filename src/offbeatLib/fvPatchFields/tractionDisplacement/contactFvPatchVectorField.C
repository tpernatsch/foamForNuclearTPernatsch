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

\*---------------------------------------------------------------------------*/

#include "contactFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "regionCoupledOFFBEATFvPatch.H"
#include "AMIInterpolation.H"
#include "PrimitivePatchInterpolation.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Private members  * * * * * * * * * * * * *//
    
Foam::tmp<Foam::scalarField> 
Foam::contactFvPatchVectorField::gapWidth() const
{
    const regionCoupledOFFBEATFvPatch& patch
        = refCast<const regionCoupledOFFBEATFvPatch>(regionCoupledPatch_);
    
    const regionCoupledOFFBEATFvPatch& nbrPatch
        = refCast<const regionCoupledOFFBEATFvPatch>(patch.nbrPatch());    

    const fvPatchVectorField& totalDispPatch = 
    ( 
        patch.boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh")
    ) ?
    patch.lookupPatchField<volVectorField, vector>("DD")
    : 
    patch.lookupPatchField<volVectorField, vector>("D"); 

    const fvPatchVectorField& totalDispNbrPatch = 
    ( 
        patch.boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh")
    ) ?
    nbrPatch.lookupPatchField<volVectorField, vector>("DD")
    :
    nbrPatch.lookupPatchField<volVectorField, vector>("D");

    vectorField nf = patch.Sf() / patch.magSf();

    vectorField Cf = patch.Cf() + totalDispPatch;                   
                   
    vectorField nbrCf = regionCoupledPatch_.regionCoupledPatch().interpolate
                      (
                             nbrPatch.Cf() + totalDispNbrPatch
                      );

    return (nbrCf - Cf) & nf;
}


Foam::scalarField contactFvPatchVectorField::boundaryStiffness() const
{
    const regionCoupledOFFBEATFvPatch& patch
        = refCast<const regionCoupledOFFBEATFvPatch>(regionCoupledPatch_);

    const fvPatchField<scalar>& threeK = 
        patch.lookupPatchField<volScalarField, scalar>("threeK");
    
    scalarField K = 1.0/3.0*threeK;

    scalarField volumes(patch.size());
            
    forAll(volumes, cellI)
    {
            volumes[cellI] = patch.boundaryMesh().mesh().V()[patch.faceCells()[cellI]];
    }

    scalarField faceAreas(patch.size());
            
    forAll(faceAreas, cellI)
    {
            faceAreas[cellI] = patch.magSf()[cellI];
    }  

    return (K)*faceAreas/volumes;
}


Foam::scalarField contactFvPatchVectorField::boundaryShearStiffness() const
{
    const regionCoupledOFFBEATFvPatch& patch
        = refCast<const regionCoupledOFFBEATFvPatch>(regionCoupledPatch_);

    const fvPatchField<scalar>& mu = 
        patch.lookupPatchField<volScalarField, scalar>("mu");

    scalarField volumes(patch.size());
            
    forAll(volumes, cellI)
    {
            volumes[cellI] = patch.boundaryMesh().mesh().V()[patch.faceCells()[cellI]];
    }

    scalarField faceAreas(patch.size());
            
    forAll(faceAreas, cellI)
    {
            faceAreas[cellI] = patch.magSf()[cellI];
            // faceAreas[cellI] = sqr(1.0/(2*patch.deltaCoeffs()[cellI]));
    }  

    return (mu)*faceAreas/volumes;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

contactFvPatchVectorField::contactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(p, iF),
    regionCoupledPatch_(refCast<const regionCoupledBaseOFFBEATFvPatch>(p)),
    penaltyScaleFact_(0.1),
    penaltyScaleFactFric_(0.1),    
    penetrationTolerance_(0.0),
    frictionCoeff_(0.0),
    gapWidth_(p.size(), 0),
    interfaceP_(p.size(), 0),
    slip_(p.size(), Foam::vector::zero),
    slip0_(p.size(), Foam::vector::zero),
    frictionTraction_(p.size(), Foam::vector::zero),
    relaxInterP_(1),
    relaxFriction_(1),
    currentTime_(0),
    offset_(0.0),
    rigidMasterNormal_(false),
    rigidMasterFriction_(false)
{}


contactFvPatchVectorField::contactFvPatchVectorField
(
    const contactFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    tractionDisplacementFvPatchVectorField(ptf, p, iF, mapper),
    regionCoupledPatch_(refCast<const regionCoupledBaseOFFBEATFvPatch>(p)),
    penaltyScaleFact_(ptf.penaltyScaleFact_),
    penaltyScaleFactFric_(ptf.penaltyScaleFactFric_),
    penetrationTolerance_(ptf.penetrationTolerance_),
    frictionCoeff_(ptf.frictionCoeff_),
    gapWidth_(ptf.gapWidth_),
    interfaceP_(ptf.interfaceP_),
    slip_(ptf.slip_),
    slip0_(ptf.slip0_),
    frictionTraction_(ptf.frictionTraction_),
    relaxInterP_(ptf.relaxInterP_),
    relaxFriction_(ptf.relaxFriction_),
    currentTime_(ptf.currentTime_),
    offset_(ptf.offset_),
    rigidMasterNormal_(ptf.rigidMasterNormal_),
    rigidMasterFriction_(ptf.rigidMasterFriction_)
{}


contactFvPatchVectorField::contactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    tractionDisplacementFvPatchVectorField(p, iF, dict, false),
    regionCoupledPatch_(refCast<const regionCoupledBaseOFFBEATFvPatch>(p)),
    penaltyScaleFact_(dict.lookupOrDefault<scalar>("penaltyFactor", 0.1)),
    penaltyScaleFactFric_(dict.lookupOrDefault<scalar>("penaltyFactorFriction", 0.1)),
    penetrationTolerance_(dict.lookupOrDefault<scalar>("relativePenetrationTolerance", 0.0)),
    frictionCoeff_(dict.lookupOrDefault<scalar>("frictionCoefficient", 0.0)),
    gapWidth_(p.size(), 0),
    interfaceP_(p.size(), 0),
    slip_(p.size(), Foam::vector::zero),
    slip0_(p.size(), Foam::vector::zero),
    frictionTraction_(p.size(), Foam::vector::zero),
    relaxInterP_(dict.lookupOrDefault<scalar>("relaxInterfacePressure", 1)),
    relaxFriction_(dict.lookupOrDefault<scalar>("relaxFriction", 1)),
    currentTime_(db().time().value()),
    offset_(dict.lookupOrDefault<scalar>("offset", 0)),
    rigidMasterNormal_(dict.lookupOrDefault<bool>("rigidMasterNormal", false)),
    rigidMasterFriction_(dict.lookupOrDefault<bool>("rigidMasterFriction", false))
{    
    if (dict.found("interfaceP"))
    {
        interfaceP_ = scalarField("interfaceP", dict, p.size());
    }
    
    if (dict.found("pressure"))
    {
        pressure_ = scalarField("pressure", dict, patch().size());
    }

    if (dict.found("slip"))
    {
        slip_ = vectorField("slip", dict, p.size());
        slip0_ = slip_;
    }

    if (dict.found("frictionTraction"))
    {
        frictionTraction_ = vectorField("frictionTraction", dict, p.size());
    }
}


contactFvPatchVectorField::contactFvPatchVectorField
(
    const contactFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(ptf, iF),
    regionCoupledPatch_(ptf.regionCoupledPatch_),
    penaltyScaleFact_(ptf.penaltyScaleFact_),
    penaltyScaleFactFric_(ptf.penaltyScaleFactFric_),
    penetrationTolerance_(ptf.penetrationTolerance_),
    frictionCoeff_(ptf.frictionCoeff_),
    gapWidth_(ptf.gapWidth_),
    interfaceP_(ptf.interfaceP_),
    slip_(ptf.slip_),
    slip0_(ptf.slip0_),
    frictionTraction_(ptf.frictionTraction_),
    relaxInterP_(ptf.relaxInterP_),
    relaxFriction_(ptf.relaxFriction_),
    currentTime_(ptf.currentTime_),
    offset_(ptf.offset_),
    rigidMasterNormal_(ptf.rigidMasterNormal_),
    rigidMasterFriction_(ptf.rigidMasterFriction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void contactFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    tractionDisplacementFvPatchVectorField::autoMap(m);
#ifdef OPENFOAMFOUNDATION    
    m(gapWidth_, gapWidth_);
    m(interfaceP_, interfaceP_);
#elif OPENFOAMESI
    gapWidth_.autoMap(m);
    interfaceP_.autoMap(m);
#endif    
}


void contactFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    tractionDisplacementFvPatchVectorField::rmap(ptf, addr);

    const contactFvPatchVectorField& dmptf =
        refCast<const contactFvPatchVectorField>(ptf);

    gapWidth_.rmap(dmptf.gapWidth_, addr);
    interfaceP_.rmap(dmptf.interfaceP_, addr);
}

void contactFvPatchVectorField::updateTraction()
{
    const regionCoupledOFFBEATFvPatch& patch
        = refCast<const regionCoupledOFFBEATFvPatch>(regionCoupledPatch_);
    
    const regionCoupledOFFBEATFvPatch& nbrPatch
        = refCast<const regionCoupledOFFBEATFvPatch>(patch.nbrPatch());

    const contactFvPatchVectorField& nbr = this->neighbour();

    if(currentTime_ < db().time().value())
    {
        slip0_ = slip_;
        currentTime_ = db().time().value();
    }

    if (!regionCoupledPatch_.owner() )
    {
        penaltyScaleFact_ = nbr.penaltyScaleFact_;
        penaltyScaleFactFric_ = nbr.penaltyScaleFactFric_;
        frictionCoeff_ = nbr.frictionCoeff_;

        const AMIInterpolation& ami(nbrPatch.regionCoupledPatch().AMI());

        //*** Calculate gap width ***//

        // Update the gap width
        gapWidth_ = 
        // 0.5*
        (
            gapWidth() + offset_
          // + patch.regionCoupledPatch().interpolate(nbr.gapWidth())
        );

        // forAll(gapWidth_, faceI)
        // {
        //     if(!ami.tgtAddress()[faceI].size())
        //     {
        //         gapWidth_[faceI] = GREAT;
        //     }
        // }   
        
        // Update the neighbour gap width
        nbr.gapWidth_ = nbrPatch.regionCoupledPatch().interpolate(gapWidth_);     

        // Debugging
        const_cast<fvPatchScalarField&>
        (
            patch.lookupPatchField<volScalarField, scalar>("gapWidth")
        ) = gapWidth_;   

        // Debugging
        const_cast<fvPatchScalarField&>
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>("gapWidth")
        ) = nbr.gapWidth_;         

        // Interpolate back to this patch, so to avoid scalarFields of different lenght
        const scalarField nbrBoundaryStiff(
            patch.regionCoupledPatch().interpolate(nbr.boundaryStiffness()));

        scalarField interfacePOld = interfaceP_;

        scalarField penaltyFact = penaltyScaleFact_*min(
            boundaryStiffness(), nbrBoundaryStiff);
            
        // Calculated new values for penalty factors
        scalarField avgDelta = 0.5*(1/patch.deltaCoeffs() + 1/max(
            patch.regionCoupledPatch().interpolate(nbrPatch.deltaCoeffs()), SMALL));

        //*** Calculate normal interface pressure ***//
        forAll(gapWidth_, faceI)
        {
            if(gapWidth_[faceI] > 0.0)
            {
                // This corresponds to open gap or gap re-opening 
                interfaceP_[faceI] = 0.0; 
            }
            else if(-gapWidth_[faceI]/avgDelta[faceI] > penetrationTolerance_)
            {
                // Calculate the interface pressure
                interfaceP_[faceI] = max
                (
                    -penaltyFact[faceI]*gapWidth_[faceI] ,
                    0.0
                );

            }
        }
        
        interfaceP_ = (relaxInterP_*interfaceP_ + 
                    (1 - relaxInterP_)*interfacePOld);

        const_cast<fvPatchScalarField&>
        (
            patch.lookupPatchField<volScalarField, scalar>("interfaceP")
        ) = interfaceP_; 

        if(!rigidMasterNormal_)
        {
            // Interpolate contact pressure on the master side
            nbr.interfaceP_ = nbrPatch.regionCoupledPatch().interpolate(interfaceP_);

            const_cast<fvPatchScalarField&>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>("interfaceP")
            ) = nbr.interfaceP_;
        }

        //*** Calculate friction stress ***//

        // Prepare slave and master DD fields for slip/friction calculation
        vectorField ownDD(patch.size(), Foam::vector::zero);
        vectorField nbrDD(patch.size(), Foam::vector::zero);

        if(db().foundObject<fvMesh>("referenceMesh"))
        {
            const volVectorField& DD =
                db().lookupObject<volVectorField>("DD");

            ownDD = DD.boundaryField()[patch.index()];
            nbrDD = patch.regionCoupledPatch().interpolate(
                DD.boundaryField()[nbrPatch.index()]);
        }
        else
        {
            const volVectorField& D =
                db().lookupObject<volVectorField>("D");

            ownDD = D.boundaryField()[patch.index()] - 
            D.oldTime().boundaryField()[patch.index()];

            nbrDD = patch.regionCoupledPatch().interpolate(
                D.boundaryField()[nbrPatch.index()] -
                D.oldTime().boundaryField()[nbrPatch.index()]);            
        }

        // Increase slip (only when penetration is > 0)
        // When gapWidth is zero total slip is equal to zero
        slip_ = slip0_ + neg(gapWidth_)*(
            (I - sqr(patch.nf())) & (ownDD - nbrDD));
        slip_ = neg(gapWidth_)* slip_;

        const scalarField nbrBoundaryShearStiff(
            patch.regionCoupledPatch().interpolate(nbr.boundaryShearStiffness()));

        scalarField penaltyFactFric = penaltyScaleFactFric_*min(
            boundaryShearStiffness(), nbrBoundaryShearStiff);

        // Calculate friction traction magnitude
        scalarField frictionMag(min(
            penaltyFactFric*mag(slip_), frictionCoeff_*interfaceP_));

        // Calculate friction traction vector and relax
        frictionTraction_ = relaxFriction_*frictionMag*(-slip_/max(mag(slip_), VSMALL)) + 
        (1 - relaxFriction_)*frictionTraction_;  

        if(!rigidMasterFriction_)
        {
            // Interpolate slip and friction on the master side (change sign)
            nbr.slip_ = -nbrPatch.regionCoupledPatch().interpolate(slip_);
            nbr.frictionTraction_ = -nbrPatch.regionCoupledPatch().interpolate(
                frictionTraction_);
        }

        // Correct for area mismatch
        interfaceP_ *= ami.tgtWeightsSum();
        frictionTraction_ *= ami.tgtWeightsSum();
        
        // Adjust to conserve total friction force
        vector Ft = gSum(frictionTraction_*patch.magSf());
        vector nbrFt = gSum(nbr.frictionTraction_*nbrPatch.magSf());
        forAll(frictionTraction_, i)
        {
            frictionTraction_[i].x() *= nbrFt.x()/max(Ft.x(), SMALL);
            frictionTraction_[i].y() *= nbrFt.y()/max(Ft.y(), SMALL);
            frictionTraction_[i].z() *= nbrFt.z()/max(Ft.z(), SMALL);
        }
        
        // Adjust to conserve total pressure force
        vector Fp = gSum(interfaceP_*patch.Sf());
        vector nbrFp = gSum(nbr.interfaceP_*nbrPatch.Sf());
        interfaceP_ *= mag(nbrFp)/max(mag(Fp), SMALL);
        
        Fp = gSum(interfaceP_*patch.Sf());
        scalar normFactor = max(gSum(interfaceP_*patch.magSf()), SMALL);
        frictionTraction_ += (nbrFp + Fp) / normFactor * interfaceP_;
        
        if (debug)
        {
            Info<< "Owner pressure (N): " << gSum(interfaceP_*patch.Sf()) << nl
                << "Neighbour pressure (N): " << gSum(nbr.interfaceP_*nbrPatch.Sf()) << nl
                << "Owner traction (N): " << gSum(frictionTraction_*patch.magSf()) << nl
                << "Neighbour traction (N): " << gSum(nbr.frictionTraction_*nbrPatch.magSf()) << nl
                << "Owner total (N): " << gSum(interfaceP_*patch.Sf() 
                                               - frictionTraction_*patch.magSf()) << nl
                << "Neighbour total (N): " << gSum(nbr.interfaceP_*nbrPatch.Sf() 
                                                   - nbr.frictionTraction_*nbrPatch.magSf()) << nl
                << endl;
        }
    }
    
    pressure_ = interfaceP_;
    traction_ = frictionTraction_;
    
    // TODO: 
    // - should we add an augemented lagrangian version of the search algorithm?
}


void contactFvPatchVectorField::write(Ostream& os) const
{
    tractionDisplacementFvPatchVectorField::write(os);

#ifdef OPENFOAMFOUNDATION        
    if (regionCoupledPatch_.owner()) 
    {    
        writeEntryIfDifferent<scalar>(os, "penaltyFactor", 0.1, penaltyScaleFact_);
        writeEntryIfDifferent<scalar>(os, "penaltyFactorFriction", 0.1, penaltyScaleFactFric_);
        writeEntryIfDifferent<scalar>(os, "frictionCoefficient", 0.0, frictionCoeff_);
        writeEntryIfDifferent<scalar>(os, "offset", 0.0, offset_);
        writeEntryIfDifferent<bool>(os, "rigidMasterNormal", false, rigidMasterNormal_);
        writeEntryIfDifferent<bool>(os, "rigidMasterFriction", false, rigidMasterFriction_);

    }

    writeEntryIfDifferent<scalar>(os, "relaxInterfacePressure", 1, relaxInterP_);
    writeEntryIfDifferent<scalar>(os, "relaxFriction", 1, relaxFriction_);
    
    writeEntry(os, "interfaceP", interfaceP_);
    writeEntry(os, "slip", slip_);
    writeEntry(os, "frictionTraction", frictionTraction_);
#elif OPENFOAMESI
    
    if (regionCoupledPatch_.owner()) 
    {    
        os.writeEntryIfDifferent<scalar>("penaltyFactor", 0.1, penaltyScaleFact_);
        os.writeEntryIfDifferent<scalar>("penaltyFactorFriction", 0.1, penaltyScaleFactFric_);
        os.writeEntryIfDifferent<scalar>("frictionCoefficient", 0.0, frictionCoeff_);
        os.writeEntryIfDifferent<scalar>("offset", 0.0, offset_);
        os.writeEntryIfDifferent<bool>("rigidMasterNormal", false, rigidMasterNormal_);
        os.writeEntryIfDifferent<bool>("rigidMasterFriction", false, rigidMasterFriction_);

    }

    os.writeEntryIfDifferent<scalar>("relaxInterfacePressure", 1, relaxInterP_);
    os.writeEntryIfDifferent<scalar>("relaxFriction", 1, relaxFriction_);
    
    interfaceP_.writeEntry("interfaceP", os);
    slip_.writeEntry("slip", os);
    frictionTraction_.writeEntry("frictionTraction", os);   
#endif      
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    contactFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
