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
        
   /* vectorField n(patch.size(), Foam::vector::zero);          
    
    if (!regionCoupledPatch_.owner() )
    {
        const AMIInterpolation& ami(nbrPatch.regionCoupledPatch().AMI()); 
                   
        vectorField nbrNf = -patch.regionCoupledPatch().interpolate(ami.srcNf());

        n += ami.tgtNf();
        n+=nbrNf;
        n /= mag(n);
    }
    else
    {
        const AMIInterpolation& ami(patch.regionCoupledPatch().AMI()); 
                   
        vectorField nbrNf = -patch.regionCoupledPatch().interpolate(ami.tgtNf());

        n += ami.srcNf();
        n+=nbrNf;
        n /= mag(n);        
    }   

    vectorField Cf = patch.Cf()
                   + totalDispPatch;                   
                   
    vectorField nbrCf = regionCoupledPatch_.regionCoupledPatch().interpolate
                      (
                             nbrPatch.Cf()
                             + totalDispNbrPatch
                      );*/

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
    gapWidth_(p.size(), 0),
    interfaceP_(p.size(), 0),
    outerPressure_(p.size(), 0),
    forceConcentricSlave(false),
    relaxInterP_(1)
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
    gapWidth_(ptf.gapWidth_),
    interfaceP_(ptf.interfaceP_),
    outerPressure_(ptf.outerPressure_),
    forceConcentricSlave(ptf.forceConcentricSlave),
    relaxInterP_(ptf.relaxInterP_)
{}


contactFvPatchVectorField::contactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    tractionDisplacementFvPatchVectorField(p, iF,dict),
    regionCoupledPatch_(refCast<const regionCoupledBaseOFFBEATFvPatch>(p)),
    penaltyScaleFact_(dict.lookupOrDefault<scalar>("penaltyFactor", 0.1)),
    gapWidth_(p.size(), 0),
    interfaceP_(p.size(), 0),
    outerPressure_(p.size(), 0),
    forceConcentricSlave(dict.lookupOrDefault<bool>("forceConcentricSlave", false)),
    relaxInterP_(dict.lookupOrDefault<scalar>("relaxInterfacePressure", 1))
{
    if (dict.found("outerPressure"))
    {
        outerPressure_ = scalarField("outerPressure", dict, p.size());
    }
    else
    {
        outerPressure_ = scalarField("pressure", dict, p.size());
    }
    
    if (dict.found("interfaceP"))
    {
        interfaceP_ = scalarField("interfaceP", dict, p.size());
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
    gapWidth_(ptf.gapWidth_),
    interfaceP_(ptf.interfaceP_),
    outerPressure_(ptf.outerPressure_),
    forceConcentricSlave(ptf.forceConcentricSlave),
    relaxInterP_(ptf.relaxInterP_)
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

void contactFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const regionCoupledOFFBEATFvPatch& patch
        = refCast<const regionCoupledOFFBEATFvPatch>(regionCoupledPatch_);
    
    const regionCoupledOFFBEATFvPatch& nbrPatch
        = refCast<const regionCoupledOFFBEATFvPatch>(patch.nbrPatch());

    const contactFvPatchVectorField& nbr = this->neighbour();

    scalarField correctionFactor(patch.size(), 1);    
    
    // if(currentTime_ < db().time().value())
    // {
    //     slip0_ = slip_;
    //     currentTime_ = db().time().value();
    // }

    if (!regionCoupledPatch_.owner() )
    {        
        penaltyScaleFact_ = nbr.penaltyScaleFact_;

        // const AMIInterpolation& ami(nbrPatch.regionCoupledPatch().AMI());
        
        // Update the gap width
        gapWidth_ = 
        // 0.5*
        (
            gapWidth() //+ offset_
          // + patch.regionCoupledPatch().interpolate(nbr.gapWidth())
        );
        
        // Update the neighbour gap width
        nbr.gapWidth_ = nbrPatch.regionCoupledPatch().interpolate(gapWidth_);  

        // forAll(gapWidth_, faceI)
        // {
        //     if(!ami.tgtAddress()[faceI].size())
        //     {
        //         gapWidth_[faceI] = GREAT;
        //     }
        // }     

        // forAll(nbr.gapWidth_, faceI)
        // {
        //     if(!ami.srcAddress()[faceI].size())
        //     {
        //         nbr.gapWidth_[faceI] = GREAT;
        //     }
        // }   

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

        // Calculate the interface pressure
        interfaceP_ = max
        (
           -penaltyFact*gapWidth_ ,
            0.0
        );

        interfaceP_ = (relaxInterP_*interfaceP_ + 
            (1 - relaxInterP_)*interfacePOld);

        const_cast<fvPatchScalarField&>
        (
            patch.lookupPatchField<volScalarField, scalar>("interfaceP")
        ) = interfaceP_;  

        // if(!rigidMasterNormal_)
        // {
            // Interpolate contact pressure on the master side
            nbr.interfaceP_ = nbrPatch.regionCoupledPatch().interpolate(interfaceP_);

            const_cast<fvPatchScalarField&>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>("interfaceP")
            ) = nbr.interfaceP_;
        // }

        // //*** Calculate friction stress ***//

        // // Prepare slave and master DD fields for slip/friction calculation
        // vectorField ownDD(patch.size(), Foam::vector::zero);
        // vectorField nbrDD(patch.size(), Foam::vector::zero);

        // if(db().foundObject<fvMesh>("referenceMesh"))
        // {
        //     const volVectorField& DD =
        //         db().lookupObject<volVectorField>("DD");

        //     ownDD = DD.boundaryField()[patch.index()];
        //     nbrDD = patch.regionCoupledPatch().interpolate(
        //         DD.boundaryField()[nbrPatch.index()]);
        // }
        // else
        // {
        //     const volVectorField& D =
        //         db().lookupObject<volVectorField>("D");

        //     ownDD = D.boundaryField()[patch.index()] - 
        //     D.oldTime().boundaryField()[patch.index()];

        //     nbrDD = patch.regionCoupledPatch().interpolate(
        //         D.boundaryField()[nbrPatch.index()] -
        //         D.oldTime().boundaryField()[nbrPatch.index()]);            
        // }

        // // Increase slip (only when penetration is > 0)
        // // When gapWidth is zero total slip is equal to zero
        // slip_ = slip0_ + neg(gapWidth_)*(
        //     (I - sqr(patch.nf())) & (ownDD - nbrDD));
        // slip_ = neg(gapWidth_)* slip_;

        // const scalarField nbrBoundaryShearStiff(
        //     patch.regionCoupledPatch().interpolate(nbr.boundaryShearStiffness()));

        // scalarField penaltyFactFric = penaltyScaleFactFric_*min(
        //     boundaryShearStiffness(), nbrBoundaryShearStiff);

        // // Calculate friction traction magnitude
        // scalarField frictionMag(min(
        //     penaltyFactFric*mag(slip_), frictionCoeff_*interfaceP_));

        // // Calculate friction traction vector and relax
        // frictionTraction_ = relaxFriction_*frictionMag*(-slip_/max(mag(slip_), VSMALL)) + 
        // (1 - relaxFriction_)*frictionTraction_;  

        // if(!rigidMasterFriction_)
        // {
        //     // Interpolate slip and friction on the master side (change sign)
        //     nbr.slip_ = -nbrPatch.regionCoupledPatch().interpolate(slip_);
        //     nbr.frictionTraction_ = -nbrPatch.regionCoupledPatch().interpolate(
        //         frictionTraction_);
        // }

        // // Update correction factors
        // correctionFactor = ami.tgtWeightsSum();

    }

    vectorField n(patch.nf());
    
    // vectorField n(patch.size(), Foam::vector::zero);    
    // if (!regionCoupledPatch_.owner() )
    // {
    //     const AMIInterpolation& ami(nbrPatch.regionCoupledPatch().AMI()); 
                   
    //     vectorField nbrNf = -patch.regionCoupledPatch().interpolate(ami.srcNf());

    //     n+=ami.tgtNf();
    //     n+=nbrNf;
    //     n /= mag(n);
    // }
    // else
    // {
    //     const AMIInterpolation& ami(patch.regionCoupledPatch().AMI()); 
                   
    //     vectorField nbrNf = -patch.regionCoupledPatch().interpolate(ami.tgtNf());

    //     n+=ami.srcNf();
    //     n+=nbrNf;
    //     n /= mag(n);        
    // }

    pressure_ = max(outerPressure_, interfaceP_)*correctionFactor;
    traction_ = vector::zero;
    // traction_ = frictionTraction_;


    tractionDisplacementFvPatchVectorField::updateCoeffs();
    
}


void contactFvPatchVectorField::write(Ostream& os) const
{
    tractionDisplacementFvPatchVectorField::write(os);
    
#ifdef OPENFOAMFOUNDATION    
    if (regionCoupledPatch_.owner())
    {
        writeEntryIfDifferent<scalar>(os, "penaltyFactor", 0.1, penaltyScaleFact_);
    }

    writeEntryIfDifferent<bool>(os, "forceConcentricSlave", false, forceConcentricSlave);
    writeEntryIfDifferent<scalar>(os, "relaxInterfacePressure", 1, relaxInterP_);
    
    writeEntry(os, "interfaceP", interfaceP_);
    writeEntry(os, "outerPressure", outerPressure_);
#elif OPENFOAMESI
    if (regionCoupledPatch_.owner())
    {
        os.writeEntryIfDifferent<scalar>("penaltyFactor", 0.1, penaltyScaleFact_);
    }

    os.writeEntryIfDifferent<bool>("forceConcentricSlave", false, forceConcentricSlave);
    os.writeEntryIfDifferent<scalar>("relaxInterfacePressure", 1, relaxInterP_);
    
    os.writeEntry("interfaceP", interfaceP_);
    os.writeEntry( "outerPressure", outerPressure_);    
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
