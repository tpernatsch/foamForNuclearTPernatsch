 /*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     |
     \\  /    A nd           | www.openfoam.com
      \\/     M anipulation  |
 -------------------------------------------------------------------------------
     Copyright (C) 2020-2021 OpenCFD Ltd.
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
 
 #include "mixedMappedOFFvPatchScalarField.H"
 #include "volFields.H"
 #include "interpolationCell.H"
 #include "mappedFixedValueFvPatchField.H"
 #include "mappedPatchFieldBase.H"
 #include "addToRunTimeSelectionTable.H"

 
 // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
 

 Foam::mixedMappedOFFvPatchScalarField::mixedMappedOFFvPatchScalarField
 (
     const fvPatch& p,
     const DimensionedField<scalar, volMesh>& iF
 )
 :
     mixedFvPatchField<scalar>(p, iF),
     mappedPatchFieldBase<scalar>
     (
         mappedFixedValueFvPatchField<scalar>::mapper(p, iF),
         *this
     ),
     thicknessLayers_(),
     kappaLayers_(),
     heatFlux_()
 {
     this->refValue() = Zero;
     this->refGrad() = Zero;
     this->valueFraction() = 0.0;
 }
 
 

 Foam::mixedMappedOFFvPatchScalarField::mixedMappedOFFvPatchScalarField
 (
     const fvPatch& p,
     const DimensionedField<scalar, volMesh>& iF,
     const dictionary& dict
 )
 :
     // Bypass dictionary constructor (all reading handled later)
     // but cannot use NO_READ since will still trigger an evaluate()
     mixedFvPatchField<scalar>(p, iF),
     mappedPatchFieldBase<scalar>
     (
         mappedFixedValueFvPatchField<scalar>::mapper(p, iF),
         *this,
         dict
     ),
     thicknessLayers_(dict.getOrDefault<scalarList>("thicknessLayers", scalarList::null())),
     kappaLayers_(dict.getOrDefault<scalarList>("kappaLayers", scalarList::null())),
     heatFlux_(scalarList::null())
 {
     fvPatchFieldBase::readDict(dict);  // Consistent with a dict constructor
 
     this->readValueEntry(dict, IOobjectOption::MUST_READ);
 
     if (this->readMixedEntries(dict))
     {
         // Full restart
     }
     else
     {
         // Start from user entered data. Assume fixedValue.
         this->refValue() = *this;
         this->refGrad() = Zero;
         this->valueFraction() = 1.0;
     }
 
 // This blocks (crashes) with more than two worlds!
 //
 }
 
 
 Foam::mixedMappedOFFvPatchScalarField::mixedMappedOFFvPatchScalarField
 (
     const mixedMappedOFFvPatchScalarField& ptf,
     const fvPatch& p,
     const DimensionedField<scalar, volMesh>& iF,
     const fvPatchFieldMapper& mapper
 )
 :
     mixedFvPatchField<scalar>(ptf, p, iF, mapper),
     mappedPatchFieldBase<scalar>
     (
         mappedFixedValueFvPatchField<scalar>::mapper(p, iF),
         *this,
         ptf
     ),
     thicknessLayers_(ptf.thicknessLayers_),
     kappaLayers_(ptf.kappaLayers_),
     heatFlux_(ptf.heatFlux_)
 
{}
 Foam::mixedMappedOFFvPatchScalarField::mixedMappedOFFvPatchScalarField
 (
     const mixedMappedOFFvPatchScalarField& ptf
 )
 :
     mixedFvPatchField<scalar>(ptf),
     mappedPatchFieldBase<scalar>(ptf),
     thicknessLayers_(ptf.thicknessLayers_),
     kappaLayers_(ptf.kappaLayers_),
     heatFlux_(ptf.heatFlux_)
 {}
 

 Foam::mixedMappedOFFvPatchScalarField::mixedMappedOFFvPatchScalarField
 (
     const mixedMappedOFFvPatchScalarField& ptf,
     const DimensionedField<scalar, volMesh>& iF
 )
 :
     mixedFvPatchField<scalar>(ptf, iF),
     mappedPatchFieldBase<scalar>
     (
         mappedFixedValueFvPatchField<scalar>::mapper(ptf.patch(), iF),
         *this,
         ptf
     ),
     thicknessLayers_(ptf.thicknessLayers_),
     kappaLayers_(ptf.kappaLayers_),
     heatFlux_(ptf.heatFlux_)
 {}
 
 
 // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
 void Foam::mixedMappedOFFvPatchScalarField::autoMap(const fvPatchFieldMapper& m)
 {
     mixedFvPatchField<scalar>::autoMap(m);
 }
 

 void Foam::mixedMappedOFFvPatchScalarField::rmap
 (
     const fvPatchField<scalar>& ptf,
     const labelList& addr
 )
 {
     mixedFvPatchField<scalar>::rmap(ptf, addr);
 }

 void Foam::mixedMappedOFFvPatchScalarField::getWeights
(
    scalarField& myWeights,
    scalarField& nbrWeights
) const
{
    
    scalarField tk(patch().lookupPatchField<volScalarField, Foam::scalar>("k"));

        // Optionally modify with explicit resistance
    if (thicknessLayers_.size())
    {
        scalarField KDelta(tk*patch().deltaCoeffs());
 
        // Harmonic averaging of kappa*deltaCoeffs
        {
            KDelta = 1.0/KDelta;
            if (thicknessLayers_.size())
            {
                forAll(thicknessLayers_, iLayer)
                {
                    KDelta += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
                }
            }
            KDelta = 1.0/KDelta;
        }
 
        // Update kappa from KDelta
        tk = KDelta/patch().deltaCoeffs();
    }
    
    
   // myWeights = patchField_.patch().deltaCoeffs()*tk;
   myWeights=tk*patchField_.patch().deltaCoeffs();

    if (mapper_.sameWorld())
    {
        // Same world so lookup
        const auto& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());
        const label nbrPatchID = mapper_.samplePolyPatch().index();
        const auto& nbrPatch = nbrMesh.boundary()[nbrPatchID];

        if (nbrMesh.foundObject<volScalarField>("htc"))
        {
            nbrWeights=nbrPatch.lookupPatchField<volScalarField, scalar>("htc").patchInternalField();
        }
        else if (nbrMesh.foundObject<volScalarField>("htc.liquid.structure"))
        {
            nbrWeights=nbrPatch.lookupPatchField<volScalarField, scalar>("htc.liquid.structure").patchInternalField()+nbrPatch.lookupPatchField<volScalarField, scalar>("htc.vapour.structure").patchInternalField();
        }
        else
        {
            WarningInFunction
            << "Weighting convective fields not found! "<<endl
            <<"Setting the weights to the delta coeffs. Have you named your phases differently from liquid and vapour?"
            << endl;
            nbrWeights=nbrPatch.deltaCoeffs();
        }

        
    }

    else
    {
    // Swap to obtain full local values of neighbour internal field
    // Different world so use my region,patch. Distribution below will
    // do the reordering
    nbrWeights = myWeights;
    }

 
    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;
 
    distribute(fieldName_ + "_weights", nbrWeights);
 
    // Restore tag
    UPstream::msgType() = oldTag;


}


 void Foam::mixedMappedOFFvPatchScalarField::updateCoeffs()
 {
     if (this->updated())
     {
         return;
     }
 
     const tmp<Field<scalar>> nbrIntFld(mapAveragedTemperature()); // give Tsol, get either T or TAvg (see GeNFoam side)

     

     //- Weighted
     scalarField myWeights;
     scalarField nbrWeights;

     getWeights(myWeights, nbrWeights); // give k/d, get back either h or hl+hg

     heatFlux_= myWeights* (*this - patchInternalField());

     //Info <<"Heat flux from offbeat to GeNFoam: "<< heatFlux_<<endl;


     Info <<"Power transmitted from offbeat to GeNFoam: "<< -1*gSum(heatFlux_*patch().magSf())<<" W" <<endl;

    
  
 
     this->refValue() = nbrIntFld;
     this->refGrad() = Zero;
     this->valueFraction() = nbrWeights/(nbrWeights + myWeights);
 
     mixedFvPatchField<scalar>::updateCoeffs();
 
     if (debug)
     {
         Info<< this->patch().boundaryMesh().mesh().name() << ':'
             << this->patch().name() << ':'
             << this->internalField().name() << " <- "
             << this->mapper_.sampleRegion() << ':'
             << this->mapper_.samplePatch() << ':'
             << this->fieldName_ << " :"
             << " value "
             << " min:" << gMin(*this)
             << " max:" << gMax(*this)
             << " avg:" << gAverage(*this)
             << endl;
     }
 }
 
 
 void Foam::mixedMappedOFFvPatchScalarField::write(Ostream& os) const
 {
     mappedPatchFieldBase<scalar>::write(os);
     mixedFvPatchField<scalar>::write(os);
 }


  Foam::tmp<Foam::Field<Foam::scalar>>  Foam::mixedMappedOFFvPatchScalarField::mapAveragedTemperature() const
 {
        // Swap to obtain full local values of neighbour internal field
    tmp<Field<scalar>> tnbrIntFld(new Field<scalar>());
    scalarField& nbrIntFld = tnbrIntFld.ref();
 
    if (mapper_.sameWorld())
    {
        // Same world so lookup
        const auto& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());
        const label nbrPatchID = mapper_.samplePolyPatch().index();
        const auto& nbrPatch = nbrMesh.boundary()[nbrPatchID];

        if (nbrMesh.foundObject<volScalarField>("htc"))
            nbrIntFld = nbrPatch.lookupPatchField<volScalarField, scalar>("T").patchInternalField(); // get the liquid side temperature
        else if (nbrMesh.foundObject<volScalarField>("htc.liquid.structure"))
        {
            scalarField htc1=nbrPatch.lookupPatchField<volScalarField, scalar>("htc.liquid.structure").patchInternalField(); 
            scalarField htc2 = nbrPatch.lookupPatchField<volScalarField, scalar>("htc.vapour.structure").patchInternalField();
            scalarField T1 = nbrPatch.lookupPatchField<volScalarField, scalar>("T.liquid").patchInternalField();
            scalarField T2 = nbrPatch.lookupPatchField<volScalarField, scalar>("T.vapour").patchInternalField();
            nbrIntFld = T1*htc1/(htc1+htc2)+T2*htc2/(htc1+htc2);
        }
    }
    else
    {
        nbrIntFld=patchField_.patchInternalField();
    }
 
    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;
    distribute(fieldName_, nbrIntFld);
 
    // Restore tag
    UPstream::msgType() = oldTag;
 
    return tnbrIntFld;
 }
 
 
 // ************************************************************************* //

  // * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(mixedMappedOFFvPatchScalarField, 0);
addToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    mixedMappedOFFvPatchScalarField
);
}

// ************************************************************************* //