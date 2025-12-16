/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2506                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.          |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "scaledDiscontinuousCyclicAMIFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "primitiveMeshTools.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scaledDiscontinuousCyclicAMIFvPatchVectorField::scaledDiscontinuousCyclicAMIFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    discontinuousCyclicAMIFvPatchField<vector>(p, iF),
    scalingFactors_(p.size(),1),
    sigma_(p.size(),0),
    fluxCorrectionFactors_(p.size(),1),
    connectedPatches_(0),
    correctFluxAndAreas_(false),
    underRelaxation_(1)
{}


Foam::scaledDiscontinuousCyclicAMIFvPatchVectorField::scaledDiscontinuousCyclicAMIFvPatchVectorField
(
    const scaledDiscontinuousCyclicAMIFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    discontinuousCyclicAMIFvPatchField<vector>(ptf, p, iF, mapper),
    scalingFactors_(ptf.scalingFactors_),
    sigma_(ptf.sigma_),
    fluxCorrectionFactors_(ptf.fluxCorrectionFactors_),
    connectedPatches_(ptf.connectedPatches_),
    correctFluxAndAreas_(ptf.correctFluxAndAreas_),
    underRelaxation_(ptf.underRelaxation_)
{}


Foam::scaledDiscontinuousCyclicAMIFvPatchVectorField::scaledDiscontinuousCyclicAMIFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    discontinuousCyclicAMIFvPatchField<vector>(p, iF, dict),
    scalingFactors_(p.size(),1),
    sigma_(p.size(),0),
    fluxCorrectionFactors_(p.size(),1),
    connectedPatches_(dict.getOrDefault<wordList>("connectedPatches", wordList(0))),
    correctFluxAndAreas_(dict.getOrDefault<bool>("correctFluxAndAreas", false)),
    underRelaxation_(dict.getOrDefault<scalar>("underRelaxation", 1))
{

    scalingFactors_.assign("scalingFactors", dict, p.size(), IOobjectOption::LAZY_READ);
    sigma_.assign("sigma", dict, p.size(), IOobjectOption::LAZY_READ);
    fluxCorrectionFactors_.assign("fluxCorrectionFactors", dict, p.size(), IOobjectOption::LAZY_READ);


    if(correctFluxAndAreas_)
    {
        scalar nBranch = connectedPatches_.size()+1;

        // Rescale polypatch areas

        discontinuousCyclicAMIPolyPatch& discPP = 
            const_cast<discontinuousCyclicAMIPolyPatch&>
            (
                this->discontinuousCyclicAMIPatch().cyclicAMIPatch()
            );

        vectorField::subField Sf = discPP.faceAreas();
        forAll(Sf, facei)
        {
            Sf[facei] = Sf[facei] / nBranch;
        }

        discPP.areaFraction(1/nBranch);

        const polyMesh& mesh = this->patch().boundaryMesh().mesh();


        primitiveMeshTools::updateCellCentresAndVols
        (
            mesh,
            mesh.faceCentres(),
            mesh.faceAreas(),                      
            uniqueSort(discPP.faceCells()), 
            mesh.cells(),
            const_cast<vectorField&>(mesh.cellCentres()),
            const_cast<scalarField&>(mesh.cellVolumes())
        );

        // // Also modify fvPatch areas

        const discontinuousCyclicAMIFvPatch& discFVP = 
        (
            this->discontinuousCyclicAMIPatch()
        );

        const_cast<vectorField&>(discFVP.Sf()) = discFVP.patch().faceAreas();
        const_cast<vectorField&>(discFVP.Cf()) = discFVP.patch().faceCentres();
        const_cast<scalarField&>(discFVP.magSf()) = mag(discFVP.patch().faceAreas());


        this->patch().boundaryMesh().mesh().V().write();
    }
}


Foam::scaledDiscontinuousCyclicAMIFvPatchVectorField::scaledDiscontinuousCyclicAMIFvPatchVectorField
(
    const scaledDiscontinuousCyclicAMIFvPatchVectorField& ptf
)
:
    discontinuousCyclicAMIFvPatchField<vector>(ptf),
    scalingFactors_(ptf.scalingFactors_),
    sigma_(ptf.sigma_),
    fluxCorrectionFactors_(ptf.fluxCorrectionFactors_),
    connectedPatches_(ptf.connectedPatches_),
    underRelaxation_(ptf.underRelaxation_)
{}


Foam::scaledDiscontinuousCyclicAMIFvPatchVectorField::scaledDiscontinuousCyclicAMIFvPatchVectorField
(
    const scaledDiscontinuousCyclicAMIFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    discontinuousCyclicAMIFvPatchField<vector>(ptf, iF),
    scalingFactors_(ptf.scalingFactors_),
    sigma_(ptf.sigma_),
    fluxCorrectionFactors_(ptf.fluxCorrectionFactors_),
    connectedPatches_(ptf.connectedPatches_),
    underRelaxation_(ptf.underRelaxation_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::scaledDiscontinuousCyclicAMIFvPatchVectorField::updateCoeffs()
{
    
    scalar nBranch = connectedPatches_.size()+1;

    // Update jump to account for differences in velocity

    const scaledDiscontinuousCyclicAMIFvPatchVectorField& nbrPatch =
            refCast<const scaledDiscontinuousCyclicAMIFvPatchVectorField>
            (
                this->neighbourPatchField()
            );

    scalingFactors_ = nbrPatch.massConservationCoeffs()/massConservationCoeffs();

    
    // Now also compute Sigma

    if(connectedPatches_.size()>0)
    {

        // Get patch index from name
        labelList connectedPatchIndeces(connectedPatches_.size());

        forAll(connectedPatches_, i)
        {
            connectedPatchIndeces[i] = this->patch().boundaryMesh().findPatchID(connectedPatches_[i]);
            if(connectedPatchIndeces[i] == -1)
            {
                FatalErrorInFunction
                    << "Cannot find patch " << connectedPatches_[i] << " in mesh "
                    <<this->patch().boundaryMesh().mesh().name() << exit(FatalError);
            } 
        }

        const surfaceScalarField flux  = this->db().lookupObject<surfaceScalarField>("alphaRhoPhi");
        
        scalarField restOfFlux(this->size(),0);

        forAll(connectedPatchIndeces, i) // sum of fluxes in all directions except the neighb patch
        {
            restOfFlux += flux.boundaryField()[connectedPatchIndeces[i]];
        }

        scalarField nbrFlux = flux.boundaryField()[this->discontinuousCyclicAMIPatch().neighbPatch().index()];


        if(correctFluxAndAreas_)
        {
            fluxCorrectionFactors_ = -nbrFlux/(nbrFlux+restOfFlux+VSMALL);
        }

        sigma_ = underRelaxation_*restOfFlux/massConservationCoeffs() + (1-underRelaxation_)*sigma_;

    }
}

Foam::scalarField
Foam::scaledDiscontinuousCyclicAMIFvPatchVectorField::massConservationCoeffs() const
{
    const scalarField& rho=
        this->patch().lookupPatchField<volScalarField, scalar>("thermo:rho");

    scalarField area(this->size(), gSum(this->patch().magSf()));

    if(correctFluxAndAreas_)
    {
        scalar nBranch = connectedPatches_.size()+1;
        area*=nBranch;
    }

    const scalarField& alpha =
        this->patch().lookupPatchField<volScalarField, scalar>("alpha");
    

    return (alpha*rho*area);
}

Foam::scalarField
Foam::scaledDiscontinuousCyclicAMIFvPatchVectorField::fluxCorrectionFactors() const
{
    if(correctFluxAndAreas_)
    {
        scalar nBranch = connectedPatches_.size()+1;

        const surfaceScalarField flux  = this->db().lookupObject<surfaceScalarField>("alphaRhoPhi");
        
        scalarField restOfFlux = sigma_*massConservationCoeffs();
    
        scalarField nbrFlux = flux.boundaryField()[this->discontinuousCyclicAMIPatch().neighbPatch().index()];
        return  (nBranch*nbrFlux/(nbrFlux+restOfFlux+VSMALL));
    }
    else
    {
        return scalarField(this->size(),1);
    }
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::scaledDiscontinuousCyclicAMIFvPatchVectorField::getNeighbourField
(
    const UList<vector>& internalData
) const
{
    // By pass polyPatch to get nbrId. Instead use discontinuousCyclicAMIFvPatch virtual
    // neighbPatch()
    const auto& neighbPatch = this->discontinuousCyclicAMIPatch().neighbPatch();
    const labelUList& nbrFaceCells = neighbPatch.faceCells();

    Field<vector> pnf(internalData, nbrFaceCells);
    Field<vector> defaultValues;

    if (this->discontinuousCyclicAMIPatch().applyLowWeightCorrection())
    {
        defaultValues = Field<vector>(internalData, this->discontinuousCyclicAMIPatch().faceCells());
    }

    tmp<Field<vector>> tpnf = this->discontinuousCyclicAMIPatch().interpolate(pnf, defaultValues);

    if (doTransform())
    {
        transform(tpnf.ref(), forwardT(), tpnf());
    }

    tpnf = tpnf*scalingFactors_; //correct by scaling factors
    tpnf = tpnf - sigma_*this->patch().nf();
    return tpnf;
}



void Foam::scaledDiscontinuousCyclicAMIFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scaledDiscontinuousCyclicAMIFvPatchVectorField::autoMap(m);
}


void Foam::scaledDiscontinuousCyclicAMIFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    scaledDiscontinuousCyclicAMIFvPatchVectorField::rmap(ptf, addr);

    const scaledDiscontinuousCyclicAMIFvPatchVectorField& tiptf =
        refCast<const scaledDiscontinuousCyclicAMIFvPatchVectorField>(ptf);
}


void Foam::scaledDiscontinuousCyclicAMIFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    scalingFactors_.writeEntry("scalingFactors", os);
    sigma_.writeEntry("sigma", os);
    fluxCorrectionFactors_.writeEntry("fluxCorrectionFactors", os);
    os.writeEntry("correctFluxAndAreas", correctFluxAndAreas_);
    connectedPatches_.writeEntry("connectedPatches", os);
    os.writeEntry("underRelaxation", underRelaxation_);
    
    fvPatchField<vector>::writeValueEntry(os);
}




namespace Foam 
{ 
    makePatchTypeField(fvPatchVectorField, scaledDiscontinuousCyclicAMIFvPatchVectorField); 
}

// ************************************************************************* //
