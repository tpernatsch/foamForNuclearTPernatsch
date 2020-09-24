/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "fluidGeometry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidGeometry::fluidGeometry
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const fluid& dispersed,
    const fluid& continuous,
    const structureModel& structure
)
:
    IOdictionary
    (
        IOobject
        (
            typeName,
            dispersed.mesh().time().timeName(),
            objReg,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    dispersed_
    (
        dispersed
    ),
    continuous_
    (
        continuous
    ),
    structure_
    (
        structure
    ),
    DhDispersed_
    (
        IOobject
        (
            "DhDispersed",
            dispersed.mesh().time().timeName(),
            objReg
        ),
        dispersed.mesh(),
        dimensionedScalar("", dimLength, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    iA_
    (
        IOobject
        (
            "fluidInterfacialArea",
            dispersed.mesh().time().timeName(),
            objReg
        ),
        dispersed.mesh(),
        dimensionedScalar("", dimArea/dimVolume, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    fKdD_
    (
        IOobject
        (
            "fKd."+dispersed.name(),
            dispersed.mesh().time().timeName(),
            objReg
        ),
        dispersed.mesh(),
        dimensionedScalar("", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    fKdC_
    (
        IOobject
        (
            "fKd."+continuous.name(),
            continuous.mesh().time().timeName(),
            objReg
        ),
        continuous.mesh(),
        dimensionedScalar("", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    fHtcD_
    (
        IOobject
        (
            "fHtc."+dispersed.name(),
            dispersed.mesh().time().timeName(),
            objReg
        ),
        dispersed.mesh(),
        dimensionedScalar("", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    fHtcC_
    (
        IOobject
        (
            "fHtc."+continuous.name(),
            continuous.mesh().time().timeName(),
            objReg
        ),
        continuous.mesh(),
        dimensionedScalar("", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    fluidDiameter_
    (
        fluidDiameterModel::New 
        (
            objReg,
            this->subDict("dispersedFluidDiameterModel"),
            dispersed_
        )
    ),
    fluidInterfacialArea_
    (
        fluidInterfacialAreaModel::New
        (
            objReg,
            this->subDict("fluidInterfacialAreaModel"),
            dispersed_,
            continuous_
        )
    )
{
    //- If the same partition model is to be used for both drag and heat
    //  and heat transfer, specify it via a single fluidPartitionModel 
    //  keyword for simplicity
    if (this->isDict("fluidPartitionModel"))
    {
        fluidPartition_.reset
        (
            fluidPartitionModel::New
            (
                objReg,
                this->subDict("fluidPartitionModel"),
                dispersed_,
                continuous_,
                structure
            )
        );
    }
    else
    {
        dragPartition_.reset
        (
            fluidPartitionModel::New
            (
                objReg,
                this->subDict("dragPartitionModel"),
                dispersed_,
                continuous_,
                structure
            )
        );

        heatTransferPartition_.reset
        (
            fluidPartitionModel::New
            (
                objReg,
                this->subDict("heatTransferPartitionModel"),
                dispersed_,
                continuous_,
                structure
            )
        );
    }
}


// * * * * * * * * * * * * * * * Static Member Functions   * * * * * * * * * //

Foam::autoPtr<Foam::fluidGeometry> Foam::fluidGeometry::makeModel
(
    const objectRegistry& srcObjReg,
    const objectRegistry& dstObjReg,
    const dictionary& dict
)
{
    word continuousName;
    word dispersedName(dict.lookup("dispersed"));
    typedef HashTable<const fluid*, word, word::hash> fluidTable;
    fluidTable fluids(srcObjReg.lookupClass<fluid>());
    forAllConstIter
    (
        fluidTable,
        fluids,
        iter
    )
    {
        if ((*iter)->name() != dispersedName)
        {
            continuousName = (*iter)->name();
        } 
    }
    const fluid& dispersed
    (
        srcObjReg.lookupObject<fluid>("alpha."+dispersedName)
    );
    const fluid& continuous
    (
        srcObjReg.lookupObject<fluid>("alpha."+continuousName)
    );

    const structureModel& structure
    (
        srcObjReg.lookupObject<structureModel>("alpha.structure")
    );

    return 
        autoPtr<fluidGeometry>
        (
            new fluidGeometry
            (
                dstObjReg,
                dict,
                dispersed,
                continuous,
                structure
            )
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluidGeometry::correct()
{
    //- It is important to update the DhDispersed BEFORE the iA_, as currently
    //  some interfacial area models depend on the DhDispersed (e.g. 
    //  sphericalTopology)
    DhDispersed_= fluidDiameter_->Dh()();
    iA_ = fluidInterfacialArea_->iA()();
    if (fluidPartition_.valid())
    {
        volScalarField fD(fluidPartition_->fracD()());
        volScalarField fC(fluidPartition_->fracC()());
        fKdD_ = fD;
        fKdC_ = fC;
        fHtcD_ = fD;
        fHtcC_ = fC;
    }
    else
    {
        fKdD_ = dragPartition_->fracD()();
        fKdC_ = dragPartition_->fracC()();
        fHtcD_ = heatTransferPartition_->fracD()();
        fHtcC_ = heatTransferPartition_->fracC()();
    }
}
    
const Foam::volScalarField& Foam::fluidGeometry::fKd(const word& name) const
{
    if (word("fKd."+name) == fKdD_.name())
    {

        Info << "A, Asked for " << name << " returning " << max(fKdD_).value() << endl;
        return fKdD_;
    }
    else if (word("fKd."+name) == fKdC_.name())
    {
        Info << "B, Asked for " << name << " returning " << max(fKdC_).value() << endl;
        return fKdC_;
    }
    FatalErrorInFunction
    << "Requested fluid " << name
    << " does not match any existing fluid!" << exit(FatalError);
    return fKdD_;
}

const Foam::volScalarField& Foam::fluidGeometry::fHtc(const word& name) const
{
    if (word("fHtc."+name) == fHtcD_.name())
        return fHtcD_;
    else if (word("fHtc."+name) == fHtcC_.name())
        return fHtcC_;
    FatalErrorInFunction
    << "Requested fluid " << name
    << " does not match any existing fluid!" << exit(FatalError);
    return fHtcD_;
}

// ************************************************************************* //
