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

#include "regime.H"
#include "fluid.H"
#include "structureModel.H"

#include "FFPair.H"
#include "FSPair.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor for non-interpolated regimes
Foam::regime::regime
(
    const fvMesh& mesh,
    const word name,
    const dictionary& dict
)
:   
    objectRegistry
    (
        IOobject
        (
            "regime."+name,
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    dict_(dict),
    name_(name),
    requiresModelCorrection_(false),
    isInterpolated_(false),
    regime1_(nullptr),
    regime2_(nullptr),
    cellField_(mesh.cells().size(), 0.0),
    cellList_(0),
    coeffs1_
    (
        mesh.cells().size(), 0.0
    ),
    coeffs2_
    (
        mesh.cells().size(), 0.0
    )
{
    Info << endl << "Constructing regime: " << name_ << endl;

    const dictionary& dragDict(dict.subDict("dragModels"));
    forAllConstIter
    (
        dictionary,
        dragDict,
        iter
    )
    {
        dragModel::makeInTable
        (
            mesh, 
            *this, 
            iter->keyword(), 
            dragDict, 
            dragModels_
        );
    }

    const dictionary& heatTransferDict
    (
        dict.subDict("heatTransferModels")
    );
    forAllConstIter
    (
        dictionary,
        heatTransferDict,
        iter
    )
    {
        heatTransferModel::makeInTable
        (
            mesh, 
            *this, 
            iter->keyword(), 
            heatTransferDict, 
            heatTransferModels_
        );
    }

    //- If this regime is used in a two-fluid solver, the
    //  fluidGeometry models need to be created
    if ((mesh.lookupClass<FFPair>()).size() != 0)
    {
        const dictionary& fluidGeometryDict
        (
            dict.subDict("fluidGeometry")
        );
        fluidGeometry_.reset
        (
            fluidGeometry::makeModel
            (
                mesh,
                *this,
                fluidGeometryDict
            )
        );

        const dictionary& twoPhaseDragMultiplierDict
        (
            dict.subDict("twoPhaseDragMultiplierModel")
        );
        twoPhaseDragMultiplier_.reset
        (
            twoPhaseDragMultiplierModel::New
            (
                *this,
                twoPhaseDragMultiplierDict,
                mesh
            )
        );
    }
}

//- Constructor for interpolated regimes
Foam::regime::regime
(
    const fvMesh& mesh,
    regime& regime1,
    regime& regime2
)
:   
    objectRegistry
    (
        IOobject
        (
            "regime."+regime1.name()+"."+regime2.name(),
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    dict_(),
    name_(regime1.name()+"."+regime2.name()),
    requiresModelCorrection_(false),
    isInterpolated_(true),
    regime1_
    (
        &regime1
    ),
    regime2_
    (
        &regime2
    ),
    cellField_(mesh.cells().size(), 0.0),
    cellList_(0),
    coeffs1_
    (
        mesh.cells().size(), 0.0
    ),
    coeffs2_
    (
        mesh.cells().size(), 0.0
    )
{
    Info<< endl << "Constructing interpolated regime: " << name_ << endl
        << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regime::initLocalTables
(
    volTensorFieldPtrTable& KdTable, 
    volScalarFieldPtrTable& htcTable
)
{
    if (!isInterpolated_)
    {
        forAllConstIter
        (
            volTensorFieldPtrTable,
            KdTable,
            iter
        )
        {
            word key(iter.key());
            if (!Kds_.found(key))
            {
                Kds_.insert
                (
                    key,
                    autoPtr<volTensorField>
                    (
                        new volTensorField(*KdTable[key])
                    )
                );
            }
            volTensorField& Kd(*Kds_[key]);
            Kd *= 0.0;
            //Kd.correctBoundaryConditions();
        }

        forAllConstIter
        (
            volScalarFieldPtrTable,
            htcTable,
            iter
        )
        {
            word key(iter.key());
            if (!htcs_.found(key))
            {
                htcs_.insert
                (
                    key,
                    autoPtr<volScalarField>
                    (
                        new volScalarField(*htcTable[key])
                    )
                );
            }
            volScalarField& htc(*htcs_[key]);
            htc *= 0.0;
            //htc.correctBoundaryConditions();
        }
    }
}

void Foam::regime::correctDragModels()
{
    if (!isInterpolated_)
    {
        forAllConstIter
        (
            volTensorFieldPtrTable,
            Kds_,
            iter
        )
        {
            word key(iter.key());
            volTensorField& Kd(*Kds_[key]);
            Kd *= 0.0;

            //- Cycles over all dragModels, the ones between a fluid and 
            //  structure are named fluidName.structure.typeName, but the
            //  dragModel.pairName() will always be just fluidName.structure
            forAllConstIter
            (
                dragModelTable,
                dragModels_,
                iter2
            )
            {
                const dragModel& dragModel(iter2());
                if (key == dragModel.pairName()) dragModel.correctKd(Kd);
            }
        }

        //- Correct the twoPhaseMultiplier and apply it to Kd
        if (twoPhaseDragMultiplier_.valid())
        {
            twoPhaseDragMultiplier_->correct();
            twoPhaseDragMultiplier_->correctKdTable(fluidGeometry_(), Kds_);
        }
    }
}


void Foam::regime::correctHeatTransferModels()
{
    if (!isInterpolated_)
    {
        forAllConstIter
        (
            volScalarFieldPtrTable,
            htcs_,
            iter
        )
        {
            word key(iter.key());
            volScalarField& htc(*htcs_[key]);
            htc *= 0.0;

            //- Cycles over all heatTransferModels, the ones between a fluid 
            //  and structure are named fluidName.structure.typeName, but the
            //  dragModel.pairName() will always be just fluidName.structure
            forAllConstIter
            (
                heatTransferModelTable,
                heatTransferModels_,
                iter2
            )
            {
                const heatTransferModel& heatTransferModel(iter2());
                if (key == heatTransferModel.pairName()) 
                {
                    heatTransferModel.correctHtc(htc);
                       
                    //- A momenti sono più coglione di quanto la Madonna
                    //  sia puttana o Maometto stronzo
                    //- Multiply by the fluid-structure contact fraction if 
                    //  this is a two-phase simulation (fluidGeometry_ only 
                    //  valid in two-phase) and, obviously, if this htc is
                    //  between a fluid and the structure
                    if 
                    (
                        fluidGeometry_.valid() 
                    and heatTransferModel.withStructure()
                    )
                        htc *= 
                            fluidGeometry_->fHtc(heatTransferModel.nameBulk());
                }
            }
        }
    }
}

void Foam::regime::correctFluidGeometryModels()
{
    if (fluidGeometry_.valid() and !isInterpolated_)
    {
        fluidGeometry_->correct();
    }
}

void Foam::regime::correctDragTable(volTensorFieldPtrTable& KdTable) const
{
    if (!isInterpolated_)
    {
        forAllIter
        (
            volTensorFieldPtrTable,
            KdTable,
            iter
        )
        {
            word key(iter.key());
            volTensorField& Kd(*KdTable[key]);
            const volTensorField& KdR(*Kds_[key]);
            forAll(cellList_, i)
            {
                label celli(cellList_[i]);
                Kd[celli] = KdR[celli];
            }
        }
    }
    else
    {
        forAllIter
        (
            volTensorFieldPtrTable,
            KdTable,
            iter
        )
        {
            word key(iter.key());
            volTensorField& Kd(*KdTable[key]);
            const volTensorField& KdR1(*(regime1_->Kds()[key]));
            const volTensorField& KdR2(*(regime2_->Kds()[key]));
            forAll(cellList_, i)
            {
                label celli(cellList_[i]);
                Kd[celli] = 
                    coeffs1_[celli]*KdR1[celli] 
                +   coeffs2_[celli]*KdR2[celli];
            }
        }
    }
}

void Foam::regime::correctHeatTransferTable(volScalarFieldPtrTable& htcTable) 
const
{
    if (!isInterpolated_)
    {
        forAllIter
        (
            volScalarFieldPtrTable,
            htcTable,
            iter
        )
        {
            word key(iter.key());
            volScalarField& htc(*htcTable[key]);
            const volScalarField& htcR(*htcs_[key]);
            forAll(cellList_, i)
            {
                label celli(cellList_[i]);
                htc[celli] = htcR[celli];
            }
        }
    }
    else
    {
        forAllIter
        (
            volScalarFieldPtrTable,
            htcTable,
            iter
        )
        {
            word key(iter.key());
            volScalarField& htc(*htcTable[key]);
            const volScalarField& htcR1(*(regime1_->htcs()[key]));
            const volScalarField& htcR2(*(regime2_->htcs()[key]));
            forAll(cellList_, i)
            {
                label celli(cellList_[i]);
                htc[celli] = 
                    coeffs1_[celli]*htcR1[celli] 
                +   coeffs2_[celli]*htcR2[celli];
            }
        }
    }
}

void Foam::regime::correctFluidGeometryFields
(
    volScalarField& iA,
    fluid& fluid1,
    fluid& fluid2,
    FSPair& FSPair1,
    FSPair& FSPair2
) const
{
    autoPtr<volScalarField>& f1Ptr(FSPair1.fPtr());
    autoPtr<volScalarField>& f2Ptr(FSPair2.fPtr());

    if(!f1Ptr.valid())
    {
        f1Ptr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "f."+FSPair1.name(),
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                scalar(0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }
    if(!f2Ptr.valid())
    {
        f2Ptr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "f."+FSPair2.name(),
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                scalar(0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }
    volScalarField& f1(f1Ptr());
    volScalarField& f2(f2Ptr());

    volScalarField& Dh1(fluid1.Dh());
    volScalarField& Dh2(fluid2.Dh());
    scalarField& dispersion1(fluid1.dispersion());
    scalarField& dispersion2(fluid2.dispersion());
    
    if (!isInterpolated_)
    {
        const volScalarField& f1R(fluidGeometry_->fHtc(fluid1.name()));
        const volScalarField& f2R(fluidGeometry_->fHtc(fluid2.name()));
        const volScalarField& DhDispersedR(fluidGeometry_->DhDispersed());
        const volScalarField& DhStructureR(fluidGeometry_->DhStructure());
        const volScalarField& iAR(fluidGeometry_->iA());
        if (fluidGeometry_->nameDispersed() == fluid1.name())
        {
            forAll(cellList_, i)
            {
                label celli(cellList_[i]);
                Dh1[celli] = DhDispersedR[celli];
                Dh2[celli] = DhStructureR[celli];
                dispersion1[celli] = 1.0;
                iA[celli] = iAR[celli];
                f1[celli] = f1R[celli];
                f2[celli] = f2R[celli];
            }
        }
        else
        {
            forAll(cellList_, i)
            {
                label celli(cellList_[i]);
                Dh2[celli] = DhDispersedR[celli];
                Dh1[celli] = DhStructureR[celli];
                dispersion2[celli] = 1.0;
                iA[celli] = iAR[celli];
                f1[celli] = f1R[celli];
                f2[celli] = f2R[celli];
            }
        }
    }
    else //- If regime is interpolated
    {
        const volScalarField& f1R1
        (
            regime1_->fluidGeometry_->fHtc(fluid1.name())
        );
        const volScalarField& f2R1
        (
            regime1_->fluidGeometry_->fHtc(fluid2.name())
        );
        const volScalarField& f1R2
        (
            regime2_->fluidGeometry_->fHtc(fluid1.name())
        );
        const volScalarField& f2R2
        (
            regime2_->fluidGeometry_->fHtc(fluid2.name())
        );
        const volScalarField& DhDispersedR1
        (
            regime1_->fluidGeometry_->DhDispersed()
        );
        const volScalarField& iAR1(regime1_->fluidGeometry_->iA());
        const volScalarField& DhDispersedR2
        (
            regime2_->fluidGeometry_->DhDispersed()
        );
        const volScalarField& iAR2(regime2_->fluidGeometry_->iA());
        const volScalarField& DhStructure
        (
            regime1_->fluidGeometry_->DhStructure()
        );

        bool fluid1DispersedInR1
        (
            regime1_->fluidGeometryRef().nameDispersed() == fluid1.name()
        );
        bool fluid1DispersedInR2
        (
            regime2_->fluidGeometryRef().nameDispersed() == fluid1.name()
        );

        //- If fluid1 is not dispersed in a regime, then fluid2 is
        bool fluid2DispersedInR1
        (
            !fluid1DispersedInR1
        );
        bool fluid2DispersedInR2
        (
            !fluid1DispersedInR2
        );
    
        forAll(cellList_, i)
        {
            label celli(cellList_[i]);
            const scalar& c1(coeffs1_[celli]);
            const scalar& c2(coeffs2_[celli]);
            const scalar& DhStructurei(DhStructure[celli]);
            const scalar& DhDispersedR1i(DhDispersedR1[celli]);
            const scalar& DhDispersedR2i(DhDispersedR2[celli]);
            Dh1[celli] =
                c1*((fluid1DispersedInR1) ? DhDispersedR1i : DhStructurei)
            +   c2*((fluid1DispersedInR2) ? DhDispersedR2i : DhStructurei);
            Dh2[celli] =
                c1*((fluid2DispersedInR1) ? DhDispersedR1i : DhStructurei)
            +   c2*((fluid2DispersedInR2) ? DhDispersedR2i : DhStructurei);
            iA[celli] = c1*iAR1[celli] + c2*iAR2[celli];
            dispersion1[celli] = 
                c1*(fluid1DispersedInR1) + c2*(fluid1DispersedInR2);
            dispersion2[celli] = 
                c1*(fluid2DispersedInR1) + c2*(fluid2DispersedInR2);
            f1[celli] = c1*f1R1[celli] + c2*f1R2[celli];
            f2[celli] = c1*f2R1[celli] + c2*f2R2[celli];
        }
    }
}


// ************************************************************************* //
