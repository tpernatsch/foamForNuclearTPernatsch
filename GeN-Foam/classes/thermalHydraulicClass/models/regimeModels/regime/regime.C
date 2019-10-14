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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor for non-interpolated regimes
Foam::regime::regime
(
    const fvMesh& mesh,
    const word name,
    const dictionary& dict
)
:   
    scalarField
    (
        mesh.cells().size(),
        scalar(0)
    ),
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
    isCurrentlyPresent_(false),
    isInterpolated_(false),
    regime1_(*this),
    regime2_(*this)
{
    Info << endl << "Constructing regime: " << name_ << endl;

    /*
    typedef HashTable<const fluid*, word, word::hash> fluidTable;
    fluidTable fluids = mesh.lookupClass<fluid>();
    wordList fluidNames;
    forAllConstIter
    (
        fluidTable,
        fluids,
        iter
    )
    {
        fluidNames.append((*iter)->name());
    }
    Info << fluidNames << endl;
    wordList regionNames
    (
        mesh.lookupObject<structureModel>("alpha.structure").regions()
    );
    */

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
    }
}

//- Constructor for interpolated regimes
Foam::regime::regime
(
    const fvMesh& mesh,
    const regime& regime1,
    const regime& regime2
)
:   
    scalarField
    (
        mesh.cells().size(),
        scalar(0)
    ),
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
    isCurrentlyPresent_(false),
    isInterpolated_(true),
    regime1_
    (
        regime1
    ),
    regime2_
    (
        regime2
    ),
    coeffs1_
    (
        mesh.cells().size()
    ),
    coeffs2_
    (
        mesh.cells().size()
    )
{
    Info<< endl << "Constructing interpolated regime: " << name_ << endl
        << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regime::correct()
{
    isCurrentlyPresent_ = (max(*this) != 0);
}

void Foam::regime::correctDragTable(volTensorFieldTable& KdTable) const
{
    if (!isInterpolated_)
    {
        forAllConstIter
        (
            dragModelTable,
            dragModels_,
            dragModelIter
        )
        {
            const dragModel& dragModel(dragModelIter());
            word pairNameInTable
            (
                (dragModel.withStructure()) ? 
                word(dragModel.name1()+".structure") :
                dragModel.pairName()
            );
            volTensorField& Kd(KdTable[pairNameInTable]);
            Kd.primitiveFieldRef() +=
                dragModel.Kd()().primitiveField()*(*this);
            Kd.correctBoundaryConditions();
        }
    }
    else
    {
        forAllConstIter
        (
            dragModelTable,
            regime1_.dragModels(),
            dragModelIter
        )
        {
            const dragModel& dragModel1
            (
                dragModelIter()
            );
            const dragModel& dragModel2
            (
                regime2_.dragModels()[dragModel1.pairName()]
            );
            word pairNameInTable
            (
                (dragModel1.withStructure()) ? 
                word(dragModel1.name1()+".structure") :
                dragModel1.pairName()
            );
            volTensorField& Kd(KdTable[pairNameInTable]);
            Kd.primitiveFieldRef() +=
                dragModel1.Kd()().primitiveField()*coeffs1_ +
                dragModel2.Kd()().primitiveField()*coeffs2_;
            Kd.correctBoundaryConditions();
        }
    }
}

void Foam::regime::correctHeatTransferTable(volScalarFieldTable& htcTable) const
{
    if (!isInterpolated_)
    {
        forAllConstIter
        (
            heatTransferModelTable,
            heatTransferModels_,
            heatTransferModelIter
        )
        {
            const heatTransferModel& heatTransferModel
            (
                heatTransferModelIter()
            );
            word pairNameInTable
            (
                (heatTransferModel.withStructure()) ? 
                word(heatTransferModel.nameBulk()+".structure") :
                heatTransferModel.pairName()
            );
            volScalarField& htc(htcTable[pairNameInTable]);
            htc.primitiveFieldRef() +=
                heatTransferModel.htc()().primitiveField()*(*this);
            htc.correctBoundaryConditions();
        }
    }
    else
    {
        forAllConstIter
        (
            heatTransferModelTable,
            regime1_.heatTransferModels(),
            heatTransferModelIter
        )
        {
            const heatTransferModel& heatTransferModel1
            (
                heatTransferModelIter()
            );
            const heatTransferModel& heatTransferModel2
            (
                regime2_.heatTransferModels()[heatTransferModel1.pairName()]
            );
            word pairNameInTable
            (
                (heatTransferModel1.withStructure()) ? 
                word(heatTransferModel1.nameBulk()+".structure") :
                heatTransferModel1.pairName()
            );
            volScalarField& htc(htcTable[pairNameInTable]);
            htc.primitiveFieldRef() +=
                heatTransferModel1.htc()().primitiveField()*coeffs1_ +
                heatTransferModel2.htc()().primitiveField()*coeffs2_;
            htc.correctBoundaryConditions();
        }
    }
}

void Foam::regime::correctFluidGeometry
(
    volScalarField& iA,
    volScalarField& frac1,
    volScalarField& frac2,
    fluid& fluid1,
    fluid& fluid2
) const
{
    volScalarField& Dh1(fluid1.Dh());
    volScalarField& Dh2(fluid2.Dh());
    scalarField& dispersion1(fluid1.dispersion());
    scalarField& dispersion2(fluid2.dispersion());
    
    if (!isInterpolated_)
    {
        if 
        (
            fluidGeometry_->nameDispersed()
            == 
            fluid1.name()
        )
        {
            Dh1.primitiveFieldRef() += 
                fluidGeometry_->DhDispersed()().primitiveField()*(*this);
            Dh2.primitiveFieldRef() +=
                fluidGeometry_->DhStructure().primitiveField()*(*this);
            dispersion1 += *this;
        }
        else
        {
            Dh2.primitiveFieldRef() += 
                fluidGeometry_->DhDispersed()().primitiveField()*(*this);
            Dh1.primitiveFieldRef() +=
                fluidGeometry_->DhStructure().primitiveField()*(*this);
            dispersion2 += *this;
        }

        iA.primitiveFieldRef() +=
            fluidGeometry_->iA()().primitiveField()*(*this);
        frac1.primitiveFieldRef() +=
            fluidGeometry_->frac1()().primitiveField()*(*this);
        frac2.primitiveFieldRef() +=
            fluidGeometry_->frac2()().primitiveField()*(*this);
    }
    else
    {
        volScalarField r1Dh1(Dh1*0.0);
        volScalarField r2Dh1(r1Dh1);
        volScalarField r1Dh2(r1Dh1);
        volScalarField r2Dh2(r1Dh1);
        
        if 
        (
            regime1_.fluidGeometryRef().nameDispersed()
            == 
            fluid1.name()
        )
        {
            r1Dh1.primitiveFieldRef() += 
                    regime1_.fluidGeometry_->DhDispersed()().primitiveField()
                *   (*this);
            r1Dh2.primitiveFieldRef() +=
                    regime1_.fluidGeometry_->DhStructure().primitiveField()
                *   (*this);
            dispersion1 += coeffs1_;
        }
        else
        {
            r1Dh2.primitiveFieldRef() += 
                    regime1_.fluidGeometry_->DhDispersed()().primitiveField()
                *   (*this);
            r1Dh1.primitiveFieldRef() +=
                    regime1_.fluidGeometry_->DhStructure().primitiveField()
                *   (*this);
            dispersion2 += coeffs1_;
        }
        if 
        (
            "Dh."+regime2_.fluidGeometryRef().nameDispersed()
            == 
            Dh1.name()
        )
        {
            r2Dh1.primitiveFieldRef() += 
                    regime2_.fluidGeometry_->DhDispersed()().primitiveField()
                *   (*this);
            r2Dh2.primitiveFieldRef() +=
                    regime2_.fluidGeometry_->DhStructure().primitiveField()
                *   (*this);
            dispersion1 += coeffs2_;
        }
        else
        {
            r2Dh2.primitiveFieldRef() += 
                    regime2_.fluidGeometry_->DhDispersed()().primitiveField()   
                *   (*this);
            r2Dh1.primitiveFieldRef() +=
                    regime2_.fluidGeometry_->DhStructure().primitiveField()
                *   (*this);
            dispersion2 += coeffs2_;
        }
        Dh1.primitiveFieldRef() +=
            coeffs1_*r1Dh1.primitiveField() + coeffs2_*r2Dh1.primitiveField();
        Dh2.primitiveFieldRef() +=
            coeffs1_*r1Dh2.primitiveField() + coeffs2_*r2Dh2.primitiveField();

        iA.primitiveFieldRef() +=
                coeffs1_*regime1_.fluidGeometryRef().iA()().primitiveField()
            +   coeffs2_*regime2_.fluidGeometryRef().iA()().primitiveField();

        frac1.primitiveFieldRef() +=
                coeffs1_*regime1_.fluidGeometryRef().frac1()().primitiveField()
            +   coeffs2_*regime2_.fluidGeometryRef().frac1()().primitiveField()
            ;
        frac2.primitiveFieldRef() +=
                coeffs1_*regime1_.fluidGeometryRef().frac2()().primitiveField()
            +   coeffs2_*regime2_.fluidGeometryRef().frac2()().primitiveField()
            ;
    }

    Dh1.correctBoundaryConditions();
    Dh2.correctBoundaryConditions();
    iA.correctBoundaryConditions();
    frac1.correctBoundaryConditions();
    frac2.correctBoundaryConditions();
}


// ************************************************************************* //