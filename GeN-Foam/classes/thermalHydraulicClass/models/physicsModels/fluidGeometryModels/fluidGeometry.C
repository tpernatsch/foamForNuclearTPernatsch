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
    const fluid& fluid1,
    const fluid& fluid2,
    const structureModel& structure
)
:
    IOdictionary
    (
        IOobject
        (
            "fluidGeometry",
            fluid1.mesh().time().timeName(),
            objReg,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    dispersed_
    (
        (word(this->lookup("dispersed")) == fluid1.name()) ?
        fluid1 :
        fluid2
    ),
    continuous_
    (
        (dispersed_.name() == fluid1.name()) ?
        fluid2 :
        fluid1
    ),
    structure_
    (
        structure
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
            fluid1,
            fluid2
        )
    ),
    structureInterfacialAreaPartition_
    (
        structureInterfacialAreaPartitionModel::New
        (
            objReg,
            this->subDict("structureInterfacialAreaPartitionModel"),
            fluid1,
            fluid2,
            structure
        )
    )
{
    word dispersedName(this->lookup("dispersed"));
    if 
    (
        dispersedName != fluid1.name() 
        and 
        dispersedName != fluid2.name()
    )
    {
        FatalErrorInFunction
            << "Dispersed fluid " << dispersedName
            << " does not match any existing fluid : " << endl
            << "- " << fluid1.name() << endl
            << "- " << fluid2.name() << exit(FatalError);
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

Foam::tmp<Foam::volScalarField> Foam::fluidGeometry::DhDispersed() const
{
    return fluidDiameter_->Dh();
}

Foam::tmp<Foam::volScalarField> Foam::fluidGeometry::iA() const
{
    return fluidInterfacialArea_->iA();
}

Foam::tmp<Foam::volScalarField> Foam::fluidGeometry::frac1() const
{
    return structureInterfacialAreaPartition_->frac1();
}

Foam::tmp<Foam::volScalarField> Foam::fluidGeometry::frac2() const
{
    return structureInterfacialAreaPartition_->frac2();
}


// ************************************************************************* //
