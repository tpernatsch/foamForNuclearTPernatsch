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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class valueType, class baseModel> 
Foam::byRegimeModel<valueType, baseModel>::byRegimeModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    baseModel
    (
        mesh,
        dict,
        objReg
    ),
    regimeMap_
    (
        regimeMapModel::lookupRegimeMap
        (
            mesh,
            (*this).template get<word>("regimeMap")
        )
    ),
    modelObjReg_(mesh.time())
{
    regimeMap_.constructModels<baseModel, fvMesh>
    (
        modelPtrs_,
        mesh,
        *this,
        modelObjReg_
    );
}

template<class valueType, class baseModel> 
template<class firstArgType>
Foam::byRegimeModel<valueType, baseModel>::byRegimeModel
(
    const firstArgType& firstArg,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    baseModel
    (
        firstArg,
        dict,
        objReg
    ),
    regimeMap_
    (
        regimeMapModel::lookupRegimeMap
        (
            firstArg.mesh(),
            (*this).template get<word>("regimeMap")
        )
    ),
    modelObjReg_(firstArg.mesh().time())
{
    regimeMap_.constructModels<baseModel, firstArgType>
    (
        modelPtrs_,
        firstArg,
        *this,
        modelObjReg_
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class valueType, class baseModel>
valueType Foam::byRegimeModel<valueType, baseModel>::value
(
    const label& celli
) const
{
    return 
        valueType
        (
            regimeMap_.interpolateValue<valueType, baseModel>
            (
                modelPtrs_,
                celli
            )
        );
}

// ************************************************************************* //
