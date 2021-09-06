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

#include "regimeMapModel.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class modelType, class firstArgType>
void Foam::regimeMapModel::constructModels
(
    List<autoPtr<modelType>>& models,
    const firstArgType& firstArg,
    const dictionary& dict,
    const objectRegistry& objReg
) const
{
    models = List<autoPtr<modelType>>(0);
    forAll(regimeLabelToName_, i)
    {
        word regimeName(regimeLabelToName_[i]);
        const dictionary& regimeDict(dict.subDict(regimeName));
        models.append
        (
            modelType::New
            (
                firstArg,
                regimeDict,
                objReg
            )
        );
    }
}

template<class valueType, class modelType>
valueType Foam::regimeMapModel::interpolateValue
(
    const List<autoPtr<modelType>>& models,
    const label& celli
) const
{
    //- The If-else appears to be useless, but it is to force a Return Value
    //  Optimization for the case n==1. I should check the actual performance 
    //  gain...
    const DynamicList<Tuple2<label,scalar>>& rlci
    (
        regimeLabelCoeffs_[celli]
    );
    int n(rlci.size());
    if (n == 1)
        return valueType(models[rlci[0].first()]->value(celli));
    //- else
    const Tuple2<label,scalar>& rlci0(rlci[0]);
    valueType value
    (
        rlci0.second()*models[rlci0.first()]->value(celli)
    );
    for(int j=1; j<n; j++)
    {
        //- The first element of the tuple is the regime label (which 
        //  corresponds to the label of the sub-model in models
        //  while the second element of the tuple corresponds to 
        //  the coefficient of the associated regime, used to weight the
        //  contribution of the correpsonding regime model value valuej in 
        //  the total value value
        const Tuple2<label,scalar>& rlcij(rlci[j]);
        value += rlcij.second()*models[rlcij.first()]->value(celli);
    }
    return value;
}
