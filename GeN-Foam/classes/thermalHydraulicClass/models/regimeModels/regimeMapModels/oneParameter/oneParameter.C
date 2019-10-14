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

#include "oneParameter.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regimeMapModels
{
    defineTypeNameAndDebug(oneParameter, 0);
    addToRunTimeSelectionTable
    (
        regimeMapModel, 
        oneParameter, 
        regimeMapModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regimeMapModels::oneParameter::oneParameter
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& physicsModelsDict
)
:
    regimeMapModel
    (
        mesh,
        dict,
        physicsModelsDict
    ),
    parameter_
    (
        mesh_.lookupObjectRef<volScalarField>
        
        (
            word(this->lookup("parameter"))
        ).field()
    ),
    thresholds_(0)
{
    //- Read the regimeNames and the bounds from the dict
    wordList unorderedRegimeNames(0);
    List<List<scalar>> bounds(0);
    forAllIter
    (
        dictionary,
        this->subDict("regimeBounds"),
        iter
    )
    {   
        unorderedRegimeNames.append(iter->keyword());  
        bounds.append
        (
            List<scalar>(iter->stream())
        );
    }

    //- Order bounds and names into ascending threshold values, and
    //  re-arrange the data so that something in a format e.g.:
    /*
            regimeBounds
            {
                "regime0" (0.5  1);
                "regime1" (1.2  1.5);
                "regime2" (0    0.37);
                "regime4" (1    1.1);
            }

        is converted into two ordered lists, ordered by threshold value, i.e.:
        
            thresholds_ = [0, 0.37, 0.5, 1, 1.1, 1.2, 1.5]
            ordredRegimeNames_ = 
            [
                "regime2", 
                "regime2.regime0",
                "regime0",
                "regime4",
                "regime4.regime1",
                "regime1"
            ]
        so that orderedRegimeNames_[i] is bounded by thresholds_[i] and
        thresholds_[i+1]. If the upper threshold of regime i and the lower
        threshold of regime i+1 mismatch, an interpolation regime (e.g. 
        "regime2.regime0") is created. 
    */
    List<scalar> tmpBound;
    word tmpName;
    bool swap(true);
    while (swap)
    {
        swap = false;
        for (int i = 0; i < bounds.size()-1; i++)
        {
            if (bounds[i][0] > bounds[i+1][0])
            {
                tmpBound = bounds[i];
                bounds[i] = bounds[i+1];
                bounds[i+1] = tmpBound;
                tmpName = unorderedRegimeNames[i];
                unorderedRegimeNames[i] = unorderedRegimeNames[i+1];
                unorderedRegimeNames[i+1] = tmpName;
                swap = true;
                break;
            }
        }
    }
    //- After the end of the while above, unorderedRegimeNames are ordered now
    wordList orderedNonInterpolatedRegimeNames(unorderedRegimeNames);
    
    //- Construct the regimeTable, non-interpolated regimes first
    forAll(orderedNonInterpolatedRegimeNames, i)
    {
        word regimeName(orderedNonInterpolatedRegimeNames[i]);
        regimes_.set
        (
            regimeName,
            autoPtr<regime>
            (
                new regime
                (
                    mesh_,
                    regimeName,
                    this->physicsModelsDict_.subDict(regimeName)
                )
            )
        );
    }

    //- Construct interpolated regimes and add them to the regimeTable
    thresholds_.append(bounds[0][0]);
    orderedRegimeNames_.append(orderedNonInterpolatedRegimeNames[0]);
    for (int i = 1; i < bounds.size(); i++)
    {
        if (bounds[i-1][1] != bounds[i][0])
        {
            //- Then this regime is interpolated
            word regime1Name(orderedNonInterpolatedRegimeNames[i-1]);
            word regime2Name(orderedNonInterpolatedRegimeNames[i]);
            word interpolatedRegimeName(regime1Name+"."+regime2Name);
            if (!regimes_.found(interpolatedRegimeName)) //- Just in case
            {
                regimes_.set
                (
                    interpolatedRegimeName,
                    autoPtr<regime>
                    (
                        new regime
                        (
                            mesh_,
                            regimes_[regime1Name],
                            regimes_[regime2Name]
                        )
                    )
                );
            }
            orderedRegimeNames_.append(interpolatedRegimeName);
            thresholds_.append(bounds[i-1][1]);
        }
        thresholds_.append(bounds[i][0]);
        orderedRegimeNames_.append(orderedNonInterpolatedRegimeNames[i]);
    }
    thresholds_.append(bounds[bounds.size()-1][1]);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regimeMapModels::oneParameter::correct()
{
    //- Adjust the marker field of each regime. Furthermore, the lowest
    //  threshold regime and highest threshold regime are extended to
    //  plus and minus infinity respectively and the regimes adjusted
    //  accordingly
    label n(orderedRegimeNames_.size());
    regime& regime0(regimes_[orderedRegimeNames_[0]]());
    regime0 = neg(parameter_-thresholds_[1]);
    regime& regimen(regimes_[orderedRegimeNames_[n-1]]());
    regimen = pos0(parameter_-thresholds_[n-1]);
    for (int i = 0; i < n-1; i++)
    {
        regime& regime(regimes_[orderedRegimeNames_[i]]());
        scalar t0(thresholds_[i]);
        scalar t1(thresholds_[i+1]);
        regime = pos0(parameter_-t0)*neg(parameter_-t1);
    }

    //- Compute the interpolation coefficients for the interpolated regimes
    forAll(orderedRegimeNames_, i)
    {
        regime& regime(regimes_[orderedRegimeNames_[i]]());
        if (regime.isInterpolated())
        {
            scalar t0(thresholds_[i]);
            scalar t1(thresholds_[i+1]);
            regime.coeffs1() = regime*(t1-parameter_)/(t1-t0);
            regime.coeffs2() = regime*(parameter_-t0)/(t1-t0);
        }

        //- Correct the isCurrentlyPresent flag of each regime
        regime.correct();
    }
}


// ************************************************************************* //
