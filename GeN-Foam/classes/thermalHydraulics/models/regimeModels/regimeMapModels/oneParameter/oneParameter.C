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
#include "myOps.H"


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

const Foam::Enum
<
    Foam::regimeMapModels::oneParameter::interpolationMode
>
Foam::regimeMapModels::oneParameter::interpolationModeNames_
(
    {
        { 
            interpolationMode::linear, 
            "linear" 
        },
        { 
            interpolationMode::quadratic, 
            "quadratic" 
        }
    }
);

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
    thresholds_(0),
    interpolationMode_
    (
        interpolationModeNames_.get
        (
            this->lookupOrDefault<word>
            (
                "interpolationMode", 
                "linear"
            )
        )
    )
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

    //- Re-arrange individual bounds so that leftmost value is smaller
    //  than the rightmost one
    for (int i = 0; i < bounds.size(); i++)
    {
        if (bounds[i][0] > bounds[i][1])
        {
            scalar tmp = bounds[i][0];
            bounds[i][0] = bounds[i][1];
            bounds[i][1] = tmp;
        }
    }
    
    //- Then do the re-arranging as described before
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
                            regimes_[regime1Name](),
                            regimes_[regime2Name]()
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

    //- Extend lower and upper bounds to infinity
    thresholds_[0] = -1e69;
    thresholds_[thresholds_.size()-1] = 1e69;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regimeMapModels::oneParameter::correct()
{
    forAll(orderedRegimeNames_, i)
    {
        regime& regime = regimes_[orderedRegimeNames_[i]]();
        labelList& cellList(regime.cellList());
        cellList = labelList(0);
        const scalar& t0(thresholds_[i]);
        const scalar& t1(thresholds_[i+1]);
        scalar dt(t1-t0);
        scalarField t1_pByDt((t1-parameter_)/dt);
        scalarField p_t0ByDt((parameter_-t0)/dt);
        bool isCurrentlyPresent = false;//- I can't use the omonymous method
                                        //  of the regime class as the len of
                                        //  its cellList is still 0 at this
                                        //  point
        //- This is equivalent to 
        //  isPresent = ((max(t1_pByDt) > 0) and max(p_t0ByDt[j]) >= 0)); 
        //  but faster
        forAll(mesh_.cells(), j)
        {
            if (t1_pByDt[j] > 0  and p_t0ByDt[j] >= 0)
            {
                isCurrentlyPresent = true;
                break;
            }
        }

        scalar f(myOps::relaxationFactor(mesh_, "regime"));
        scalarField& cellField(regime.cellField());
        if (isCurrentlyPresent)
        {
            cellField = pos0(p_t0ByDt)*pos(t1_pByDt);
            if (regime.isInterpolated())
            {
                switch (interpolationMode_)
                {
                    case interpolationMode::linear :
                    {
                        if (f == 1.0)
                        {
                            forAll(mesh_.cells(), j)
                            {
                                if (cellField[j] == 1)
                                {
                                    cellList.append(j);
                                    regime.coeffs1()[j] = t1_pByDt[j];
                                    regime.coeffs2()[j] = p_t0ByDt[j];
                                }
                            }
                        }
                        else
                        {
                            scalarField& C1(regime.coeffs1());
                            scalarField& C2(regime.coeffs2());
                            forAll(mesh_.cells(), j)
                            {
                                if (cellField[j] == 1)
                                {
                                    cellList.append(j);
                                    scalar& c1(C1[j]);
                                    scalar& c2(C2[j]);
                                    c1 = f*t1_pByDt[j] + (1.0-f)*c1;
                                    c2 = 1.0-c1;
                                }
                            }
                        }
                        break;
                    }
                    case interpolationMode::quadratic :
                    {
                        scalar mid(t0+dt/2.0);
                        if (f == 1.0)
                        {
                            forAll(mesh_.cells(), j)
                            {
                                if (cellField[j] == 1)
                                {
                                    cellList.append(j);
                                    scalar& c1(regime.coeffs1()[j]);
                                    scalar& c2(regime.coeffs2()[j]);
                                    c1 = 
                                    (
                                        (parameter_[j] <= mid) ?
                                        1.0-2.0*sqr(p_t0ByDt[j]) :
                                        2.0*sqr(t1_pByDt[j])
                                    );
                                    c2 = 1.0-c1;
                                }
                            }
                        }
                        else
                        {
                            forAll(mesh_.cells(), j)
                            {
                                if (cellField[j] == 1)
                                {
                                    cellList.append(j);
                                    scalar& c1(regime.coeffs1()[j]);
                                    scalar& c2(regime.coeffs2()[j]);
                                    c1 = 
                                    (
                                        f*
                                        (
                                            (parameter_[j] <= mid) ?
                                            1.0-2.0*sqr(p_t0ByDt[j]) :
                                            2.0*sqr(t1_pByDt[j])
                                        ) 
                                    +   (1.0-f)*c1
                                    );
                                    c2 = 1.0-c1;
                                }
                            }
                        }
                        break;
                    }
                }
            }
            else
            {
                forAll(mesh_.cells(), j)
                {
                    if (cellField[j] == 1)
                    {
                        cellList.append(j);
                    }
                }
            }
        }
        else
        {
            cellField *= 0;
        }
    }

    //- Update the requiresModelCorrection flag. The way this is done is
    //  not related to the specific run-time-selectable regimeMapModel,
    //  so its wrapped in a function in the base class
    this->setRequiresModelCorrection();
}


// ************************************************************************* //
