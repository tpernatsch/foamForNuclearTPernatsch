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

#include "nucleateBoiling.H"
#include "addToRunTimeSelectionTable.H"
#include "myStringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(nucleateBoiling, 0);
    addToRunTimeSelectionTable
    (
        heatTransferModel, 
        nucleateBoiling, 
        FSHeatTransferModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::nucleateBoiling::
nucleateBoiling
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FSPair& FSPair,
    const wordList& regions
)
:
    heatTransferModel
    (
        objReg,
        dict,
        FSPair,
        regions
    ),
    flowFactor_
    (
        flowFactorModel::New
        (
            *this,
            objReg,
            FSPair
        )
    ),
    suppressionFactor_
    (
        suppressionFactorModel::New
        (
            *this,
            objReg,
            FSPair
        )
    ),
    convectionHeatTransfer_
    (
        heatTransferModel::New
        (
            objReg,
            dict.subDict("convectionModel"),
            FSPair,
            //- Use the regions list obtained from the top level  
            //  heatTransferModel dictionary name via myStringOps
            myStringOps::split<word>(this->dictName(), ':')
        )
    ),
    poolBoilingHeatTransfer_
    (
        poolBoilingModel::New
        (
            *this,
            objReg,
            FSPair
        )
    ),
    F_
    (
        IOobject
        (
            "flowFactor."+objReg.name(),
            FSPair.fluidRef().mesh().time().timeName(),
            FSPair.fluidRef().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        FSPair.fluidRef().mesh(),
        dimensionedScalar("", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    S_
    (
        IOobject
        (
            "suppressionFactor."+objReg.name(),
            FSPair.fluidRef().mesh().time().timeName(),
            FSPair.fluidRef().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        FSPair.fluidRef().mesh(),
        dimensionedScalar("", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    dmdt_
    (
        mesh_.lookupObject<volScalarField>("dmdt")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heatTransferModels::nucleateBoiling::correctHtc
(
    volScalarField& htc
) const
{    
    //- Correct htc by first computing its convective component. Now htc
    //  consists only of the purely convective component
    convectionHeatTransfer_->correctHtc(htc);

    //- Only correct nucleate boiling models if the normalized vapour alpha is
    //  larger than 1.0-alpha0
    //scalar alpha0(0.99999);
    //const volScalarField& alphaNorm(FSPair_->fluidRef().normalized());
    const scalarField& phaseChangeSign(FSPair_->fluidRef().phaseChangeSign());
    if (min(phaseChangeSign) < 0.0)
    {
        //- Correct the convection enhancement flow factor
        flowFactor_->correct();
        const scalarField& F(flowFactor_->flowFactor());

        //- Correct the nucleate boiling suppression factor
        suppressionFactor_->correct();
        const scalarField& S(suppressionFactor_->suppressionFactor());

        forAll(cellList_, i)
        {
            const label& celli(cellList_[i]);
            F_[celli] = F[i];
            S_[celli] = S[i];
        }
        F_.correctBoundaryConditions();
        S_.correctBoundaryConditions();

        //- Correct the pool boiling heat transfer coefficient
        poolBoilingHeatTransfer_->correct();
        const scalarField& htcPB(poolBoilingHeatTransfer_->htc());

        forAll(cellList_, i)
        {
            const label& celli(cellList_[i]);

            //- Correct only in those cells where I am boiling
            //- Recall that htc[celli] on the right hand side consists only of
            //  the convective component before the assignment is done
            if (phaseChangeSign[celli] < 0.0)
            {
                scalar& htci(htc[celli]);
                htci = F[i]*htci + S[i]*htcPB[i];
                /*
                Info<< celli << " " << htc[celli] << " " << htcC << " " 
                    << htcPB[i] << " " << min(F[i], 50) << " " << S[i] 
                    << endl;
                */
            }
        }
    }    
}


// ************************************************************************* //
