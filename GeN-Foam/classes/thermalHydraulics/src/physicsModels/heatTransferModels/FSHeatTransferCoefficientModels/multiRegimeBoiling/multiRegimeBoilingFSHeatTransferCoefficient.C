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

#include "FSPair.H"
#include "FFPair.H"
#include "multiRegimeBoilingFSHeatTransferCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FSHeatTransferCoefficientModels
{
    defineTypeNameAndDebug(multiRegimeBoiling, 0);
    addToRunTimeSelectionTable
    (
        FSHeatTransferCoefficientModel, 
        multiRegimeBoiling, 
        FSHeatTransferCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSHeatTransferCoefficientModels::multiRegimeBoiling::
multiRegimeBoiling
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    FSHeatTransferCoefficientModel
    (
        pair,
        dict,
        objReg
    ),
    Tw_(pair.structureRef().Twall()),
    Tf_(pair.fluidRef().thermo().T()),
    FFPairPtr_(nullptr),
    htcFCPtr_
    (
        FSHeatTransferCoefficientModel::New
        (
            pair,
            this->subDict("forcedConvectionModel"),
            pair.mesh()
        )
    ),
    htcPBPtr_
    (
        FSHeatTransferCoefficientModel::New
        (
            pair,
            this->subDict("poolBoilingModel"),
            pair.mesh()
        )
    ),
    FPtr_
    (
        flowEnhancementFactorModel::New
        (
            pair,
            this->subDict("flowEnhancementFactorModel"),
            pair.mesh()
        )
    ),
    SPtr_
    (
        suppressionFactorModel::New
        (
            pair,
            this->subDict("suppressionFactorModel"),
            pair.mesh()
        )
    )
{
    if (this->found("nucleateBoilingOnsetModel"))
    {
        TONBPtr_.reset
        (
            TONBModel::New
            (
                pair,
                this->subDict("nucleateBoilingOnsetModel"),
                pair.mesh()
            )
        );
    }
    if (this->found("subCooledBoilingFractionModel"))
    {
        SCBFPtr_.reset
        (
            subCooledBoilingFractionModel::New
            (
                pair,
                this->subDict("subCooledBoilingFractionModel"),
                pair.mesh()
            )
        );
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar 
Foam::FSHeatTransferCoefficientModels::multiRegimeBoiling::value
(
    const label& celli
) const
{
    if (FFPairPtr_ == nullptr)
    {
        HashTable<const FFPair*> FFPairs(pair_.mesh().lookupClass<FFPair>());
        FFPairPtr_ = FFPairs[FFPairs.toc()[0]];
    }
    if (dmdtWPtr_ == nullptr)
    {
        dmdtWPtr_ = 
            &
            (
                pair_.mesh().lookupObjectRef<volScalarField>
                (
                    "dmdtW."+FFPairPtr_->name()
                )
            );
    }

    const scalar& Twi(Tw_[celli]);
    const scalar& Tfi(Tf_[celli]);
    const scalar& Tsati(FFPairPtr_->iT()[celli]);

    scalar htc2pFCi(htcFCPtr_->value(celli)*FPtr_->value(celli));
    scalar qFCi(htc2pFCi*(Twi-Tfi));

    if (TONBPtr_.valid())
    {
        scalar TONBi(TONBPtr_->value(celli, htc2pFCi));
        if (Twi > TONBi)
        {
            //- Calc hPB at boiling onset, i.e. the hPB when the wall temperature
            //  is TONBi. This is done by caching Twi, setting it to TONB, 
            //  using the htcPB model to calc the hPB with Tw=TONB, then 
            //  re-setting the Twi. Forgive the const_cast, Stroustroup...
            scalar Twi0(Twi);
            const_cast<scalar&>(Twi) = TONBi;
            scalar hPBONBi(htcPBPtr_->value(celli));
            const_cast<scalar&>(Twi) = Twi0;
            scalar qBIi(hPBONBi*(TONBi-Tsati));

            //-
            scalar hPBi(htcPBPtr_->value(celli));
            scalar qPBi(hPBi*(Twi-Tsati));

            scalar qNBi = pow(pow(qFCi, 3) + pow(qPBi-qBIi, 3), 1.0/3.0);

            /*
            Info<< celli << " Tw=" << Twi << " Tf=" << Tfi 
                << " Ts=" << Tsati << " TONB=" << TONBi
                << " hFC =" << htc2pFCi
                << " hPB =" << hPBi
                << " hPBONB =" << hPBONBi
                << " qNB= " << qNBi 
                << " qPB= " << qPBi 
                << " qBI= " << qBIi; 
            */

            if (Tfi < Tsati)
            {
                //- Sub-cooled boiling
                scalar f(1.0);
                if (SCBFPtr_.valid())
                {
                    f = SCBFPtr_->value(celli, qNBi);
                }
                scalar qSCDmdti(f*(qNBi-qFCi));

                //- If fluid1 is liquid and 2 is vapour then L > 0 and this
                //  sub-cooled boiling term is also > 0. If fluid2 is liquid
                //  and fluid1 is vapour then L < 0 and everything still 
                //  works out in terms of the dmdtW sign, as I recall that
                //  it is positive for phase changes from fluid1 to fluid2 and
                //  negative vice-versa
                (*dmdtWPtr_)[celli] = 
                    pair_.structureRef().iAact()[celli]*qSCDmdti/
                    mag(FFPairPtr_->L()[celli]);
            }

            scalar dT(Twi-Tfi);
            dT = (dT >= 0.0) ? max(dT, 1e-3) : min(-dT, -1e-3);
            return qNBi/dT;
        }
        else
        {
            return htc2pFCi;
        }
    }
    else //- No TONB model ever specified, i.e. boiling starts only after Tf >
         //  Tsat, i.e. no sub-cooled boiling
    {
        if ((Tfi > Tsati) and (Twi > Tsati))
        {
            scalar hPBi(htcPBPtr_->value(celli));
            scalar qPBi(hPBi*(Twi-Tsati));
            scalar qNBi = pow(pow(qFCi, 3) + pow(qPBi, 3), 1.0/3.0);
            scalar dT(Twi-Tfi);
            dT = (dT >= 0.0) ? max(dT, 1e-3) : min(-dT, -1e-3);
            return qNBi/dT;
        }
        else
        {
            return htc2pFCi;
        }
    }
}

// ************************************************************************* //
