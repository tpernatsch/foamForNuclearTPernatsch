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

#include "DarcyReynoldsPower.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(DarcyReynoldsPower, 0);
    addToRunTimeSelectionTable(dragModel, DarcyReynoldsPower, FSDragModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::DarcyReynoldsPower::DarcyReynoldsPower
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FSPair& FSPair
)
:
    dragModel
    (
        objReg,
        dict,
        FSPair
    ),
    coeff_(dict.get<scalar>("coeff")),
    exp_(dict.get<scalar>("exp"))
{
    /*
    const entry* eCoeff(this->findEntry("coeff"));
    const entry* eExp(this->findEntry("exp"));
    if (eCoeff->stream().size() == 1 and eExp->stream().size() == 1)
    {
        scalar coeff(this->get<scalar>("coeff"));
        scalar exp(this->get<scalar>("exp"));
        coeff_ = vector(coeff, coeff, coeff);
        exp_ = vector(exp, exp, exp);
        isotropic_ = true;
    }
    else if (eCoeff->stream().size() != 1 and eExp->stream().size() != 1)
    {
        coeff_ = this->get<vector>("coeff");
        exp_ = this->get<vector>("exp");
        isotropic_ = false;
    }
    else
    {
        FatalErrorInFunction
            << "coeff and exp must both be either scalars or vectors!" << endl
            << exit(FatalError);
    }
    */
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::dragModels::DarcyReynoldsPower::fd(const scalar& Re) const
{
    return coeff_*pow(Re, exp_);
}

void Foam::dragModels::DarcyReynoldsPower::correctKd(volTensorField& Kd) const
{   
    #include "calcKdFromFd.H"
    /*
    const volScalarField& alpha(FSPair_->fluidRef()); 
    const volScalarField& rho(FSPair_->fluidRef().thermo().rho());
    const volScalarField& magU(FSPair_->fluidRef().magU());
    
    if (isotropic_)
    {
        const volScalarField& Dh(FSPair_->structure().Dh());
        
        const volScalarField& Re(FSPair_->Re());

        forAll(cellList_, i)
        {
            label celli(cellList_[i]);
            scalar value
            (
                0.5*alpha[celli]*rho[celli]*magU[celli]*coeff_[0]*
                pow(Re[celli], exp_[0])/Dh[celli]
            );
            tensor& Kdi(Kd[celli]);
            Kdi[0] = value;
            Kdi[4] = value;
            Kdi[8] = value;
        }
    }
    else
    {
        const volVectorField& lDh(FSPair_->structure().lDh());
        const volVectorField& lRe(FSPair_->lRe());

        //- lRe should be limited
        forAll(cellList_, i)
        {
            label celli(cellList_[i]);
            const vector& lRei(lRe[celli]);
            const vector& lDhi(lDh[celli]);
            tensor& Kdi(Kd[celli]);
            forAll(coeff_, j)
            {
                Kdi[j*4] = 
                    0.5*alpha[celli]*rho[celli]*magU[celli]*coeff_[j]*
                    pow(lRei[j], exp_[j])/lDhi[j];
            }
        }
    }
    */
}


// ************************************************************************* //
