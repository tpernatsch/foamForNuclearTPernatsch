/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "porouskEpsilon.H"
#include "bound.H"
#include "byZoneCorrelationPorousMedium.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


template<class BasicTurbulenceModel>
void porouskEpsilon<BasicTurbulenceModel>::correctNut()
{
    volScalarField clearYesNo = 1.0 - mag(porousMedium_.kepsilonConvergenceRate())/(mag(porousMedium_.kepsilonConvergenceRate())+dimensionedScalar("", dimensionSet(0,0,-1,0,0,0,0), SMALL) );

     volScalarField nuStab = pow(mag(this->U_),2)*porousMedium_.hydraulicDiameterStructure() / 100.0; //100 e' un Reynolds laminare 

     this->nut_ = (this->Cmu_*sqr(this->k_)/this->epsilon_)*clearYesNo + (nuStab)*(1-clearYesNo);
     this->nut_.correctBoundaryConditions();
     //fv::options::New(this->mesh_).correct(this->nut_);
 
     BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> porouskEpsilon<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> porouskEpsilon<BasicTurbulenceModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
porouskEpsilon<BasicTurbulenceModel>::porouskEpsilon
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            -0.33
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    porousMedium_(this->mesh_.objectRegistry::lookupObject<correlationPorousMedium>("porousMediumProperties"))

{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool porouskEpsilon<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void porouskEpsilon<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;

    eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G(this->GName(), nut*(dev(twoSymm(tgradU())) && tgradU()));
    tgradU.clear();

    // Update epsilon and G at the wall
    this->epsilon_.boundaryField().updateCoeffs();

    volScalarField clearYesNo = 1.0 - mag(porousMedium_.kepsilonConvergenceRate())/(mag(porousMedium_.kepsilonConvergenceRate())+dimensionedScalar("", dimensionSet(0,0,-1,0,0,0,0), SMALL) );

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, this->epsilon_)
      + fvm::div(alphaRhoPhi, this->epsilon_)
      - fvm::laplacian(alpha*rho*this->DepsilonEff(), this->epsilon_)
     ==
        this->C1_*alpha*rho*G*this->epsilon_/this->k_ * clearYesNo
      - fvm::SuSp(((2.0/3.0)*this->C1_ + this->C3_)*alpha*rho*divU * clearYesNo, this->epsilon_)
      - fvm::Sp(this->C2_*alpha*rho*this->epsilon_/this->k_ * clearYesNo, this->epsilon_)
      + epsilonSource()
      - fvm::Sp(alpha*rho* porousMedium_.kepsilonConvergenceRate(), this->epsilon_)
      + alpha*rho* porousMedium_.kepsilonConvergenceRate() * porousMedium_.equilibriumEpsilon()
    );
    epsEqn.ref().relax();
    epsEqn().boundaryManipulate(this->epsilon_.boundaryField());
    solve(epsEqn);
    bound(this->epsilon_, this->epsilonMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, this->k_)
      + fvm::div(alphaRhoPhi, this->k_)
      - fvm::laplacian(alpha*rho*this->DkEff(), this->k_)
     ==
        alpha*rho*G * clearYesNo
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU * clearYesNo, this->k_)
      - fvm::Sp(alpha*rho*this->epsilon_/this->k_ * clearYesNo, this->k_)
      - fvm::Sp(alpha*rho* porousMedium_.kepsilonConvergenceRate(), this->k_)
      + alpha*rho* porousMedium_.kepsilonConvergenceRate() * porousMedium_.equilibriumK()
      + kSource()
    );
    kEqn.ref().relax();
    solve(kEqn);
    bound(this->k_, this->kMin_);

    this->correctNut();


}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam




// ************************************************************************* //
