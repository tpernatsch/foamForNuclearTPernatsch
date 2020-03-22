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

#include "Schor.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidInterfacialAreaModels
{
    defineTypeNameAndDebug(Schor, 0);
    addToRunTimeSelectionTable
    (
        fluidInterfacialAreaModel, 
        Schor, 
        fluidInterfacialAreaModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidInterfacialAreaModels::Schor::Schor
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const fluid& dispersed,
    const fluid& continuous
)
:
    fluidInterfacialAreaModel
    (
        objReg,
        dict,
        dispersed,
        continuous
    ),
    alpha_
    (
        (
            fvc::domainIntegrate(dispersed_.thermo().hc()).value() 
        >   fvc::domainIntegrate(continuous_.thermo().hc()).value()
        ) ?
        //(dispersed_.name() == this->get<word>("alphaName")) ?
        dispersed_
        :
        continuous_
    ),
    D_("D", dimLength, dict),
    P_("P", dimLength, dict),
    PD_(P_/D_),
    A_
    (
        2.0*Foam::sqrt(3.0)*Foam::sqr(PD_)
    ),
    PI_(3.1415927),
    alpha1_("", dimless, 0.55),
    alpha2_("", dimless, 0.65),
    iA1_
    (
        (4.0/D_)*Foam::sqrt
        (
            PI_*alpha1_/(3.0*(A_-PI_))
        )
    ),
    iA2_
    (
        (4.0/D_)*Foam::sqrt(PI_)
        /
        (
            3.0*(A_-PI_)
        )*
        Foam::sqrt
        (
            (1.0-alpha2_)*A_ + PI_*alpha2_
        )
    )
{
    /*
    word alphaName(this->get<word>("alphaName"));
    if 
    (
        alphaName != dispersed_.name() 
        and 
        alphaName != continuous_.name()
    )
    {
        FatalErrorInFunction
            << "Fluid " << alphaName << " specified for the Schorr model does "
            << "not match any existing fluid : " << endl
            << "- " << dispersed_.name() << endl
            << "- " << continuous_.name() << exit(FatalError);
    }
    */
    if 
    (
        fvc::domainIntegrate(dispersed_.thermo().hc()).value() 
    ==  fvc::domainIntegrate(continuous_.thermo().hc()).value()
    )
    {
        FatalErrorInFunction
            << "Fluids must have different enthalpies of formation "
            << "(thermophysicalProperties.Hf) in order to determine which "
            << "fluid is the vapour and which is the liquid" 
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidInterfacialAreaModels::Schor::~Schor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::fluidInterfacialAreaModels::Schor::iA() const
{
    volScalarField alphaSum(dispersed_+continuous_);
    
    volScalarField alpha(alpha_/alphaSum);

    volScalarField weight((alpha2_-alpha)/(alpha2_-alpha1_));

    tmp<volScalarField> tiA
    (
        new volScalarField
        (
            neg(alpha-alpha1_)*
            (
                (4.0/D_)*Foam::sqrt
                (
                    PI_*alpha/(3.0*(A_-PI_))
                )
            ) 

            +
            
            pos0(alpha-alpha2_)*
            (
                (4.0/D_)*Foam::sqrt(PI_)
                /
                (
                    3.0*(A_-PI_)
                )*
                Foam::sqrt
                (
                    (1.0-alpha)*A_ + PI_*alpha
                )*
                min
                (
                    (1-alpha)/(0.043),
                    1.0
                )
            )

            + 

            pos0(alpha-alpha1_)*neg(alpha-alpha2_)*
            (
                weight*iA1_ + (1.0-weight)*iA2_
            )
        )
    );
    
    return tiA;
}


// ************************************************************************* //
