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

#include "Autruffe.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(Autruffe, 0);
    addToRunTimeSelectionTable(dragModel, Autruffe, FFDragModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::Autruffe::Autruffe
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FFPair& FFPair
)
:
    dragModel
    (
        objReg,
        dict,
        FFPair
    ),
    vapour_
    (
        ( 
            fvc::domainIntegrate(FFPair.fluid1().thermo().hc()).value() >
            fvc::domainIntegrate(FFPair.fluid2().thermo().hc()).value()  
        ) ?
        FFPair.fluid1() : FFPair.fluid2()
    ),
    liquid_
    (
        (FFPair.fluid1().name() == vapour_.name()) ?
        FFPair.fluid2() : FFPair.fluid1()
    )
{
    if 
    (
        fvc::domainIntegrate(FFPair.fluid1().thermo().hc()).value()
    ==  fvc::domainIntegrate(FFPair.fluid2().thermo().hc()).value() 
    )
    {
        FatalErrorInFunction
            << "Fluids must have different enthalpies of formation "
            << "(thermophysicalProperties.Hf) in order to determine which "
            << "fluid is the vapour and which is the liquid" 
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dragModels::Autruffe::correctKd(volTensorField& Kd) const
{
    const volScalarField& DhContinuous(FFPair_->DhContinuous());
    const volScalarField& magUr(FFPair_->magUr());
    const volScalarField& rhoVap(vapour_.rho());

    //- I don't care about mesh boundaries, cell-by-cell is faster
    forAll(mesh_.cells(), i)
    {
        tensor& Kdi(Kd[i]);
        scalar alpha(vapour_[i]/(vapour_[i]+liquid_[i]));
        scalar value
        (
            4.31/(2*DhContinuous[i])*magUr[i]*rhoVap[i]*
            pow
            (
                (1.0-alpha)*
                (1.0+75.0*(1.0-alpha)),
                0.95
            )
        );

        Kdi[0] = value;
        Kdi[4] = value;
        Kdi[8] = value;
    }
}


// ************************************************************************* //
