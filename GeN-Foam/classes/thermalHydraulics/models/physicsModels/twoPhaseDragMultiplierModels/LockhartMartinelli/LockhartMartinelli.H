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

Class
    Foam::twoPhaseDragMultiplierModels::LockhartMartinelli

Description
    Model by Lockart and Martinelli for the multiplier:

    phi2 = 1 + C/X + 1/X^2

    with C being a regime dependent parameter that varies for different 
    combinations of laminar and turbulent flow of each of the two phases. It
    defaults to universally used value of 20, but can be user-specified.
    X is the Lockhart Martinelli parameter defined as:

    X = (mu_I/mu_O)^exp * (flowQuality_I/flowQuality_O)^(1-exp) * rho_O/rho_I

    with I being the phase whose single-phase drag is used as basis to compute
    the two-phase multiplier and O being the other phase. The exponent exp
    defaults to 0.2 but can be user-specified. The choise of the exponent
    should reflect the general behaviour of the single phase friction factor
    Reynolds (e.g., if you are using Blasius, i.e. f ~ 0.316/Re^(0.25), then
    exp = 0.25). If you are not using Blasius, or something with a clear
    dependence of f on Re through a single exponent, just use the generally
    accepted value of 0.2

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    LockhartMartinelli.C

DictEntryExample ("" is the placeholder for the variable to be inserted)

    twoPhaseDragMultiplierModel
    {
        type            LockhartMartinelli;
        multiplierFluid liquid;
        C               20;         //- Not mandatory, defaults to 20
        exp             0.2;        //- Not mandatory, defaults to 0.2
    }


\*---------------------------------------------------------------------------*/

#ifndef LockhartMartinelli_H
#define LockhartMartinelli_H

#include "twoPhaseDragMultiplierModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace twoPhaseDragMultiplierModels
{

/*---------------------------------------------------------------------------*\
                           Class LockhartMartinelli Declaration
\*---------------------------------------------------------------------------*/

class LockhartMartinelli
:
    public twoPhaseDragMultiplierModel
{
private:

    //- C in phi2 = 1 + C/X + 1/X^2
    scalar C_;

public:

    //- Runtime type information
    TypeName("LockhartMartinelli");

    //- Constructors

        LockhartMartinelli
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~LockhartMartinelli();

    //- Member Functions

        //- Calculate phi2_;
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseDragMultiplierModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //