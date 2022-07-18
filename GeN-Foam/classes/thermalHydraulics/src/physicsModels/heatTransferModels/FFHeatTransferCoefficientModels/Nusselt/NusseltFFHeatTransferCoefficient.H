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
    Foam::FFHeatTransferCoefficientModels::Nusselt

Description
    Heat transfer coefficient for forced convective flows that is computed
    starting from the Nusselt, which is turn obtained via correlation in the
    form 

    Nu = A_+B_*(Re^C_)*(Pr^D_)

    with Re, Pr being the Reynolds (superficial in two-phase) and Prandlt
    of the fluid, and A, B, C, D selectable consts

Author
    Stefan Radman (sradman@pm.me / stefan.radman@epfl.ch)

Dictionary entries example (mandatory)
    const   4;
    coeff   0.01;
    expRe   0.8;
    expPr   0.8;

SourceFiles
    NusseltFFHeatTransferCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef NusseltFFHeatTransferCoefficient_H
#define NusseltFFHeatTransferCoefficient_H

#include "FFHeatTransferCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class FFPair;

namespace FFHeatTransferCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class Nusselt Declaration
\*---------------------------------------------------------------------------*/

class Nusselt
:
    public FFHeatTransferCoefficientModel
{
private:

    //- Ref to Reynolds of pair
    const volScalarField& Re_;

    //- Bulk fluid thermal conductivity
    const volScalarField& kappa_;

    //- Ref to Prandtl of bulk fluid
    const volScalarField& Pr_;

    //- Ref to characteristic dimension of bulk fluid
    const volScalarField& Dh_;

    //- Coefficients in Nu = A_ + B_*Re^(C_)*Pr^(D_)
    const scalar A_;
    const scalar B_;
    const scalar C_;
    const scalar D_;

    //- If C_ == D_, no need to do pow(Re, C_)*pow(Pr, D_), you can just do
    //  pow(Re*Pr, C_) which is faster. Re*Pr = Pe (Peclet)
    const bool usePeclet_;

public:

    TypeName("NusseltReynoldsPrandtlPower");

    //- Constructors

        //- Construct for a fluid-structure pair
        Nusselt
        (
            const FFPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~Nusselt(){}


    //- Member Functions

        //- Drag coefficient
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FFHeatTransferCoefficientModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
