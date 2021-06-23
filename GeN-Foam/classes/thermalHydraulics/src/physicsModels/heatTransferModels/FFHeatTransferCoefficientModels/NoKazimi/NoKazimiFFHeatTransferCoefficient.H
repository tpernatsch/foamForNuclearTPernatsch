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
    Foam::FFHeatTransferCoefficientModels::NoKazimi

Description
    NoKazimi model for the liquid side of a liquid-vapour heat transfer
    coefficient. Can be only used if a phaseChangeModel is specified as this
    models relies on latent heat.

Author
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    NoKazimi.C

\*---------------------------------------------------------------------------*/

#ifndef NoKazimiFFHeatTransferCoefficient_H
#define NoKazimiFFHeatTransferCoefficient_H

#include "FFHeatTransferCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class FFPair;

namespace FFHeatTransferCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class NoKazimi Declaration
\*---------------------------------------------------------------------------*/

class NoKazimi
:
    public FFHeatTransferCoefficientModel
{
private:

    const volScalarField& alpha_;
    
    const volScalarField& rho_;

    const volScalarField& kappa_;

    const volScalarField& Cp_;

    const volScalarField& Pr_;
    
    const volScalarField& magU_;

    //- Ref to pressure
    const volScalarField& p_;

    //- Ref to fluid-fluid interfacial area
    const volScalarField& iA_;

    //- Ref to fluid-fluid interfacial temperature
    const volScalarField& iT_;

    //- A pre-computed coefficients
    scalar coeff_;

public:

    TypeName("NoKazimi");

    //- Constructors

        //- Construct for a fluid-structure pair
        NoKazimi
        (
            const FFPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~NoKazimi(){}


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