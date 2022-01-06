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
    Foam::FSHeatTransferCoefficientModels::multiRegimeBoilingVapourTRACE

Description
    Heat transfer coefficient for the transition Boiling Regime and the 
    POST CHF Regimes. It uses the correlation from TRACE.

Dictionary entries

Author
    Stefan Radman (sradman@pm.me / stefan.radman@epfl.ch)
    Gauthier Lazare (gauthier.lazare@gmail.com)

SourceFiles
    multiRegimeBoilingVapourTRACEFSHeatTransferCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef multiRegimeBoilingVapourTRACEFSHeatTransferCoefficient_H
#define multiRegimeBoilingVapourTRACEFSHeatTransferCoefficient_H

#include "FSHeatTransferCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class FFPair;

namespace FSHeatTransferCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class multiRegimeBoilingVapourTRACE Declaration
\*---------------------------------------------------------------------------*/

class multiRegimeBoilingVapourTRACE
:
    public FSHeatTransferCoefficientModel
{
private:
    //- Ref to pressure 
    const volScalarField& p_;

    //- Ref to wall temperature
    const volScalarField& Tw_;

    // - Boolean values to know if we are after or before CHF 
    const volScalarField& CHFBool_;

    // - Weighted fraction for the transition boiling
    const volScalarField& wfTB_;

    //- Temperature of the Leidenfrost model 
    const volScalarField& TLF_;

    //- Ptr to FFPair (breaking encapsulation as usual)
    mutable const FFPair* FFPairPtr_;

    //- Annular Flow heat transfer coefficient
    autoPtr<FSHeatTransferCoefficientModel> htcAFPtr_;


public:

    TypeName("multiRegimeBoilingVapourTRACE");

    //- Constructors

        //- Construct for a fluid-structure pair
        multiRegimeBoilingVapourTRACE
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~multiRegimeBoilingVapourTRACE(){}


    //- Member Functions

        //- Drag coefficient
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FSHeatTransferCoefficientModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //