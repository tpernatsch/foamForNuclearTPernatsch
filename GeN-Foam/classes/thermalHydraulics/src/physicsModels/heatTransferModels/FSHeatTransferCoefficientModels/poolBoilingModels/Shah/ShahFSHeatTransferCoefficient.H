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
    Foam::FSHeatTransferCoefficientModels::Shah

Description
    Drag coefficient in the form fd = coeff*Re^exp with 
    Kd = 0.5*fd*alpha*rho*magU/Dh

Dictionary entries (mandatory)
    coeff   0.316;
    exp     -0.25;

SourceFiles
    Shah.C

\*---------------------------------------------------------------------------*/

#ifndef ShahFSHeatTransferCoefficient_H
#define ShahFSHeatTransferCoefficient_H

#include "FSHeatTransferCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace FSHeatTransferCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class Shah Declaration
\*---------------------------------------------------------------------------*/

class Shah
:
    public FSHeatTransferCoefficientModel
{
private:

    bool useExplicitHeatFlux_;

    //- Ptr to structure wall temperature (only used if !useExplicitHeatFlux_)
    const volScalarField* Twall_;
 
    //- Ptr to fluid temperature (only used if !useExplicitHeatFlux_)
    const volScalarField* Tf_;

    //- Ptr to structure heat flux (only used if useExplicitHeatFlux_)
    const volScalarField* q_;

    //- Ref to pressure
    const volScalarField& p_;

    //- Critical pressure in Pa
    scalar pCrit_;

    //- Model coefficients
    scalar C0_;
    scalar C1_;
    scalar m0_;
    scalar m1_;
    scalar n_;
    scalar exp_;

    //- Interpolation limits (I've chones them arbitrarily but and abrupt
    //  change in the model coefficients is for sure not helpful). The 
    //  abrupt threhsold provided by Shah is at a reduced pressure (i.e.
    //  p_/pCrit_) of 1e-3
    scalar pR0_;
    scalar deltaPR_;

public:

    TypeName("Shah");

    //- Constructors

        //- Construct for a fluid-structure pair
        Shah
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~Shah(){}


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