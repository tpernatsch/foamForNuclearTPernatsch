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
    Foam::FSHeatTransferCoefficientModels::superposition

Description
    Heat transfer coefficient for that covers many regimes, from single phase
    convection, two-phase convection, sub-cooled and nucleate boiling, to
    film boiling (the latter, hopefully).

Dictionary entries

Author
    Stefan Radman (sradman@pm.me / stefan.radman@epfl.ch)

SourceFiles
    multiRegimeBoiling.C

\*---------------------------------------------------------------------------*/

#ifndef multiRegimeBoilingFSHeatTransferCoefficient_H
#define multiRegimeBoilingFSHeatTransferCoefficient_H

#include "FSHeatTransferCoefficientModel.H"
#include "subCooledBoilingFractionModel.H"
#include "flowEnhancementFactorModel.H"
#include "suppressionFactorModel.H"
#include "TONBModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class FFPair;

namespace FSHeatTransferCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class superposition Declaration
\*---------------------------------------------------------------------------*/

class multiRegimeBoiling
:
    public FSHeatTransferCoefficientModel
{
protected:

    //- Ref to wall temperature
    const volScalarField& Tw_;

    //- Ref to fluid temperature
    const volScalarField& Tf_;

    //- Ptr to FFPair (breaks encapsulation but so does the objectRegistry)
    mutable const FFPair* FFPairPtr_;

    //- Forced convection heat transfer coefficient
    autoPtr<FSHeatTransferCoefficientModel> htcFCPtr_;

    //- Pool boiling heat transfer coefficient
    autoPtr<FSHeatTransferCoefficientModel> htcPBPtr_;

    //- Film condensation heat transfer coefficient
    autoPtr<FSHeatTransferCoefficientModel> htcCndPtr_;

    //- Flow enhcement factor model: increases the convective heat trasfer
    //  coefficient to account for bubble turbulence, increase in flow velocity
    //  due to the vapour generation
    autoPtr<flowEnhancementFactorModel> FPtr_;

    //- Suppression factor model: "turns-off" the pool boiling component
    //  as the boiling transition to film boiling
    autoPtr<suppressionFactorModel> SPtr_;

    //- Fraction of sub-cooled heat transfer that results in vapour generation
    autoPtr<subCooledBoilingFractionModel> SCBFPtr_;

    //- Temperature of the onset of nucleate boiling model
    autoPtr<TONBModel> TONBPtr_;

    //- Ptr to mass transfer term at the wall (of the phaseChange class) so
    //  that it can be set here to account e.g. for sub-cooled boiling.
    //  Please, never forget that encapsulation leads to total madness if 
    //  followed religiously. Still I encourage people to find more
    //  elegant and simpler solutions, if any
    mutable volScalarField* dmdtWPtr_;

    //- Superposition exponent
    scalar exp_;
    scalar oneByExp_;

    //- Superposition is performed on heat fluxes rather than heat transfer 
    //  coefficients
    bool qMode_;

    //- Cache (yikes!)
    mutable scalar htc2pFCi_;
    mutable scalar TONBi_;
    mutable scalar htcNBi_;

//- Protected methods

    //- Set the subcooled boiling mass transfer term for celli 
    void setDmdtW
    (
        const label& celli, 
        const scalar& qNBi, 
        const scalar& qFCi
    ) const;

    //- Return the nucleate boiling heat flux in celli for a given value of the
    //  wall temperature (will be used to numerically calculate the derivative 
    //  of qNB, which is required to find the critical heat flux temperature 
    //  given the critical flux value)
    scalar qNB
    (
        const label& celli, 
        const scalar& Twi1
    ) const; 

public:

    TypeName("multiRegimeBoiling");

    //- Constructors

        //- Construct for a fluid-structure pair
        multiRegimeBoiling
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~multiRegimeBoiling(){}


    //- Member Functions

        //- Return htc for cell celli
        virtual scalar value(const label& celli) const;

        
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FSHeatTransferCoefficientModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //