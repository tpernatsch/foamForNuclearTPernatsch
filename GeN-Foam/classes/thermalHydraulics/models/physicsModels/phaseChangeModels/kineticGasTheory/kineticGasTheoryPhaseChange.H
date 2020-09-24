/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::phaseChangeModels::kineticGasTheoryPhaseChange

Description
    Phase change model based on the Kinetic Gas Theory. This is a mix of the
    approach found in the NATOF-2D code (see report MIT-EL 81-030 by R. G.
    Zielinksi and M. S. Kazimi, 09/1981) and the SABENA code (see H. Ninokata,
    T. Okano, "SABENA: Subassembly Boiling Evolution Numerical Analysis", 
    Nuclear Engineering and Design vol. 120, pp. 349-367, 1990). It does not
    work well with sodium due to its explicitness, but I did not have the time
    to thoroughly test it. Do not use with constant saturation models, it won't
    work.

Author
    Stefan Radman (stefanradman92@gmail.com / stefan.radman@epfl.ch)

SourceFiles
    kineticGasTheoryPhaseChange.C

\*---------------------------------------------------------------------------*/

#ifndef kineticGasTheoryPhaseChange_H
#define kineticGasTheoryPhaseChange_H

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeModels
{

/*---------------------------------------------------------------------------*\
                         Class kineticGasTheoryPhaseChange Declaration
\*---------------------------------------------------------------------------*/

class kineticGasTheoryPhaseChange
:
    public phaseChangeModel
{
protected:

    //- Sigma. This represents the fraction of molecules from either phase 
    //  that actually move through the interface to join the other phase
    //  rather than "bouncing back" thus not resulting in a change of state
    const dimensionedScalar sigma_;

    //- Coeff of the interfacial mass flux computed once as 
    //      (2*sigma/(2-sigma))*sqrt(M/2*pi*R)
    //  with M being the molar mass in kg/mol and R the universal gas constant
    const dimensionedScalar coeff_;

    //- Phase fraction that establishes wether the mass transfer is governed
    //  by the liquid temperature or vapour temperature
    //const dimensionedScalar alphaDryout_;

    //- Refs to liquid, vapour temperature and phase fraction
    const fluid& liquid_;
    const fluid& vapour_;

public:

    TypeName("kineticGasTheory");

    //- Constructors

        //- Construct from a dictionary
        kineticGasTheoryPhaseChange
        (
            const dictionary& dict,
            const pimpleControl& pimple, 
            const fluid& fluid1,
            const fluid& fluid2,
            const volScalarField& p,
            const volScalarFieldPtrTable& htcs,
            volScalarField& dmdt,
            volScalarField& iT,
            volScalarField& iA
        );

        //- Destructor
        virtual ~kineticGasTheoryPhaseChange(){}

    //- Member Functions

        //- Return mass transfer term
        virtual void correctMassTransfer();
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseChangeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
