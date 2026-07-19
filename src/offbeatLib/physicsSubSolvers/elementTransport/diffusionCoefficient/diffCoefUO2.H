/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::diffCoefUO2

Description
    Built-in diffusion coefficient correlations for fission products in a UO2
    kernel.  Dispatches internally per species (Cs, Ag, Kr, Sr).

    From Parfume manual and INL report.

SourceFiles
    diffCoefUO2.C

\*---------------------------------------------------------------------------*/

#ifndef diffCoefUO2_H
#define diffCoefUO2_H

#include "diffCoefModel.H"
#include "primitiveFields.H"
#include "fvMesh.H"
#include "volFields.H"

#include "userParameters.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class diffCoefUO2 Declaration
\*---------------------------------------------------------------------------*/

class diffCoefUO2
:
    public diffCoefModel
{
    // Private Member Functions

    diffCoefUO2(const diffCoefUO2&);
    void operator=(const diffCoefUO2&);

protected:

public:

    //- Runtime type information
    TypeName("UO2Arrhenius");

    // Constructors

    //- Construct from dictionary
    diffCoefUO2
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word defaultModel
    );

    //- Destructor
    virtual ~diffCoefUO2();


    // Member Functions

    virtual void updateCsCoef(volScalarField& sf, const volScalarField& T, const labelList& addr);
    virtual void updateKrCoef(volScalarField& sf, const volScalarField& T, const labelList& addr);
    virtual void updateAgCoef(volScalarField& sf, const volScalarField& T, const labelList& addr);
    virtual void updateSrCoef(volScalarField& sf, const volScalarField& T, const labelList& addr);

    //- Update diffusion coefficient
    virtual void updateCoef(volScalarField& sf, const volScalarField& T, const labelList& addr, const word FpName);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
