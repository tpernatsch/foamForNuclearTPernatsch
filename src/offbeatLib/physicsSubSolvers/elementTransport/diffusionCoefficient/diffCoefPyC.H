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
    Foam::diffCoefPyC

Description
    Diffusion coefficient for PyC material from Parfume manual and INL report.


SourceFiles
    diffCoefPyC.C

\mainauthor
    F. Xiang - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour) & XJTU (Xi'an JiaoTong
    University, China, Nuclear Thermal-hydraulic Laboratory)

\contribution
    A. Scolaro - EPFL
    E. Brunetto - EPFL

\date
    June 2023

\*---------------------------------------------------------------------------*/

#ifndef diffCoefPyC_H
#define diffCoefPyC_H

#include "diffCoefModel.H"
#include "primitiveFields.H"
#include "fvMesh.H"
#include "volFields.H"

#include "userParameters.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class diffCoefPyC Declaration
\*---------------------------------------------------------------------------*/

class diffCoefPyC
:
    public diffCoefModel
{
    // Private data
    
    // Private Member Functions

    //- Disallow default bitwise copy construct
    diffCoefPyC(const diffCoefPyC&);

    //- Disallow default bitwise assignment
    void operator=(const diffCoefPyC&);

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("PyCArrhenius");

    // Constructors

    //- Construct from dictionary
    diffCoefPyC
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word defaultModel
    );

    //- Destructor
    virtual ~diffCoefPyC();


    // Member Functions

    virtual void updateCsCoef(volScalarField& sf, const volScalarField& T, const labelList& addr);
    virtual void updateKrCoef(volScalarField& sf, const volScalarField& T, const labelList& addr);
    virtual void updateAgCoef(volScalarField& sf, const volScalarField& T, const labelList& addr);
    virtual void updateSrCoef(volScalarField& sf, const volScalarField& T, const labelList& addr);

    //- Update diffusion Coefficient
    virtual void updateCoef(volScalarField& sf, const volScalarField& T, const labelList& addr, const word FpName);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
