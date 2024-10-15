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
    Foam::Zy_CSED_RIA_EPRI

Description
    Clad failure criterion during the PCMI phase of a RIA transient based on the strain energy density (SED).
    The complete description of this model is available at : https://www.epri.com/research/products/1021036
    If the average SED on an axial slice goes over a critical value, the clad is considered to have failed.
    The SED field corresponds to the strain energy density.
    The CFP field corresponds to the SED field over the critical value. Thus if CFP > 1, the clad is considered failed.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        cladding
        {
            material zircaloy;

            ...

            failureModel         ZyCSED;
            CH                   400; //hydrogen content in wpmm. The critical value is dependant on the hydrogen content.
            stopIfFailed         false;
        }

        ...
    }
    \endverbatim

SourceFiles
    Zy_CSED_RIA_EPRI.C

\mainauthor
    M. Reymond - EPFL

\contribution
    A. Scolaro - EPFL

\date 
    February 2023

\*---------------------------------------------------------------------------*/

#ifndef Zy_CSED_RIA_EPRI_H
#define Zy_CSED_RIA_EPRI_H

#include "failureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Zy_CSED_ROA_EPRI Declaration
\*---------------------------------------------------------------------------*/

class Zy_CSED_RIA_EPRI
:
    public failureModel
{
    // Private data
    volScalarField  SED_;

    volScalarField CFP_;

    volScalarField K1_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Zy_CSED_RIA_EPRI(const Zy_CSED_RIA_EPRI&);

        //- Disallow default bitwise assignment
        void operator=(const Zy_CSED_RIA_EPRI&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("ZyCSED");

    //Concentration Hydrogen (ppm)
    scalar CH_Value_;

    //Crack length (m) for stress intensity facotr (K1) field
    scalar crackLengthValue_;

    // Declare run-time constructor selection table
    // Constructors

        //- Construct from dictionary
        Zy_CSED_RIA_EPRI
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Selectors


    //- Destructor
    virtual ~Zy_CSED_RIA_EPRI();


    // Member Functions
    
        //- Return true if material is failed
        virtual bool isFailed(const labelList& addr);

        //- Return name of failure criterion
        virtual word criterionName()
        {
            return typeName_();
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
