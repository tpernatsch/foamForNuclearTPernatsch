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
    Foam::ZyOverstrainRIA

Description
    Overstrain criterion to asses clad failure during the PCMI phase of a RIA transient.

    If the hoop strain is superior to a critical value, the clad is considered failed.
    Critical value of hoop strain is a function of Temperature, strain rate, hydrogen content and fast neutron fluence (E > 1 MeV)

    From: https://inis.iaea.org/collection/NCLCollectionStore/_Public/42/022/42022542.pdf

Usage
    In solverDict file:
    \verbatim
    materials
    {
        cladding
        {
            material zircaloy;

            ...

            failureModel         ZyOverstrainRIA;
            CH       400;
            phi       1; //Fast neutron fluence in 1e25 n/m2
            stopIfFailed         false;
        }

        ...
    }
    \endverbatim

SourceFiles
    ZyOverstrainBISON.C

\mainauthor
    Matthieu Reymond - EPFL

\contribution
    A. Scolaro - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef ZyOverstrainRIA_H
#define ZyOverstrainRIA_H

#include "failureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ZyOverstrainRIA Declaration
\*---------------------------------------------------------------------------*/

class ZyOverstrainRIA
:
    public failureModel
{
    // Private data
    volScalarField hoopStrainLimit_;

    volScalarField epsDotYY_;

    volScalarField K1_;

    volScalarField CFP_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ZyOverstrainRIA(const ZyOverstrainRIA&);

        //- Disallow default bitwise assignment
        void operator=(const ZyOverstrainRIA&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("ZyOverstrainRIA");

    //Concentration Hydrogen (ppm)
    scalar CH_Value_;

    //Fast neutron fluence
    scalar phiValue_;

    //Crack length (m) for stress intensity facotr (K1) field
    scalar crackLengthValue_;


    // Declare run-time constructor selection table
    // Constructors

        //- Construct from dictionary
        ZyOverstrainRIA
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Selectors


    //- Destructor
    virtual ~ZyOverstrainRIA();


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
